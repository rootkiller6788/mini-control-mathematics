/**
 * @file example_digital_filter.c
 * @brief Digital Lowpass Filter Design — Bilinear Transform Method
 *
 * L6: Design a Butterworth IIR lowpass filter from analog prototype
 *     using the bilinear (Tustin) transform.
 * L7: Application — audio anti-aliasing filter for 48kHz sampling.
 */
#include "../include/digital_filter.h"
#include "../include/bilinear_transform.h"
#include "../include/z_transform.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    printf("=== Digital Butterworth Lowpass Filter Design ===\n");
    printf("Application: Audio anti-aliasing for 48 kHz ADC\n\n");

    /* Specifications */
    FilterSpec spec;
    spec.type               = FILTER_LOWPASS;
    spec.sampling_freq      = 48000.0;    /* 48 kHz */
    spec.passband_edge      = 18000.0;    /* Pass up to 18 kHz */
    spec.stopband_edge      = 24000.0;    /* Stop by Nyquist (24 kHz) */
    spec.passband_ripple_db = 0.5;        /* 0.5 dB ripple */
    spec.stopband_atten_db  = 60.0;       /* 60 dB attenuation */

    printf("Specifications:\n");
    printf("  Fs      = %.0f Hz\n", spec.sampling_freq);
    printf("  F_pass  = %.0f Hz\n", spec.passband_edge);
    printf("  F_stop  = %.0f Hz\n", spec.stopband_edge);
    printf("  A_pass  = %.1f dB\n", spec.passband_ripple_db);
    printf("  A_stop  = %.1f dB\n", spec.stopband_atten_db);

    /* Estimate required order */
    int N_butter = butterworth_order(&spec);
    int N_cheby  = chebyshev1_order(&spec);
    printf("\nRequired filter order:\n");
    printf("  Butterworth: N = %d\n", N_butter);
    printf("  Chebyshev I: N = %d\n", N_cheby);

    /* Design Butterworth IIR */
    FilterCoeff coeff;
    int ret = iir_filter_design(&spec, IIR_BUTTERWORTH, &coeff);
    if (ret != 0) {
        printf("Filter design failed.\n");
        return 1;
    }

    printf("\nDesigned IIR filter:\n");
    printf("  Numerator order M = %d\n", coeff.M);
    printf("  Denominator order N = %d\n", coeff.N);
    printf("  FIR: %s\n", coeff.is_fir ? "yes" : "no (IIR)");

    printf("\nNumerator coefficients b[k]:\n");
    for (int i = 0; i <= coeff.M; i++)
        printf("  b[%d] = % .8f\n", i, coeff.b[i]);

    printf("\nDenominator coefficients a[k]:\n");
    for (int i = 0; i <= coeff.N; i++)
        printf("  a[%d] = % .8f\n", i, coeff.a[i]);

    /* Frequency response at key points */
    double test_freqs[] = {1000.0, 10000.0, 18000.0, 24000.0};
    int nf = sizeof(test_freqs) / sizeof(test_freqs[0]);

    printf("\nFrequency response:\n");
    printf("  Freq (Hz)    |H| (dB)    Phase (rad)\n");
    printf("  ---------    --------    -----------\n");
    for (int i = 0; i < nf; i++) {
        double omega = 2.0 * M_PI * test_freqs[i] / spec.sampling_freq;
        double mag_db, phase;
        filter_freq_response(&coeff, omega, &mag_db, &phase);
        printf("  %8.0f    %8.2f    %8.4f\n", test_freqs[i], mag_db, phase);
    }

    /* Bilinear transform verification: discretize a simple 1/(s+1) */
    printf("\n=== Bilinear Transform Verification ===\n");
    TransferFunction ct;
    ct.num_order = 0; ct.num_coeff[0] = 1.0;
    ct.den_order = 1; ct.den_coeff[0] = 1.0; ct.den_coeff[1] = 1.0;
    ct.gain = 1.0; ct.sampling_period = 0.0;

    TransferFunction dt;
    double T_verify = 1.0 / 1000.0;  /* 1 ms sampling */
    bilinear_discretize(&ct, T_verify, 0.0, &dt);

    printf("Continuous G(s) = 1/(s+1)\n");
    printf("Discrete G(z) at T=%.3f ms:\n", T_verify * 1e3);
    printf("  b[0] = %.6f, b[1] = %.6f\n", dt.num_coeff[0], dt.num_coeff[1]);
    printf("  a[0] = %.6f, a[1] = %.6f\n", dt.den_coeff[0], dt.den_coeff[1]);

    /* Analytical: G(z) = (T/2)·(1+z⁻¹) / (1 + T/2 + (T/2 - 1)z⁻¹) */
    double T = T_verify;
    double b0_expected = T / 2.0;
    double b1_expected = T / 2.0;
    double a0_expected = 1.0 + T / 2.0;
    double a1_expected = T / 2.0 - 1.0;

    printf("Expected:\n");
    printf("  b[0] = %.6f, b[1] = %.6f\n", b0_expected, b1_expected);
    printf("  a[0] = %.6f, a[1] = %.6f\n", a0_expected, a1_expected);

    filter_coeff_free(&coeff);
    printf("\n=== Design Complete ===\n");
    return 0;
}
