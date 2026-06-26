/**
 * @file transfer_function.c
 * @brief Transfer Function Analysis — Evaluation, Frequency Response, Specs
 *
 * Knowledge points: DC gain, system type, 2nd-order params, frequency response,
 *                   Bode sweep, resonant peak, bandwidth, pole-zero map,
 *                   minimum-phase check, relative degree.
 */

#include "transfer_function.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*----------------------------------------------------------------------
 * L2 — Transfer Function Evaluation
 *----------------------------------------------------------------------*/

double complex tf_evaluate(const TransferFunction *tf, double complex s_z)
{
    if (!tf) { errno = EINVAL; return NAN + I * NAN; }

    /* Evaluate numerator */
    double complex N = tf->num_coeff[tf->num_order];
    for (int k = tf->num_order - 1; k >= 0; k--)
        N = N * s_z + tf->num_coeff[k];

    /* Evaluate denominator */
    double complex D = tf->den_coeff[tf->den_order];
    for (int k = tf->den_order - 1; k >= 0; k--)
        D = D * s_z + tf->den_coeff[k];

    if (cabs(D) < DBL_EPSILON * 1e6) { errno = EDOM; return NAN + I * NAN; }

    return tf->gain * N / D;
}

/*----------------------------------------------------------------------
 * L2 — DC Gain
 *----------------------------------------------------------------------*/

double tf_dc_gain(const TransferFunction *tf)
{
    if (!tf) { errno = EINVAL; return NAN; }

    if (tf->sampling_period <= 0.0) {
        /* Continuous: G(0) */
        double complex val = tf_evaluate(tf, 0.0 + 0.0 * I);
        return creal(val);
    } else {
        /* Discrete: G(1) */
        double complex val = tf_evaluate(tf, 1.0 + 0.0 * I);
        return creal(val);
    }
}

/*----------------------------------------------------------------------
 * L2 — System Type (number of integrators)
 *----------------------------------------------------------------------*/

SystemType tf_system_type(const TransferFunction *tf)
{
    if (!tf) return SYSTEM_TYPE_0;

    int type = 0;

    if (tf->sampling_period <= 0.0) {
        /* Continuous: count leading zero coefficients in denominator */
        /* D(s) = a₀ + a₁·s + ... + aₙ·sⁿ. Poles at s=0 ↔ a₀ = a₁ = ... = a_{k-1} = 0 */
        for (int i = 0; i <= tf->den_order; i++) {
            if (fabs(tf->den_coeff[i]) < DBL_EPSILON * 1e3)
                type++;
            else
                break;
        }
    } else {
        /* Discrete: poles at z=1 → D(1)=0. Count multiplicity by D(z) at z=1 */
        /* D(z) = Σ a_i·z^(-i). D(1) = Σ a_i = 0 means at least type 1. */
        /* For higher multiplicity, check successive derivatives at z=1. */
        double sum = 0.0;
        for (int i = 0; i <= tf->den_order; i++)
            sum += tf->den_coeff[i];
        if (fabs(sum) < 1e-10) {
            type = 1;
            /* Check D'(1): derivative of D(z^{-1}) w.r.t z^{-1} at z=1 = Σ i·a_i */
            double dsum = 0.0;
            for (int i = 1; i <= tf->den_order; i++)
                dsum += i * tf->den_coeff[i];
            if (fabs(dsum) < 1e-10) {
                type = 2;
                /* Check D''(1) */
                double d2sum = 0.0;
                for (int i = 2; i <= tf->den_order; i++)
                    d2sum += i * (i - 1) * tf->den_coeff[i];
                if (fabs(d2sum) < 1e-10) type = 3;
            }
        }
    }

    if (type > 3) type = 3;
    return (SystemType)type;
}

/*----------------------------------------------------------------------
 * L5 — Extract 2nd-Order Parameters ωₙ, ζ
 *----------------------------------------------------------------------*/

int tf_extract_second_order(const TransferFunction *tf,
                             SecondOrderParams *params)
{
    if (!tf || !params) { errno = EINVAL; return -1; }
    params->omega_n = 0.0;
    params->zeta   = 0.0;
    params->is_standard_form = 0;

    if (tf->den_order < 2) return -1;
    if (tf->sampling_period > 0.0) return -1; /* continuous only */

    /* D(s) = d₀ + d₁·s + d₂·s² + ... → normalize to d₂=1 */
    double d2 = tf->den_coeff[2];
    if (fabs(d2) < DBL_MIN) return -1;

    double a0 = tf->den_coeff[0] / d2;
    double a1 = tf->den_coeff[1] / d2;

    /* Standard form: s² + 2ζωₙ·s + ωₙ² = 0 */
    /* → ωₙ² = a₀, 2ζωₙ = a₁ */
    if (a0 <= 0.0) return -1;  /* ωₙ² must be positive */

    params->omega_n = sqrt(a0);
    params->zeta    = a1 / (2.0 * params->omega_n);
    params->is_standard_form = 1;
    return 0;
}

/*----------------------------------------------------------------------
 * L5 — Frequency Response at Single ω
 *----------------------------------------------------------------------*/

int tf_frequency_response(const TransferFunction *tf,
                           double omega,
                           FrequencyResponse *resp)
{
    if (!tf || !resp || omega < 0.0) { errno = EINVAL; return -1; }

    double complex jw = 0.0 + omega * I;
    double complex G = tf_evaluate(tf, jw);

    resp->omega       = omega;
    resp->real_part   = creal(G);
    resp->imag_part   = cimag(G);

    double mag = cabs(G);
    resp->magnitude_db = (mag > 1e-15) ? 20.0 * log10(mag) : -300.0;
    resp->phase_deg    = atan2(cimag(G), creal(G)) * 180.0 / M_PI;
    return 0;
}

/*----------------------------------------------------------------------
 * L5 — Frequency Sweep (Bode Data)
 *----------------------------------------------------------------------*/

int tf_frequency_sweep(const TransferFunction *tf,
                        double omega_start,
                        double omega_end,
                        int n_points,
                        FrequencyResponse *responses)
{
    if (!tf || !responses || n_points < 2 ||
        omega_start <= 0.0 || omega_end <= omega_start)
    { errno = EINVAL; return -1; }

    double log_start = log10(omega_start);
    double log_end   = log10(omega_end);

    for (int i = 0; i < n_points; i++) {
        double frac = (double)i / (double)(n_points - 1);
        double w = pow(10.0, log_start + frac * (log_end - log_start));
        tf_frequency_response(tf, w, &responses[i]);
    }
    return 0;
}

/*----------------------------------------------------------------------
 * L5 — Resonant Peak
 *----------------------------------------------------------------------*/

int tf_resonant_peak(const TransferFunction *tf,
                      double *omega_r, double *M_r_db)
{
    if (!tf || !omega_r || !M_r_db) { errno = EINVAL; return -1; }

    /* Try analytic for 2nd-order */
    SecondOrderParams sop;
    if (tf_extract_second_order(tf, &sop) == 0 && sop.zeta < 0.707 && sop.zeta > 0.0) {
        *omega_r = sop.omega_n * sqrt(1.0 - 2.0 * sop.zeta * sop.zeta);
        if (*omega_r <= 0.0) return -1;
        double M_r = 1.0 / (2.0 * sop.zeta * sqrt(1.0 - sop.zeta * sop.zeta));
        *M_r_db = 20.0 * log10(M_r);
        return 0;
    }

    /* Numerical search */
    double best_mag = -1e300;
    double best_w   = 0.0;
    double w_low  = 1e-3;
    double w_high = 1e4;
    int N = 500;

    for (int i = 0; i < N; i++) {
        double w = w_low * pow(w_high / w_low, (double)i / (N - 1));
        FrequencyResponse fr;
        if (tf_frequency_response(tf, w, &fr) == 0) {
            if (fr.magnitude_db > best_mag) {
                best_mag = fr.magnitude_db;
                best_w   = w;
            }
        }
    }

    if (best_w <= 0.0) return -1;
    *omega_r = best_w;
    *M_r_db  = best_mag;
    return 0;
}

/*----------------------------------------------------------------------
 * L5 — Bandwidth
 *----------------------------------------------------------------------*/

int tf_bandwidth(const TransferFunction *tf, double *omega_bw)
{
    if (!tf || !omega_bw) { errno = EINVAL; return -1; }

    double dc = tf_dc_gain(tf);
    double dc_db = 20.0 * log10(fabs(dc));
    double target_db = dc_db - 3.0;

    /* Binary search for -3dB frequency */
    double low  = 1e-6;
    double high = 1e8;

    for (int iter = 0; iter < 80; iter++) {
        double mid = sqrt(low * high);
        FrequencyResponse fr;
        if (tf_frequency_response(tf, mid, &fr) != 0) return -1;
        if (fr.magnitude_db > target_db)
            low = mid;
        else
            high = mid;
    }

    *omega_bw = (low + high) / 2.0;
    return 0;
}

/*----------------------------------------------------------------------
 * L3/L5 — Pole-Zero Map
 *----------------------------------------------------------------------*/

int tf_pole_zero_map(const TransferFunction *tf,
                      LaplacePoleZero *pz)
{
    if (!tf || !pz) { errno = EINVAL; return -1; }

    /* Convert to LaplaceRational for reuse */
    LaplaceRational L;
    memset(&L, 0, sizeof(L));
    L.numerator.order   = tf->num_order;
    L.denominator.order = tf->den_order;
    for (int i = 0; i <= tf->num_order; i++)
        L.numerator.coeff[i] = tf->num_coeff[i];
    for (int i = 0; i <= tf->den_order; i++)
        L.denominator.coeff[i] = tf->den_coeff[i];

    memset(pz, 0, sizeof(*pz));
    pz->gain = tf->gain;

    if (laplace_find_poles(&L, pz) != 0) return -1;

    LaplacePoleZero pz_copy = *pz;
    laplace_find_zeros(&L, pz);
    /* Restore poles (find_zeros overwrites them) */
    pz->num_poles = pz_copy.num_poles;
    for (int i = 0; i < pz_copy.num_poles; i++)
        pz->poles[i] = pz_copy.poles[i];

    return 0;
}

/*----------------------------------------------------------------------
 * L8 — Minimum-Phase Check
 *----------------------------------------------------------------------*/

int tf_is_minimum_phase(const TransferFunction *tf)
{
    if (!tf) return 0;

    LaplacePoleZero pz;
    if (tf_pole_zero_map(tf, &pz) != 0) return 0;

    /* All zeros must have negative real part (CT) or |z| < 1 (DT) */
    for (int i = 0; i < pz.num_zeros; i++) {
        if (tf->sampling_period <= 0.0) {
            /* Continuous: check Re < 0 */
            if (creal(pz.zeros[i]) >= 0.0) return 0;
        } else {
            /* Discrete: check magnitude < 1 */
            if (cabs(pz.zeros[i]) >= 1.0) return 0;
        }
    }
    return 1;
}

/*----------------------------------------------------------------------
 * L2 — Relative Degree
 *----------------------------------------------------------------------*/

int tf_relative_degree(const TransferFunction *tf)
{
    if (!tf) return 0;
    return tf->den_order - tf->num_order;
}
