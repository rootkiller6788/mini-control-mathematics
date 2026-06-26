/**
 * @file example_motor_control.c
 * @brief DC Motor Speed Control — Laplace Domain Analysis
 *
 * L6: Classic second-order system analysis for DC motor speed control.
 *     Transfer function: G(s) = K / (Js + B)(Ls + R)
 *     Approximated as 2nd-order: G(s) = ωₙ²/(s² + 2ζωₙs + ωₙ²)
 *
 * L7: Application — physical motor parameters from a real DC motor.
 *      Parameters based on a typical 12V DC motor (Maxon RE 25).
 */
#include "../include/transfer_function.h"
#include "../include/system_response.h"
#include "../include/stability.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    printf("=== DC Motor Speed Control Analysis ===\n");
    printf("Motor: typical 12V brushed DC (Maxon RE 25 class)\n\n");

    /* Motor parameters */
    double J = 1.13e-6;   /* rotor inertia (kg·m²) */
    double B = 1.0e-6;    /* viscous friction (N·m·s/rad) */
    double Kt = 0.0234;   /* torque constant (N·m/A) */
    double Ke = 0.0234;   /* back EMF constant (V·s/rad) */
    double R = 0.85;      /* armature resistance (Ω) */
    double L = 0.16e-3;   /* armature inductance (H) */

    /* Electrical time constant: τ_e = L/R */
    double tau_e = L / R;
    printf("Electrical time constant τ_e = %.2f ms\n", tau_e * 1e3);

    /* Mechanical time constant: τ_m = JR/(Kt·Ke) */
    double tau_m = J * R / (Kt * Ke);
    printf("Mechanical time constant τ_m = %.1f ms\n", tau_m * 1e3);

    /* Reduced 2nd-order: G(s) = K / ((τ_e·s+1)(τ_m·s+1)) */
    /* → G(s) = K/(τ_e·τ_m) / (s² + (1/τ_e + 1/τ_m)s + 1/(τ_e·τ_m)) */

    double omega_n = 1.0 / sqrt(tau_e * tau_m);
    double zeta   = 0.5 * (tau_e + tau_m) / sqrt(tau_e * tau_m);

    printf("\nSecond-order parameters:\n");
    printf("  ωₙ = %.1f rad/s  (%.2f Hz)\n", omega_n, omega_n / (2.0 * M_PI));
    printf("  ζ  = %.3f\n", zeta);

    /* Step response analysis */
    TransientSpecs specs;
    second_order_specs_closed_form(omega_n, zeta, &specs);

    printf("\nStep response specifications:\n");
    printf("  Rise time (10%%-90%%): %.3f ms\n", specs.rise_time * 1e3);
    printf("  Peak time:            %.3f ms\n", specs.peak_time * 1e3);
    printf("  Settling time (2%%):   %.3f ms\n", specs.settling_time * 1e3);
    printf("  Overshoot:            %.1f %%\n", specs.overshoot_percent);

    /* Build transfer function */
    TransferFunction motor;
    motor.num_order = 0;
    motor.den_order = 2;
    motor.num_coeff[0] = 1.0 / (tau_e * tau_m);
    motor.den_coeff[0] = 1.0 / (tau_e * tau_m);
    motor.den_coeff[1] = 2.0 * zeta * omega_n;
    motor.den_coeff[2] = 1.0;
    motor.gain = 1.0;
    motor.sampling_period = 0.0;

    /* Poles */
    LaplacePoleZero pz;
    tf_pole_zero_map(&motor, &pz);
    printf("\nPoles:\n");
    for (int i = 0; i < pz.num_poles; i++) {
        printf("  s%d = %.3f %+.3fj\n", i+1,
               creal(pz.poles[i]), cimag(pz.poles[i]));
    }

    /* Stability */
    LaplacePolynomial den;
    den.order = 2;
    den.coeff[0] = motor.den_coeff[0];
    den.coeff[1] = motor.den_coeff[1];
    den.coeff[2] = motor.den_coeff[2];
    RouthResult routh;
    routh_stability(&den, &routh);
    printf("\nStability (Routh-Hurwitz): %s\n",
           routh.is_stable ? "STABLE" : "UNSTABLE");

    /* Frequency response */
    FrequencyResponse fr;
    tf_frequency_response(&motor, omega_n, &fr);
    printf("\nFrequency response at ω = ωₙ = %.1f rad/s:\n", omega_n);
    printf("  Magnitude: %.1f dB\n", fr.magnitude_db);
    printf("  Phase:     %.1f deg\n", fr.phase_deg);

    printf("\n=== Analysis Complete ===\n");
    return 0;
}
