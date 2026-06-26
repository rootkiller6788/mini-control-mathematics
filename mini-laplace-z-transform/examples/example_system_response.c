/**
 * @file example_system_response.c
 * @brief Second-Order System — Step and Impulse Response Analysis
 *
 * L6: Canonical second-order system response analysis.
 *     Demonstrates the complete transient specification extraction.
 *     Covers underdamped, critically damped, and overdamped cases.
 *
 * Tesla motor control: real-world underdamped second-order response analysis
 * for servo positioning.
 */
#include "../include/system_response.h"
#include "../include/transfer_function.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    printf("=== Second-Order System Response Analysis ===\n");
    printf("Tesla servo positioning system — 2nd-order model\n\n");

    /* Three damping cases */
    double zeta_cases[] = {0.3, 0.7, 1.0, 2.0};
    double omega_n = 10.0;  /* natural frequency = 10 rad/s (~1.6 Hz) */
    const char *labels[] = {"Underdamped (ζ=0.3)", "Underdamped (ζ=0.7)",
                            "Critically damped (ζ=1.0)", "Overdamped (ζ=2.0)"};
    int n_cases = 4;

    for (int c = 0; c < n_cases; c++) {
        double zeta = zeta_cases[c];

        printf("--- %s ---\n", labels[c]);
        printf("ωₙ = %.1f rad/s, ζ = %.1f\n", omega_n, zeta);

        /* Closed-form 2nd-order specs */
        TransientSpecs specs;
        second_order_specs_closed_form(omega_n, zeta, &specs);

        printf("  Rise time (0-100%%):  %.4f s\n", specs.rise_time);
        printf("  Peak time:           %.4f s\n", specs.peak_time);
        printf("  Settling time (2%%):  %.4f s\n", specs.settling_time);
        printf("  Overshoot:           %.1f %%\n", specs.overshoot_percent);
        printf("  Delay time:          %.4f s\n", specs.delay_time);

        /* Step response at key points */
        double t_points[] = {0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.5};
        int n_t = sizeof(t_points) / sizeof(t_points[0]);

        printf("  Step response y(t):\n");
        for (int i = 0; i < n_t; i++) {
            double t = t_points[i];
            double y = second_order_step_response_at(omega_n, zeta, t);
            printf("    t=%.2f s → y=%.4f\n", t, y);
        }
        printf("\n");
    }

    /* Build transfer function and compute time response for ζ=0.7 */
    printf("--- Numerical Time Response (ζ=0.7) ---\n");
    TransferFunction tf;
    tf.num_order = 0;
    tf.den_order = 2;
    double z = 0.7;
    double wn = omega_n;
    tf.num_coeff[0] = wn * wn;
    tf.den_coeff[0] = wn * wn;
    tf.den_coeff[1] = 2.0 * z * wn;
    tf.den_coeff[2] = 1.0;
    tf.gain = 1.0;
    tf.sampling_period = 0.0;

    TimeResponse resp;
    int ret = compute_time_response(&tf, RESPONSE_STEP, 0.0, 1.5, 200, &resp);
    if (ret == 0) {
        TransientSpecs num_specs;
        extract_transient_specs(&resp, &num_specs);

        printf("Numerical step response (200 samples):\n");
        printf("  Rise time:           %.4f s\n", num_specs.rise_time);
        printf("  Peak time:           %.4f s\n", num_specs.peak_time);
        printf("  Settling time:       %.4f s\n", num_specs.settling_time);
        printf("  Overshoot:           %.1f %%\n", num_specs.overshoot_percent);
        printf("  Steady-state value:  %.4f\n", num_specs.steady_state_value);

        /* Compare with closed-form */
        TransientSpecs closed;
        second_order_specs_closed_form(omega_n, 0.7, &closed);
        printf("\n  Comparison (numerical vs closed-form):\n");
        printf("  Overshoot: %.1f%% vs %.1f%%\n",
               num_specs.overshoot_percent, closed.overshoot_percent);
        printf("  Settling:  %.4f vs %.4f s\n",
               num_specs.settling_time, closed.settling_time);

        time_response_free(&resp);
    }

    /* Steady-state error analysis */
    printf("\n--- Steady-State Error Analysis ---\n");
    TransferFunction plant;
    plant.num_order = 0; plant.num_coeff[0] = 1.0;
    plant.den_order = 1; plant.den_coeff[0] = 0.0; plant.den_coeff[1] = 1.0;
    plant.gain = 2.0; plant.sampling_period = 0.0;

    double ess_step, ess_ramp;
    steady_state_error(&plant, RESPONSE_STEP, &ess_step);
    steady_state_error(&plant, RESPONSE_RAMP, &ess_ramp);

    printf("Plant G(s) = 2/s (Type-1):\n");
    printf("  Step error e_ss = %.4f (expect 0 for Type-1)\n", ess_step);
    printf("  Ramp error e_ss = %.4f (expect 1/Kv = 0.5)\n", ess_ramp);

    printf("\n=== Analysis Complete ===\n");
    return 0;
}
