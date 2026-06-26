/**
 * @example example_harmonic_oscillator.c
 * @brief  End-to-end example: damped harmonic oscillator analysis.
 *
 * Demonstrates:
 *   - IVP creation and fixed-step RK4 integration
 *   - Frequency response computation
 *   - Resonance detection
 *   - Stability analysis via eigenvalues
 *
 * Physics: m·ẍ + c·ẋ + k·x = 0
 * With m=1, k=10, and varying c (damping).
 */

#include "ode_types.h"
#include "ode_ivp.h"
#include "ode_numerical.h"
#include "ode_nonlinear.h"
#include "ode_applications.h"
#include "ode_linear.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(void) {
    printf("=== Harmonic Oscillator Analysis ===\n\n");

    /* Parameters: m=1, k=10, test three damping cases */
    double m = 1.0, k = 10.0;
    double damping_cases[] = {0.0, 0.5, 6.5};  /* undamped, underdamped, overdamped */
    const char *labels[] = {"Undamped (c=0)", "Underdamped (c=0.5)", "Overdamped (c=6.32)"};
    double omega_n = sqrt(k / m);
    printf("Natural frequency omega_n = %.4f rad/s\n", omega_n);
    printf("Critical damping c_crit = %.4f\n\n", 2.0 * sqrt(m * k));

    for (int ic = 0; ic < 3; ic++) {
        OscillatorParams p = {m, damping_cases[ic], k, 0.0, 0.0};
        double zeta = p.c / (2.0 * sqrt(m * k));

        printf("--- Case %d: %s (zeta=%.4f) ---\n", ic+1, labels[ic], zeta);

        /* Build IVP: ẏ = f(t,y) with y = [x, v] */
        double y0[2] = {1.0, 0.0};  /* Initial displacement = 1 */
        ODEIVP ivp = ode_ivp_create(ode_rhs_harmonic_oscillator, NULL, &p,
                                      0.0, 5.0, y0, 2, ODE_LINEAR_CONSTANT_COEFF);

        /* Integrate */
        ODESolution sol;
        memset(&sol, 0, sizeof(sol));
        int ret = ode_integrate_fixed_step(&ivp, METHOD_RK4_CLASSICAL, 0.01, &sol);
        if (ret != 0) {
            printf("  Integration failed\n");
            ode_ivp_destroy(&ivp);
            continue;
        }

        /* Print trajectory */
        printf("  t       x(t)        v(t)\n");
        int step = sol.n_points / 5;
        for (int i = 0; i < sol.n_points; i += step) {
            printf("  %.2f   %+.6f   %+.6f\n",
                   sol.t_vals[i], sol.y_vals[i*2], sol.y_vals[i*2+1]);
        }
        printf("  ... (%d total points)\n", sol.n_points);

        /* Frequency response */
        double freqs[] = {0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0};
        double M[7], phase[7];
        ode_oscillator_frequency_response(&p, freqs, 7, M, phase);
        printf("  Frequency response (omega/omega_n=%.2f):\n", omega_n);
        printf("  omega    M(omega)   phase(deg)\n");
        for (int i = 0; i < 7; i++) {
            printf("  %5.2f   %8.4f   %+8.2f\n", freqs[i], M[i], phase[i]*180.0/M_PI);
        }

        /* Resonance */
        double omega_res, M_max;
        if (ode_find_resonance(&p, &omega_res, &M_max)) {
            printf("  Resonance at omega=%.4f, M_max=%.4f\n", omega_res, M_max);
        } else {
            printf("  No resonance peak (zeta=%.4f >= 0.707)\n", zeta);
        }

        ode_solution_free(&sol);
        ode_ivp_destroy(&ivp);
        printf("\n");
    }

    return 0;
}
