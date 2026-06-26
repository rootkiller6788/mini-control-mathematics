/**
 * @example example_lorenz_chaos.c
 * @brief  End-to-end example: Lorenz system simulation showing chaos.
 *
 * Demonstrates:
 *   - Lorenz attractor integration
 *   - Sensitivity to initial conditions (hallmark of chaos)
 *   - Lyapunov exponent estimation
 *   - Equilibrium analysis
 *
 * Reference: Lorenz "Deterministic Nonperiodic Flow" (1963)
 *            Strogatz "Nonlinear Dynamics and Chaos" (2015), Chap. 9
 *
 * Application: Atmospheric convection model, chaotic systems (L7, L8).
 */

#include "ode_types.h"
#include "ode_ivp.h"
#include "ode_numerical.h"
#include "ode_nonlinear.h"
#include "ode_applications.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int lorenz_rhs(double t, const double *y, double *f, int dim, void *ctx) {
    (void)t;
    LorenzParams *p = (LorenzParams *)ctx;
    f[0] = p->sigma * (y[1] - y[0]);
    f[1] = y[0] * (p->rho - y[2]) - y[1];
    f[2] = y[0] * y[1] - p->beta * y[2];
    return 0;
}

int main(void) {
    printf("=== Lorenz Attractor: Chaos Analysis ===\n\n");

    /* Classic parameters (chaotic regime) */
    LorenzParams params = {10.0, 28.0, 8.0/3.0};
    printf("Parameters: sigma=%.1f, rho=%.1f, beta=%.4f\n",
           params.sigma, params.rho, params.beta);

    /* Equilibrium points */
    double eq[9];
    int n_eq;
    ode_lorenz_equilibria(&params, eq, &n_eq);
    printf("\nEquilibrium points (n=%d):\n", n_eq);
    for (int i = 0; i < n_eq; i++) {
        printf("  C%d = (%+.4f, %+.4f, %+.4f)\n",
               i, eq[i*3], eq[i*3+1], eq[i*3+2]);
    }

    /* Simulate from a point near the origin */
    double y0[3] = {1.0, 1.0, 1.0};
    double t_end = 40.0;
    int n_steps = 4000;

    double *x = (double *)malloc((size_t)n_steps * sizeof(double));
    double *y_t = (double *)malloc((size_t)n_steps * sizeof(double));
    double *z = (double *)malloc((size_t)n_steps * sizeof(double));
    double *t_vals = (double *)malloc((size_t)n_steps * sizeof(double));

    int ret = ode_lorenz_simulate(&params, y0[0], y0[1], y0[2],
                                    t_end, n_steps, x, y_t, z, t_vals);
    if (ret != 0) {
        printf("Simulation failed!\n");
        free(x); free(y_t); free(z); free(t_vals);
        return 1;
    }

    printf("\nSimulation: t=[0, %.0f], %d steps\n", t_end, n_steps);
    printf("Trajectory snapshot:\n");
    printf("  t       x(t)         y(t)         z(t)\n");
    for (int i = 0; i < n_steps; i += n_steps/8) {
        printf("  %5.1f   %+9.4f   %+9.4f   %+9.4f\n",
               t_vals[i], x[i], y_t[i], z[i]);
    }

    /* Sensitivity to initial conditions */
    printf("\n--- Sensitivity to Initial Conditions ---\n");
    double y0_pert[3] = {1.0 + 1e-6, 1.0, 1.0};  /* Tiny perturbation */
    double *x2 = (double *)malloc((size_t)n_steps * sizeof(double));
    double *y2 = (double *)malloc((size_t)n_steps * sizeof(double));
    double *z2 = (double *)malloc((size_t)n_steps * sizeof(double));
    double *t2 = (double *)malloc((size_t)n_steps * sizeof(double));

    ode_lorenz_simulate(&params, y0_pert[0], y0_pert[1], y0_pert[2],
                          t_end, n_steps, x2, y2, z2, t2);

    /* Report divergence */
    printf("Initial separation: |delta| = 1e-6\n");
    int check_pts[] = {10, 50, 100, 200, 500, 1000, 2000, 3999};
    for (int i = 0; i < 8; i++) {
        int idx = check_pts[i];
        double dx = x[idx] - x2[idx];
        double dy = y_t[idx] - y2[idx];
        double dz = z[idx] - z2[idx];
        double dist = sqrt(dx*dx + dy*dy + dz*dz);
        printf("  t=%.1f: separation = %.2e\n", t_vals[idx], dist);
    }

    /* Lyapunov exponent estimation */
    printf("\n--- Lyapunov Exponent ---\n");
    VectorField vf;
    vf.f = lorenz_rhs;
    vf.jac = NULL;
    vf.ctx = &params;
    vf.dim = 3;

    double lambda_max;
    ret = ode_lyapunov_exponent(&vf, y0, 1e-6, 50.0, 500, &lambda_max);
    if (ret == 0) {
        printf("Maximal Lyapunov exponent λ_max ≈ %+.4f\n", lambda_max);
        if (lambda_max > 0.0) {
            printf("  λ_max > 0 → CHAOTIC (sensitive dependence on IC)\n");
        } else if (lambda_max < 0.0) {
            printf("  λ_max < 0 → STABLE (attracting fixed point/cycle)\n");
        }
    }

    free(x); free(y_t); free(z); free(t_vals);
    free(x2); free(y2); free(z2); free(t2);
    return 0;
}
