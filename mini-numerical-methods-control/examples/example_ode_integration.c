#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/numerical_core.h"
#include "../include/integration.h"

static void vanderpol(double t, const double *s, double *d, void *p, size_t n) {
    (void)t; (void)n;
    double mu = *((double*)p);
    d[0] = s[1];
    d[1] = mu * (1.0 - s[0]*s[0]) * s[1] - s[0];
}

int main(void) {
    printf("=== ODE Integration: Van der Pol Oscillator ===\n\n");

    double mu = 2.0;
    ODESystem sys = {vanderpol, &mu, 2};
    double x0[] = {2.0, 0.0};

    printf("Integrating with RK4, h=0.01, tf=20.0...\n");
    ODESolution *sol = solve_ode_fixed_step(&sys, 0.0, 20.0, x0, 0.01);

    if (sol) {
        size_t last = sol->n_steps - 1;
        printf("Completed %zu steps\n", sol->n_steps);
        printf("Final state: t=%.1f, x=[%.4f, %.4f]\n",
               sol->time[last],
               sol->states[last*2],
               sol->states[last*2+1]);
        printf("\nTrajectory samples:\n");
        for (size_t i = 0; i < sol->n_steps; i += sol->n_steps / 10) {
            printf("  t=%.2f: [%.4f, %.4f]\n",
                   sol->time[i], sol->states[i*2], sol->states[i*2+1]);
        }
        odesolution_free(sol);
    }

    printf("\nApplication: Quadrotor attitude simulation uses RK4 (X-Plane, JSBSim).\n");
    return 0;
}
