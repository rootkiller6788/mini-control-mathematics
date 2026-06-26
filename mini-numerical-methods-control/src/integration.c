/**
 * integration.c — ODE Integration Methods for Control Systems
 *
 * Implements: L4 Conservation (ODE as dynamical system), L5 Methods (RK family),
 *             L8 Advanced Methods (adaptive stepping, multistep)
 *
 * L7 Applications: quadrotor flight simulation (RK4 in X-Plane/JSBSim),
 * SpaceX Falcon 9 re-entry trajectory (adaptive RKF45), Mars entry
 * descent landing (EDL) simulation, Tesla Autopilot trajectory prediction.
 *
 * For control systems dx/dt = f(t,x,u), numerical integration is the bridge
 * between continuous-time models and discrete-time computation. This module
 * implements the algorithms that power simulation, MPC prediction, and
 * trajectory optimization.
 */

#include "integration.h"
#include <string.h>
#include <stdio.h>
#include <float.h>

/* ==========================================================================
 * L5: Forward Euler (Explicit, 1st order)
 * ========================================================================== */

int euler_forward_step(ODESystem *sys, double t, const double *x,
                       double h, double *x_next) {
    if (!sys || !sys->rhs || !x || !x_next) return -1;
    double *f = (double*)malloc(sys->n * sizeof(double));
    if (!f) return -1;
    sys->rhs(t, x, f, sys->params, sys->n);
    for (size_t i = 0; i < sys->n; i++)
        x_next[i] = x[i] + h * f[i];
    free(f);
    return 0;
}

/* ==========================================================================
 * L5: Backward Euler (Implicit, 1st order, A-stable)
 * ========================================================================== */

int euler_backward_step(ODESystem *sys, double t, const double *x,
                        double h, double *x_next, size_t max_newton, double tol) {
    if (!sys || !sys->rhs || !x || !x_next) return -1;
    size_t n = sys->n;
    double *f = (double*)malloc(n * sizeof(double));
    double *F = (double*)malloc(n * sizeof(double));
    if (!f || !F) { free(f); free(F); return -1; }

    /* Initial guess: forward Euler */
    sys->rhs(t, x, f, sys->params, n);
    for (size_t i = 0; i < n; i++)
        x_next[i] = x[i] + h * f[i];

    /* Newton iteration to solve x_next - h*f(t+h, x_next) - x = 0 */
    for (size_t k = 0; k < max_newton; k++) {
        sys->rhs(t + h, x_next, f, sys->params, n);
        double norm = 0.0;
        for (size_t i = 0; i < n; i++) {
            F[i] = x_next[i] - h * f[i] - x[i];
            norm += F[i] * F[i];
        }
        if (sqrt(norm) < tol) break;
        /* Simplified: use Jacobian = I - h*∂f/∂x ≈ I (ignoring off-diagonals) */
        for (size_t i = 0; i < n; i++)
            x_next[i] -= F[i];
    }
    free(f); free(F);
    return 0;
}

/* ==========================================================================
 * L5: Heun's Method (2nd order, explicit)
 * ========================================================================== */

int heun_step(ODESystem *sys, double t, const double *x,
              double h, double *x_next) {
    if (!sys || !sys->rhs || !x || !x_next) return -1;
    size_t n = sys->n;
    double *k1 = (double*)malloc(n * sizeof(double));
    double *k2 = (double*)malloc(n * sizeof(double));
    double *xt = (double*)malloc(n * sizeof(double));
    if (!k1 || !k2 || !xt) { free(k1); free(k2); free(xt); return -1; }

    sys->rhs(t, x, k1, sys->params, n);
    for (size_t i = 0; i < n; i++) xt[i] = x[i] + h * k1[i];
    sys->rhs(t + h, xt, k2, sys->params, n);
    for (size_t i = 0; i < n; i++)
        x_next[i] = x[i] + 0.5 * h * (k1[i] + k2[i]);

    free(k1); free(k2); free(xt);
    return 0;
}

/* ==========================================================================
 * L5: Classical RK4 (4th order, 4 stages)
 * ========================================================================== */

int rk4_step(ODESystem *sys, double t, const double *x,
             double h, double *x_next) {
    if (!sys || !sys->rhs || !x || !x_next) return -1;
    size_t n = sys->n;
    double *k1 = (double*)malloc(4 * n * sizeof(double));
    double *k2 = k1 + n, *k3 = k2 + n, *k4 = k3 + n;
    double *xt = (double*)malloc(n * sizeof(double));
    if (!k1 || !xt) { free(k1); free(xt); return -1; }

    /* k1 = f(t, x) */
    sys->rhs(t, x, k1, sys->params, n);

    /* k2 = f(t + h/2, x + h*k1/2) */
    for (size_t i = 0; i < n; i++) xt[i] = x[i] + 0.5 * h * k1[i];
    sys->rhs(t + 0.5 * h, xt, k2, sys->params, n);

    /* k3 = f(t + h/2, x + h*k2/2) */
    for (size_t i = 0; i < n; i++) xt[i] = x[i] + 0.5 * h * k2[i];
    sys->rhs(t + 0.5 * h, xt, k3, sys->params, n);

    /* k4 = f(t + h, x + h*k3) */
    for (size_t i = 0; i < n; i++) xt[i] = x[i] + h * k3[i];
    sys->rhs(t + h, xt, k4, sys->params, n);

    /* x_{n+1} = x_n + (h/6)*(k1 + 2k2 + 2k3 + k4) */
    for (size_t i = 0; i < n; i++)
        x_next[i] = x[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);

    free(k1); free(xt);
    return 0;
}

/* ==========================================================================
 * L5/L8: Runge-Kutta-Fehlberg 5(4) — Adaptive Step
 * ========================================================================== */

/* Butcher tableau coefficients for RKF45 (Dormand-Prince 5(4) variant) */
static const double rkf_a21 = 1.0/4.0;
static const double rkf_a31 = 3.0/32.0,  rkf_a32 = 9.0/32.0;
static const double rkf_a41 = 1932.0/2197.0, rkf_a42 = -7200.0/2197.0, rkf_a43 = 7296.0/2197.0;
static const double rkf_a51 = 439.0/216.0, rkf_a52 = -8.0, rkf_a53 = 3680.0/513.0, rkf_a54 = -845.0/4104.0;
static const double rkf_a61 = -8.0/27.0, rkf_a62 = 2.0, rkf_a63 = -3544.0/2565.0, rkf_a64 = 1859.0/4104.0, rkf_a65 = -11.0/40.0;

/* 5th-order weights */
static const double rkf_b5_1 = 16.0/135.0, rkf_b5_3 = 6656.0/12825.0, rkf_b5_4 = 28561.0/56430.0, rkf_b5_5 = -9.0/50.0, rkf_b5_6 = 2.0/55.0;
/* 4th-order weights */
static const double rkf_b4_1 = 25.0/216.0, rkf_b4_3 = 1408.0/2565.0, rkf_b4_4 = 2197.0/4104.0, rkf_b4_5 = -1.0/5.0;

int rkf45_step(ODESystem *sys, double t, const double *x,
               double h, double *x_next, double *error_est,
               double *h_next, double rtol, double atol) {
    if (!sys || !sys->rhs || !x || !x_next) return -1;
    size_t n = sys->n;
    double *k1 = (double*)malloc(6 * n * sizeof(double));
    double *k2 = k1 + n, *k3 = k2 + n, *k4 = k3 + n;
    double *k5 = k4 + n, *k6 = k5 + n;
    double *xt = (double*)malloc(n * sizeof(double));
    if (!k1 || !xt) { free(k1); free(xt); return -1; }

    /* Stage 1 */
    sys->rhs(t, x, k1, sys->params, n);

    /* Stage 2 */
    for (size_t i = 0; i < n; i++) xt[i] = x[i] + h * rkf_a21 * k1[i];
    sys->rhs(t + 0.25 * h, xt, k2, sys->params, n);

    /* Stage 3 */
    for (size_t i = 0; i < n; i++) xt[i] = x[i] + h * (rkf_a31 * k1[i] + rkf_a32 * k2[i]);
    sys->rhs(t + 0.375 * h, xt, k3, sys->params, n);

    /* Stage 4 */
    for (size_t i = 0; i < n; i++) xt[i] = x[i] + h * (rkf_a41 * k1[i] + rkf_a42 * k2[i] + rkf_a43 * k3[i]);
    sys->rhs(t + 12.0/13.0 * h, xt, k4, sys->params, n);

    /* Stage 5 */
    for (size_t i = 0; i < n; i++) xt[i] = x[i] + h * (rkf_a51 * k1[i] + rkf_a52 * k2[i] + rkf_a53 * k3[i] + rkf_a54 * k4[i]);
    sys->rhs(t + h, xt, k5, sys->params, n);

    /* Stage 6 */
    for (size_t i = 0; i < n; i++) xt[i] = x[i] + h * (rkf_a61 * k1[i] + rkf_a62 * k2[i] + rkf_a63 * k3[i] + rkf_a64 * k4[i] + rkf_a65 * k5[i]);
    sys->rhs(t + 0.5 * h, xt, k6, sys->params, n);

    /* Compute 4th and 5th order solutions */
    double max_err = 0.0;
    for (size_t i = 0; i < n; i++) {
        double x5 = x[i] + h * (rkf_b5_1 * k1[i] + rkf_b5_3 * k3[i] + rkf_b5_4 * k4[i] + rkf_b5_5 * k5[i] + rkf_b5_6 * k6[i]);
        double x4 = x[i] + h * (rkf_b4_1 * k1[i] + rkf_b4_3 * k3[i] + rkf_b4_4 * k4[i] + rkf_b4_5 * k5[i]);
        x_next[i] = x5;  /* Use 5th-order as output */
        double err = fabs(x5 - x4);
        double scale = atol + rtol * fmax(fabs(x[i]), fabs(x5));
        double scaled_err = err / scale;
        if (scaled_err > max_err) max_err = scaled_err;
    }
    *error_est = max_err;

    /* Step size control: PI controller (Gustafsson) */
    if (max_err > 0.0) {
        double factor = 0.9 * pow(1.0 / max_err, 0.2);  /* 1/5 for 5th order */
        if (factor > 4.0) factor = 4.0;
        if (factor < 0.1) factor = 0.1;
        *h_next = h * factor;
    } else {
        *h_next = h * 2.0;
    }

    free(k1); free(xt);
    return 0;
}

/* ==========================================================================
 * L5: Adams-Bashforth Multistep Methods
 * ========================================================================== */

int adams_bashforth2_step(ODESystem *sys, double t, double h,
                          const double *x_prev, const double *f_prev,
                          const double *x_curr, double *x_next) {
    if (!sys || !sys->rhs || !x_prev || !f_prev || !x_curr || !x_next) return -1;
    size_t n = sys->n;
    double *f_curr = (double*)malloc(n * sizeof(double));
    if (!f_curr) return -1;
    sys->rhs(t, x_curr, f_curr, sys->params, n);
    /* AB2: x_{n+2} = x_{n+1} + h * (3/2*f_{n+1} - 1/2*f_n) */
    for (size_t i = 0; i < n; i++)
        x_next[i] = x_curr[i] + h * (1.5 * f_curr[i] - 0.5 * f_prev[i]);
    free(f_curr);
    return 0;
}

int adams_bashforth3_step(ODESystem *sys, double t, double h,
                          const double *f_n, const double *f_n1,
                          const double *f_n2, const double *x_curr,
                          double *x_next) {
    if (!sys || !sys->rhs || !f_n || !f_n1 || !f_n2 || !x_curr || !x_next) return -1;
    size_t n = sys->n;
    double *f_curr = (double*)malloc(n * sizeof(double));
    if (!f_curr) return -1;
    sys->rhs(t, x_curr, f_curr, sys->params, n);
    (void)f_n; (void)f_n1; (void)f_n2;
    /* AB3: coefficients 23/12, -16/12, 5/12 using actual stored derivatives */
    for (size_t i = 0; i < n; i++)
        x_next[i] = x_curr[i] + h * (23.0/12.0 * f_curr[i] - 16.0/12.0 * f_n2[i] + 5.0/12.0 * f_n1[i]);
    free(f_curr);
    return 0;
}

int adams_bashforth4_step(ODESystem *sys, double t, double h,
                          const double *f_n, const double *f_n1,
                          const double *f_n2, const double *f_n3,
                          const double *x_curr, double *x_next) {
    if (!sys || !sys->rhs || !f_n || !f_n1 || !f_n2 || !f_n3 || !x_curr || !x_next) return -1;
    size_t n = sys->n;
    double *f_curr = (double*)malloc(n * sizeof(double));
    if (!f_curr) return -1;
    sys->rhs(t, x_curr, f_curr, sys->params, n);
    (void)f_n; (void)f_n3;
    /* AB4: x_{n+1} = x_n + h*(55/24*f_n - 59/24*f_{n-1} + 37/24*f_{n-2} - 9/24*f_{n-3}) */
    for (size_t i = 0; i < n; i++)
        x_next[i] = x_curr[i] + h * (55.0/24.0 * f_curr[i] - 59.0/24.0 * f_n2[i] + 37.0/24.0 * f_n1[i] - 9.0/24.0 * f_n[i]);
    free(f_curr);
    return 0;
}

/* ==========================================================================
 * L5: Fixed-Step Integration Driver
 * ========================================================================== */

ODESolution* solve_ode_fixed_step(ODESystem *sys, double t0, double tf,
                                  const double *x0, double h) {
    if (!sys || !sys->rhs || !x0 || h <= 0.0 || tf <= t0) return NULL;
    if (h > tf - t0) h = tf - t0;

    size_t n = sys->n;
    size_t N = (size_t)((tf - t0) / h) + 1;
    if (N < 2) N = 2;

    ODESolution *sol = (ODESolution*)malloc(sizeof(ODESolution));
    if (!sol) return NULL;
    sol->n_states = n;
    sol->n_steps = N;
    sol->time = (double*)malloc(N * sizeof(double));
    sol->states = (double*)malloc(N * n * sizeof(double));
    if (!sol->time || !sol->states) { odesolution_free(sol); return NULL; }

    /* Initial condition */
    sol->time[0] = t0;
    memcpy(sol->states, x0, n * sizeof(double));

    double t = t0;
    double *x_curr = (double*)malloc(n * sizeof(double));
    double *x_next = (double*)malloc(n * sizeof(double));
    if (!x_curr || !x_next) { free(x_curr); free(x_next); odesolution_free(sol); return NULL; }
    memcpy(x_curr, x0, n * sizeof(double));

    for (size_t k = 1; k < N; k++) {
        t = t0 + (k - 1) * h;
        double h_use = (k == N - 1 && t + h > tf) ? (tf - t) : h;
        rk4_step(sys, t, x_curr, h_use, x_next);
        sol->time[k] = t + h_use;
        memcpy(sol->states + k * n, x_next, n * sizeof(double));
        memcpy(x_curr, x_next, n * sizeof(double));
    }

    free(x_curr); free(x_next);
    return sol;
}

/* ==========================================================================
 * L5/L8: Adaptive Integration Driver
 * ========================================================================== */

ODESolution* solve_ode_adaptive(ODESystem *sys, double t0, double tf,
                                const double *x0, double h_init,
                                double rtol, double atol,
                                ODEStatistics *stats) {
    if (!sys || !sys->rhs || !x0 || h_init <= 0.0 || tf <= t0) return NULL;
    size_t n = sys->n;

    /* Pre-allocate with reasonable initial capacity */
    size_t cap = 1024;
    ODESolution *sol = (ODESolution*)malloc(sizeof(ODESolution));
    if (!sol) return NULL;
    sol->n_states = n;
    sol->n_steps = 0;
    sol->time = (double*)malloc(cap * sizeof(double));
    sol->states = (double*)malloc(cap * n * sizeof(double));
    if (!sol->time || !sol->states) { odesolution_free(sol); return NULL; }

    double *x_curr = (double*)malloc(n * sizeof(double));
    double *x_next = (double*)malloc(n * sizeof(double));
    if (!x_curr || !x_next) { free(x_curr); free(x_next); odesolution_free(sol); return NULL; }
    memcpy(x_curr, x0, n * sizeof(double));

    if (stats) memset(stats, 0, sizeof(ODEStatistics));

    double t = t0, h = h_init;
    size_t step_count = 0, rej_count = 0;

    /* Store initial point */
    sol->time[0] = t0;
    memcpy(sol->states, x0, n * sizeof(double));
    sol->n_steps = 1;

    while (t < tf) {
        if (t + h > tf) h = tf - t;
        if (h < MNC_TINY) break;

        double err_est, h_next;
        rkf45_step(sys, t, x_curr, h, x_next, &err_est, &h_next, rtol, atol);
        step_count++;

        if (err_est <= 1.0 || h <= 1e-12) {
            /* Accept step */
            t += h;
            memcpy(x_curr, x_next, n * sizeof(double));

            /* Grow arrays if needed */
            if (sol->n_steps >= cap) {
                cap *= 2;
                sol->time = (double*)realloc(sol->time, cap * sizeof(double));
                sol->states = (double*)realloc(sol->states, cap * n * sizeof(double));
                if (!sol->time || !sol->states) break;
            }
            sol->time[sol->n_steps] = t;
            memcpy(sol->states + sol->n_steps * n, x_curr, n * sizeof(double));
            sol->n_steps++;
        } else {
            rej_count++;
        }
        h = h_next;
    }

    if (stats) {
        stats->n_steps_total = step_count;
        stats->n_steps_rejected = rej_count;
    }

    free(x_curr); free(x_next);
    return sol;
}

void odesolution_free(ODESolution *sol) {
    if (!sol) return;
    free(sol->time);
    free(sol->states);
    free(sol);
}

int odesolution_interpolate(const ODESolution *sol, double t, double *x_out) {
    if (!sol || !x_out || sol->n_steps < 2) return -1;
    /* Find interval */
    if (t < sol->time[0] || t > sol->time[sol->n_steps - 1]) return -1;
    size_t k = 0;
    for (size_t i = 0; i < sol->n_steps - 1; i++) {
        if (sol->time[i] <= t && t <= sol->time[i + 1]) { k = i; break; }
    }
    double *xk = sol->states + k * sol->n_states;
    double *xk1 = sol->states + (k + 1) * sol->n_states;
    double dt = sol->time[k + 1] - sol->time[k];
    if (dt < MNC_TINY) { memcpy(x_out, xk, sol->n_states * sizeof(double)); return 0; }
    double theta = (t - sol->time[k]) / dt;
    for (size_t i = 0; i < sol->n_states; i++)
        x_out[i] = (1.0 - theta) * xk[i] + theta * xk1[i];
    return 0;
}
