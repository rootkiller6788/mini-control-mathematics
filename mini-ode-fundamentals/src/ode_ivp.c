/**
 * @file    ode_ivp.c
 * @brief   Initial value problem implementation: creation, validation,
 *          Lipschitz estimation, Picard iteration, Gronwall inequality.
 *
 * L1-L4: Definitions through Fundamental Laws for IVPs.
 *
 * Each function implements an independent knowledge point:
 *   ode_ivp_create          → L1: IVP struct initialization
 *   ode_ivp_validate        → L1: Well-posedness check
 *   ode_lipschitz_estimate  → L2: Lipschitz continuity (Picard-Lindelöf prerequisite)
 *   ode_one_sided_lipschitz → L2: One-sided Lipschitz (contractivity)
 *   ode_picard_iteration    → L4: Picard existence proof construction
 *   ode_picard_solve        → L4: Full Picard iterative solver
 *   ode_gronwall_bound      → L4: Gronwall's inequality (growth bounds)
 *   ode_continuous_dependence_bound → L4: Solution continuous dependence on IC
 */

#include "ode_ivp.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>

/* ──────────────────────── Helper: Vector Math ──────────────────────── */

static double vec_norm_inf(const double *v, int n) {
    double max_val = 0.0;
    for (int i = 0; i < n; i++) {
        double abs_i = fabs(v[i]);
        if (abs_i > max_val) max_val = abs_i;
    }
    return max_val;
}

static double vec_norm_2(const double *v, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += v[i] * v[i];
    return sqrt(sum);
}

static void vec_copy(const double *src, double *dst, int n) {
    memcpy(dst, src, (size_t)n * sizeof(double));
}

static void vec_sub(const double *a, const double *b, double *c, int n) {
    for (int i = 0; i < n; i++) c[i] = a[i] - b[i];
}

static void vec_scale(double *v, double s, int n) {
    for (int i = 0; i < n; i++) v[i] *= s;
}

static void vec_axpy(double a, const double *x, double *y, int n) {
    for (int i = 0; i < n; i++) y[i] += a * x[i];
}

static double vec_dot(const double *a, const double *b, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

/* ──────────────────────── L1: IVP Construction ──────────────────────── */

ODEIVP ode_ivp_create(ODERHSFunc f, ODEJacobianFunc jac, void *ctx,
                       double t0, double t_end, const double *y0,
                       int dim, ODEClass ode_class) {
    ODEIVP ivp;
    ivp.f = f;
    ivp.jac = jac;
    ivp.event = NULL;
    ivp.ctx = ctx;
    ivp.t0 = t0;
    ivp.t_end = t_end;
    ivp.dim = dim;
    ivp.ode_class = ode_class;
    ivp.ode_order = (dim == 1) ? ODE_FIRST_ORDER : ODE_FIRST_ORDER_SYSTEM;

    /* Deep copy initial condition */
    ivp.y0 = (double *)malloc((size_t)dim * sizeof(double));
    if (ivp.y0 && y0) {
        vec_copy(y0, ivp.y0, dim);
    } else if (ivp.y0) {
        memset(ivp.y0, 0, (size_t)dim * sizeof(double));
    }
    return ivp;
}

bool ode_ivp_validate(const ODEIVP *ivp) {
    if (!ivp || !ivp->f || !ivp->y0) return false;
    if (ivp->dim <= 0) return false;

    /* Check time direction */
    if (ivp->t_end < ivp->t0) return false;

    /* Check initial condition for NaN/Inf */
    for (int i = 0; i < ivp->dim; i++) {
        if (!isfinite(ivp->y0[i])) return false;
    }

    /* Minimal functional test: evaluate at initial point */
    double *ftest = (double *)malloc((size_t)ivp->dim * sizeof(double));
    if (!ftest) return false;
    int ret = ivp->f(ivp->t0, ivp->y0, ftest, ivp->dim, ivp->ctx);
    if (ret != 0) { free(ftest); return false; }
    for (int i = 0; i < ivp->dim; i++) {
        if (!isfinite(ftest[i])) { free(ftest); return false; }
    }
    free(ftest);
    return true;
}

void ode_ivp_destroy(ODEIVP *ivp) {
    if (ivp && ivp->y0) {
        free(ivp->y0);
        ivp->y0 = NULL;
    }
}

ODEIVP ode_ivp_clone(const ODEIVP *ivp) {
    ODEIVP clone = *ivp;
    clone.y0 = (double *)malloc((size_t)ivp->dim * sizeof(double));
    if (clone.y0 && ivp->y0) {
        vec_copy(ivp->y0, clone.y0, ivp->dim);
    }
    return clone;
}

void ode_solution_free(ODESolution *sol) {
    if (sol) {
        free(sol->t_vals);
        free(sol->y_vals);
        sol->t_vals = NULL;
        sol->y_vals = NULL;
    }
}

/* ──────────────────────── L2: Lipschitz Condition ───────────────────── */

double ode_lipschitz_estimate(const ODEIVP *ivp, const double *y_center,
                                double radius, int n_samples) {
    if (!ivp || !y_center || n_samples <= 0 || ivp->dim <= 0) return 0.0;

    int dim = ivp->dim;
    double *y1 = (double *)malloc((size_t)dim * sizeof(double));
    double *y2 = (double *)malloc((size_t)dim * sizeof(double));
    double *f1 = (double *)malloc((size_t)dim * sizeof(double));
    double *f2 = (double *)malloc((size_t)dim * sizeof(double));
    double *diff_f = (double *)malloc((size_t)dim * sizeof(double));
    if (!y1 || !y2 || !f1 || !f2 || !diff_f) {
        free(y1); free(y2); free(f1); free(f2); free(diff_f);
        return 0.0;
    }

    double max_ratio = 0.0;

    for (int s = 0; s < n_samples; s++) {
        /* Generate two random points in the ball */
        for (int i = 0; i < dim; i++) {
            double r1 = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
            double r2 = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
            y1[i] = y_center[i] + r1 * radius;
            y2[i] = y_center[i] + r2 * radius;
        }

        if (ivp->f(ivp->t0, y1, f1, dim, ivp->ctx) != 0) continue;
        if (ivp->f(ivp->t0, y2, f2, dim, ivp->ctx) != 0) continue;

        vec_sub(f1, f2, diff_f, dim);
        double norm_diff_f = vec_norm_2(diff_f, dim);
        double norm_diff_y = vec_norm_2(y1, dim);
        double tmp_norm = vec_norm_2(y2, dim);
        vec_sub(y1, y2, diff_f, dim);
        double dist = vec_norm_2(diff_f, dim);

        double ratio = (dist > 1e-14) ? (norm_diff_f / dist) : 0.0;
        if (ratio > max_ratio) max_ratio = ratio;
    }

    free(y1); free(y2); free(f1); free(f2); free(diff_f);
    return max_ratio;
}

bool ode_one_sided_lipschitz(const ODEIVP *ivp, const double *y_center,
                               double radius, int n_samples, double *nu_out) {
    if (!ivp || !y_center || n_samples <= 0 || !nu_out) return false;

    int dim = ivp->dim;
    double *y1 = (double *)malloc((size_t)dim * sizeof(double));
    double *y2 = (double *)malloc((size_t)dim * sizeof(double));
    double *f1 = (double *)malloc((size_t)dim * sizeof(double));
    double *f2 = (double *)malloc((size_t)dim * sizeof(double));
    double *diff_y = (double *)malloc((size_t)dim * sizeof(double));
    double *diff_f = (double *)malloc((size_t)dim * sizeof(double));
    if (!y1 || !y2 || !f1 || !f2 || !diff_y || !diff_f) {
        free(y1); free(y2); free(f1); free(f2); free(diff_y); free(diff_f);
        return false;
    }

    double max_nu = -INFINITY;

    for (int s = 0; s < n_samples; s++) {
        for (int i = 0; i < dim; i++) {
            y1[i] = y_center[i] + ((double)rand() / RAND_MAX - 0.5) * 2.0 * radius;
            y2[i] = y_center[i] + ((double)rand() / RAND_MAX - 0.5) * 2.0 * radius;
        }

        if (ivp->f(ivp->t0, y1, f1, dim, ivp->ctx) != 0) continue;
        if (ivp->f(ivp->t0, y2, f2, dim, ivp->ctx) != 0) continue;

        vec_sub(y1, y2, diff_y, dim);
        vec_sub(f1, f2, diff_f, dim);
        double dot = vec_dot(diff_f, diff_y, dim);
        double norm_sq = vec_dot(diff_y, diff_y, dim);
        double nu = (norm_sq > 1e-14) ? (dot / norm_sq) : 0.0;
        if (nu > max_nu) max_nu = nu;
    }

    *nu_out = max_nu;
    free(y1); free(y2); free(f1); free(f2); free(diff_y); free(diff_f);
    return true;
}

/* ──────────────────────── L4: Picard Iteration ──────────────────────── */

double ode_picard_iteration(const ODEIVP *ivp,
                             const double *y_prev,
                             const double *t_vals,
                             int n_points,
                             double *y_next) {
    if (!ivp || !y_prev || !t_vals || !y_next || n_points < 2) return INFINITY;

    int dim = ivp->dim;
    double *f_val = (double *)malloc((size_t)dim * sizeof(double));
    double *integral = (double *)calloc((size_t)dim, sizeof(double));
    if (!f_val || !integral) { free(f_val); free(integral); return INFINITY; }

    /* y_{k+1}(t₀) = y₀ */
    vec_copy(ivp->y0, &y_next[0], dim);

    /* Trapezoidal integration of f(s, y_k(s)) from t₀ to each t_j */
    for (int j = 1; j < n_points; j++) {
        /* Integrate using trapezoidal rule step from t_{j-1} to t_j */
        int idx_prev = (j - 1) * dim;
        int idx_curr = j * dim;

        /* Evaluate f at previous point and current point of y_prev */
        double *f_prev = f_val;
        double *f_curr = (double *)malloc((size_t)dim * sizeof(double));
        if (!f_curr) { free(f_val); free(integral); return INFINITY; }

        ivp->f(t_vals[j-1], &y_prev[idx_prev], f_prev, dim, ivp->ctx);
        ivp->f(t_vals[j],   &y_prev[idx_curr], f_curr, dim, ivp->ctx);

        /* Trapezoidal: Δt/2 * (f_{j-1} + f_j) */
        double dt = t_vals[j] - t_vals[j-1];
        for (int i = 0; i < dim; i++) {
            integral[i] += 0.5 * dt * (f_prev[i] + f_curr[i]);
            y_next[idx_curr + i] = ivp->y0[i] + integral[i];
        }
        free(f_curr);
    }

    /* Compute change from y_prev */
    double max_change = 0.0;
    for (int j = 0; j < n_points * dim; j++) {
        double diff = fabs(y_next[j] - y_prev[j]);
        if (diff > max_change) max_change = diff;
    }

    free(f_val); free(integral);
    return max_change;
}

int ode_picard_solve(const ODEIVP *ivp, const double *t_vals, int n_points,
                      int max_iter, double rtol, double *y_out) {
    if (!ivp || !t_vals || !y_out || n_points < 2 || max_iter < 1) return 0;

    int dim = ivp->dim;
    int total = n_points * dim;

    /* Initialize with constant y₀ */
    double *y_prev = (double *)malloc((size_t)total * sizeof(double));
    double *y_next = (double *)malloc((size_t)total * sizeof(double));
    if (!y_prev || !y_next) { free(y_prev); free(y_next); return 0; }

    for (int j = 0; j < n_points; j++) {
        vec_copy(ivp->y0, &y_prev[j * dim], dim);
    }

    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        double change = ode_picard_iteration(ivp, y_prev, t_vals,
                                               n_points, y_next);
        if (!isfinite(change)) { free(y_prev); free(y_next); return 0; }

        /* Check convergence */
        double max_y = 0.0;
        for (int i = 0; i < total; i++) {
            if (fabs(y_next[i]) > max_y) max_y = fabs(y_next[i]);
        }
        double atol = rtol * (max_y > 1.0 ? max_y : 1.0);
        if (change <= atol) {
            vec_copy(y_next, y_out, total);
            free(y_prev); free(y_next);
            return iter + 1;
        }

        /* Swap for next iteration */
        double *tmp = y_prev; y_prev = y_next; y_next = tmp;
    }

    /* Max iterations reached */
    vec_copy(y_prev, y_out, total);
    free(y_prev); free(y_next);
    return 0; /* 0 = didn't converge */
}

/* ──────────────────────── L4: Gronwall's Inequality ──────────────────── */

int ode_gronwall_bound(double u0, const double *beta_vals,
                        const double *t_vals, int n_points,
                        double *u_bound) {
    if (!beta_vals || !t_vals || !u_bound || n_points < 1) return -1;
    if (u0 < 0.0) return -1;

    /* u(t) ≤ u(0)·exp(∫₀ᵗ β(s) ds) */
    double integral = 0.0;
    u_bound[0] = u0;

    for (int j = 1; j < n_points; j++) {
        double dt = t_vals[j] - t_vals[j-1];
        /* Trapezoidal integration of β */
        integral += 0.5 * dt * (beta_vals[j-1] + beta_vals[j]);
        /* u_bound[j] = u0 * exp(integral) */
        u_bound[j] = u0 * exp(integral);
    }
    return 0;
}

double ode_continuous_dependence_bound(const ODEIVP *ivp, double L,
                                         const double *y0_a,
                                         const double *y0_b, double t) {
    if (!ivp || !y0_a || !y0_b || L < 0.0) return INFINITY;

    int dim = ivp->dim;
    double *diff = (double *)malloc((size_t)dim * sizeof(double));
    if (!diff) return INFINITY;

    vec_sub(y0_a, y0_b, diff, dim);
    double init_dist = vec_norm_2(diff, dim);
    double bound = init_dist * exp(L * fabs(t - ivp->t0));
    free(diff);
    return bound;
}

/* ──────────────────────── Default Config ────────────────────────────── */

ODESimConfig ode_simconfig_default(const ODEIVP *ivp) {
    ODESimConfig cfg;
    cfg.method = METHOD_RK45_DORMAND_PRINCE;
    cfg.h_init = 0.01 * (ivp->t_end - ivp->t0);
    if (cfg.h_init <= 0.0) cfg.h_init = 1e-6;
    cfg.h_min = 1e-14;
    cfg.h_max = 0.1 * (ivp->t_end - ivp->t0);
    if (cfg.h_max <= cfg.h_init) cfg.h_max = cfg.h_init * 10.0;
    cfg.rtol = 1e-6;
    cfg.atol = 1e-8;
    cfg.max_steps = 1000000;
    cfg.dense_output = false;
    cfg.dense_dt = cfg.h_init;
    return cfg;
}

/* ──────────────────────── Solution Utilities ────────────────────────── */

void ode_solution_print(const ODESolution *sol) {
    if (!sol) { printf("ODESolution: NULL\n"); return; }
    printf("ODESolution: n_points=%d  dim=%d  steps=%d  rejected=%d\n",
           sol->n_points, sol->dim, sol->total_steps, sol->rejected_steps);
    printf("  fevals=%d  final_err=%.2e  success=%s\n",
           sol->total_fevals, sol->final_err, sol->success ? "yes" : "no");
    printf("  message: %s\n", sol->message);
}

int ode_solution_interpolate(const ODESolution *sol, double t, double *y_out) {
    if (!sol || !y_out || sol->n_points < 2) return -1;
    if (t < sol->t_vals[0] - 1e-12 || t > sol->t_vals[sol->n_points-1] + 1e-12)
        return -1;

    /* Find interval via binary search */
    int lo = 0, hi = sol->n_points - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (sol->t_vals[mid] <= t) lo = mid;
        else hi = mid;
    }

    int i0 = lo, i1 = hi;
    double t0 = sol->t_vals[i0], t1 = sol->t_vals[i1];
    double alpha = (t1 - t0 > 1e-14) ? (t - t0) / (t1 - t0) : 0.0;
    alpha = alpha < 0.0 ? 0.0 : (alpha > 1.0 ? 1.0 : alpha);

    int dim = sol->dim;
    for (int d = 0; d < dim; d++) {
        double y0 = sol->y_vals[i0 * dim + d];
        double y1 = sol->y_vals[i1 * dim + d];
        y_out[d] = y0 + alpha * (y1 - y0);
    }
    return 0;
}

double ode_solution_error_norm(const ODESolution *computed,
                                 const ODESolution *reference) {
    if (!computed || !reference || computed->dim != reference->dim) return INFINITY;
    if (computed->n_points < 2 || reference->n_points < 2) return INFINITY;

    int dim = computed->dim;
    double max_err = 0.0;

    /* Compare at computed time points via interpolation of reference */
    for (int j = 0; j < computed->n_points; j++) {
        double *y_ref = (double *)malloc((size_t)dim * sizeof(double));
        if (!y_ref) continue;
        if (ode_solution_interpolate(reference, computed->t_vals[j], y_ref) == 0) {
            for (int d = 0; d < dim; d++) {
                double err = fabs(computed->y_vals[j * dim + d] - y_ref[d]);
                if (err > max_err) max_err = err;
            }
        }
        free(y_ref);
    }
    return max_err;
}
