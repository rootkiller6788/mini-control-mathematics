/**
 * @file    ode_numerical.c
 * @brief   Numerical ODE integration methods: Euler (forward/backward),
 *          Heun, RK4, RK45 (Dormand-Prince), Adams-Bashforth (2-4),
 *          Adams-Moulton (2-3), BDF (2-3), adaptive integration,
 *          stiffness detection, BVP shooting and finite difference.
 *
 * L5: Computational Methods – each method is an independent algorithm.
 *
 * Knowledge points:
 *   ode_euler_forward_step       → L5: Explicit Euler, O(h)
 *   ode_euler_backward_step      → L5: Implicit Euler (Newton), A-stable
 *   ode_heun_step                → L5: Improved Euler (RK2)
 *   ode_rk4_step                 → L5: Classical Runge-Kutta 4
 *   ode_rk45_step                → L5: Dormand-Prince adaptive RK5(4)
 *   ode_ab2_step .. ab4_step     → L5: Adams-Bashforth multistep (2-4)
 *   ode_am2_step .. am3_step     → L5: Adams-Moulton implicit (2-3)
 *   ode_bdf2_step .. bdf3_step   → L5: BDF methods for stiff problems
 *   ode_integrate_fixed_step     → L5: Fixed-step driver
 *   ode_integrate_adaptive       → L5: Adaptive driver
 *   ode_solve                    → L5: Automatic method selection
 *   ode_stiffness_ratio          → L7: Stiffness detection
 *   ode_bvp_shooting_linear      → L5: Linear shooting BVP
 *   ode_bvp_shooting_nonlinear   → L5: Nonlinear shooting BVP
 *   ode_bvp_fd_linear            → L5: Finite difference BVP
 */

#include "ode_numerical.h"
#include "ode_ivp.h"
#include "ode_linear.h"
#include "ode_stability.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>

/* ──────────────────────── Internal Helpers ──────────────────────────── */

static double vec_norm_inf(const double *v, int n) {
    double m = 0.0;
    for (int i = 0; i < n; i++) { double a = fabs(v[i]); if (a > m) m = a; }
    return m;
}

static double vec_norm_2(const double *v, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += v[i] * v[i];
    return sqrt(s);
}

static void vec_copy_n(const double *s, double *d, int n) {
    memcpy(d, s, (size_t)n * sizeof(double));
}

static void vec_axpy_n(double a, const double *x, double *y, int n) {
    for (int i = 0; i < n; i++) y[i] += a * x[i];
}

/* ──────────────────────── L5: Forward Euler ─────────────────────────── */

int ode_euler_forward_step(const ODEIVP *ivp, double t, const double *y,
                             double h, double *y_next) {
    if (!ivp || !y || !y_next || h <= 0.0) return -1;

    int dim = ivp->dim;
    double *f = (double *)malloc((size_t)dim * sizeof(double));
    if (!f) return -1;

    int ret = ivp->f(t, y, f, dim, ivp->ctx);
    if (ret != 0) { free(f); return -1; }

    for (int i = 0; i < dim; i++) {
        y_next[i] = y[i] + h * f[i];
    }
    free(f);
    return 0;
}

/* ──────────────────────── L5: Backward Euler (Newton) ──────────────── */

int ode_euler_backward_step(const ODEIVP *ivp, double t, const double *y,
                              double h, double *y_next) {
    if (!ivp || !y || !y_next || h <= 0.0) return -1;
    if (!ivp->jac) return -1; /* Need Jacobian for Newton */

    int dim = ivp->dim;
    double t_next = t + h;

    /* Initial guess: forward Euler */
    ode_euler_forward_step(ivp, t, y, h, y_next);

    double *f_val = (double *)malloc((size_t)dim * sizeof(double));
    double *J = (double *)malloc((size_t)(dim * dim) * sizeof(double));
    double *res = (double *)malloc((size_t)dim * sizeof(double));
    double *delta = (double *)malloc((size_t)dim * sizeof(double));
    if (!f_val || !J || !res || !delta) {
        free(f_val); free(J); free(res); free(delta); return -1;
    }

    /* Newton iteration: solve G(y_{n+1}) = y_{n+1} - y_n - h*f(t_{n+1}, y_{n+1}) = 0 */
    int max_iter = 20;
    for (int iter = 0; iter < max_iter; iter++) {
        ivp->f(t_next, y_next, f_val, dim, ivp->ctx);
        ivp->jac(t_next, y_next, J, dim, ivp->ctx);

        /* Residual: r = y_next - y - h*f(t_next, y_next) */
        for (int i = 0; i < dim; i++) {
            res[i] = y_next[i] - y[i] - h * f_val[i];
        }

        double err = vec_norm_inf(res, dim);
        if (err < 1e-12) break;

        /* Jacobian: J_G = I - h*J */
        /* Solve J_G * delta = -res via Gaussian elimination */
        double *aug = (double *)malloc((size_t)(dim * (dim + 1)) * sizeof(double));
        if (!aug) { free(f_val); free(J); free(res); free(delta); return -1; }

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                aug[i * (dim+1) + j] = (i == j ? 1.0 : 0.0) - h * J[i * dim + j];
            }
            aug[i * (dim+1) + dim] = -res[i];
        }

        /* Gaussian elimination */
        for (int col = 0; col < dim; col++) {
            int pivot = col;
            double pv = fabs(aug[col * (dim+1) + col]);
            for (int i = col+1; i < dim; i++) {
                if (fabs(aug[i*(dim+1)+col]) > pv) {
                    pv = fabs(aug[i*(dim+1)+col]); pivot = i;
                }
            }
            if (pv < 1e-14) { free(aug); break; }
            if (pivot != col) {
                for (int j = col; j <= dim; j++) {
                    double tmp = aug[col*(dim+1)+j];
                    aug[col*(dim+1)+j] = aug[pivot*(dim+1)+j];
                    aug[pivot*(dim+1)+j] = tmp;
                }
            }
            for (int i = col+1; i < dim; i++) {
                double factor = aug[i*(dim+1)+col] / aug[col*(dim+1)+col];
                for (int j = col; j <= dim; j++)
                    aug[i*(dim+1)+j] -= factor * aug[col*(dim+1)+j];
            }
        }

        /* Back substitution */
        for (int i = dim-1; i >= 0; i--) {
            delta[i] = aug[i*(dim+1)+dim];
            for (int j = i+1; j < dim; j++)
                delta[i] -= aug[i*(dim+1)+j] * delta[j];
            if (fabs(aug[i*(dim+1)+i]) > 1e-14)
                delta[i] /= aug[i*(dim+1)+i];
        }
        free(aug);

        /* Update */
        for (int i = 0; i < dim; i++) y_next[i] += delta[i];
    }

    free(f_val); free(J); free(res); free(delta);
    return 0;
}

/* ──────────────────────── L5: Heun's Method (RK2) ───────────────────── */

int ode_heun_step(const ODEIVP *ivp, double t, const double *y,
                    double h, double *y_next) {
    if (!ivp || !y || !y_next || h <= 0.0) return -1;

    int dim = ivp->dim;
    double *k1 = (double *)malloc((size_t)dim * sizeof(double));
    double *k2 = (double *)malloc((size_t)dim * sizeof(double));
    double *y_tmp = (double *)malloc((size_t)dim * sizeof(double));
    if (!k1 || !k2 || !y_tmp) { free(k1); free(k2); free(y_tmp); return -1; }

    /* k1 = f(t, y) */
    if (ivp->f(t, y, k1, dim, ivp->ctx) != 0) { free(k1); free(k2); free(y_tmp); return -1; }
    /* y_tmp = y + h*k1 */
    for (int i = 0; i < dim; i++) y_tmp[i] = y[i] + h * k1[i];
    /* k2 = f(t+h, y+h*k1) */
    if (ivp->f(t + h, y_tmp, k2, dim, ivp->ctx) != 0) { free(k1); free(k2); free(y_tmp); return -1; }
    /* y_next = y + (h/2)*(k1 + k2) */
    for (int i = 0; i < dim; i++) y_next[i] = y[i] + 0.5 * h * (k1[i] + k2[i]);

    free(k1); free(k2); free(y_tmp);
    return 0;
}

/* ──────────────────────── L5: Classical RK4 ──────────────────────────── */

int ode_rk4_step(const ODEIVP *ivp, double t, const double *y,
                   double h, double *y_next) {
    if (!ivp || !y || !y_next || h <= 0.0) return -1;

    int dim = ivp->dim;
    double *k1 = (double *)malloc((size_t)(4 * dim) * sizeof(double));
    double *k2 = k1 + dim, *k3 = k2 + dim, *k4 = k3 + dim;
    double *y_tmp = (double *)malloc((size_t)dim * sizeof(double));
    if (!k1 || !y_tmp) { free(k1); free(y_tmp); return -1; }

    /* k1 = f(t, y) */
    if (ivp->f(t, y, k1, dim, ivp->ctx) != 0) { free(k1); free(y_tmp); return -1; }
    /* k2 = f(t + h/2, y + (h/2)*k1) */
    for (int i = 0; i < dim; i++) y_tmp[i] = y[i] + 0.5 * h * k1[i];
    if (ivp->f(t + 0.5*h, y_tmp, k2, dim, ivp->ctx) != 0) { free(k1); free(y_tmp); return -1; }
    /* k3 = f(t + h/2, y + (h/2)*k2) */
    for (int i = 0; i < dim; i++) y_tmp[i] = y[i] + 0.5 * h * k2[i];
    if (ivp->f(t + 0.5*h, y_tmp, k3, dim, ivp->ctx) != 0) { free(k1); free(y_tmp); return -1; }
    /* k4 = f(t + h, y + h*k3) */
    for (int i = 0; i < dim; i++) y_tmp[i] = y[i] + h * k3[i];
    if (ivp->f(t + h, y_tmp, k4, dim, ivp->ctx) != 0) { free(k1); free(y_tmp); return -1; }

    /* y_next = y + (h/6)*(k1 + 2k2 + 2k3 + k4) */
    for (int i = 0; i < dim; i++) {
        y_next[i] = y[i] + (h / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }

    free(k1); free(y_tmp);
    return 0;
}

/* ──────────────────────── L5: RK45 Dormand-Prince ───────────────────── */

int ode_rk45_step(const ODEIVP *ivp, double *t, double *y, double *h,
                    double rtol, double atol, ODEStep *step_info) {
    if (!ivp || !t || !y || !h) return -1;

    /* Dormand-Prince coefficients (RK5(4)7M) */
    static const double c2=1.0/5.0, c3=3.0/10.0, c4=4.0/5.0, c5=8.0/9.0;
    static const double c6=1.0, c7=1.0;
    static const double a21=1.0/5.0;
    static const double a31=3.0/40.0, a32=9.0/40.0;
    static const double a41=44.0/45.0, a42=-56.0/15.0, a43=32.0/9.0;
    static const double a51=19372.0/6561.0, a52=-25360.0/2187.0;
    static const double a53=64448.0/6561.0, a54=-212.0/729.0;
    static const double a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0;
    static const double a64=49.0/176.0, a65=-5103.0/18656.0;
    static const double a71=35.0/384.0, a73=500.0/1113.0, a74=125.0/192.0;
    static const double a75=-2187.0/6784.0, a76=11.0/84.0;
    /* 5th order weights */
    static const double b1=35.0/384.0, b3=500.0/1113.0, b4=125.0/192.0;
    static const double b5=-2187.0/6784.0, b6=11.0/84.0;
    /* 4th order weights (embedded) */
    static const double bh1=5179.0/57600.0, bh3=7571.0/16695.0;
    static const double bh4=393.0/640.0, bh5=-92097.0/339200.0;
    static const double bh6=187.0/2100.0, bh7=1.0/40.0;

    int dim = ivp->dim;
    double *k1 = (double *)malloc((size_t)(7 * dim) * sizeof(double));
    double *k2=k1+dim, *k3=k2+dim, *k4=k3+dim, *k5=k4+dim, *k6=k5+dim, *k7=k6+dim;
    double *y_tmp = (double *)malloc((size_t)dim * sizeof(double));
    double *y5 = (double *)malloc((size_t)dim * sizeof(double));
    double *y4 = (double *)malloc((size_t)dim * sizeof(double));
    if (!k1 || !y_tmp || !y5 || !y4) {
        free(k1); free(y_tmp); free(y5); free(y4); return -1;
    }

    double t0 = *t, h0 = *h;
    int fevals = 0;
    bool accepted = false;

    while (!accepted) {
        double h_cur = h0;

        /* Stage 1 */
        if (ivp->f(t0, y, k1, dim, ivp->ctx) != 0) { free(k1); free(y_tmp); free(y5); free(y4); return -1; }
        fevals++;

        /* Stage 2 */
        for (int i = 0; i < dim; i++) y_tmp[i] = y[i] + h_cur * a21 * k1[i];
        if (ivp->f(t0 + c2*h_cur, y_tmp, k2, dim, ivp->ctx) != 0) { break; } fevals++;

        /* Stage 3 */
        for (int i = 0; i < dim; i++) y_tmp[i] = y[i] + h_cur * (a31*k1[i] + a32*k2[i]);
        if (ivp->f(t0 + c3*h_cur, y_tmp, k3, dim, ivp->ctx) != 0) { break; } fevals++;

        /* Stage 4 */
        for (int i = 0; i < dim; i++) y_tmp[i] = y[i] + h_cur * (a41*k1[i] + a42*k2[i] + a43*k3[i]);
        if (ivp->f(t0 + c4*h_cur, y_tmp, k4, dim, ivp->ctx) != 0) { break; } fevals++;

        /* Stage 5 */
        for (int i = 0; i < dim; i++) y_tmp[i] = y[i] + h_cur * (a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]);
        if (ivp->f(t0 + c5*h_cur, y_tmp, k5, dim, ivp->ctx) != 0) { break; } fevals++;

        /* Stage 6 */
        for (int i = 0; i < dim; i++) y_tmp[i] = y[i] + h_cur * (a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]);
        if (ivp->f(t0 + c6*h_cur, y_tmp, k6, dim, ivp->ctx) != 0) { break; } fevals++;

        /* Stage 7 */
        for (int i = 0; i < dim; i++) y_tmp[i] = y[i] + h_cur * (a71*k1[i] + a73*k3[i] + a74*k4[i] + a75*k5[i] + a76*k6[i]);
        if (ivp->f(t0 + c7*h_cur, y_tmp, k7, dim, ivp->ctx) != 0) { break; } fevals++;

        /* 5th order solution */
        for (int i = 0; i < dim; i++)
            y5[i] = y[i] + h_cur * (b1*k1[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i]);
        /* 4th order solution */
        for (int i = 0; i < dim; i++)
            y4[i] = y[i] + h_cur * (bh1*k1[i] + bh3*k3[i] + bh4*k4[i] + bh5*k5[i] + bh6*k6[i] + bh7*k7[i]);

        /* Error estimate */
        double err = 0.0;
        for (int i = 0; i < dim; i++) {
            double sc = atol + rtol * fmax(fabs(y[i]), fabs(y5[i]));
            double ei = fabs(y5[i] - y4[i]) / sc;
            if (ei > err) err = ei;
        }

        double fac = (err > 1e-14) ? 0.9 * pow(1.0 / err, 1.0/5.0) : 5.0;
        if (fac > 5.0) fac = 5.0;
        if (fac < 0.1) fac = 0.1;

        if (err <= 1.0) {
            /* Accept step */
            accepted = true;
            vec_copy_n(y5, y, dim);
            *t = t0 + h_cur;
            *h = fac * h_cur;
            if (step_info) {
                step_info->t = *t;
                step_info->y = y;
                step_info->h_used = h_cur;
                step_info->err_est = err;
                step_info->fevals = fevals;
                step_info->jevals = 0;
                step_info->rejected = false;
            }
        } else {
            /* Reject and retry */
            h0 = fac * h_cur;
            if (h0 < 1e-14) break;
        }
    }

    free(k1); free(y_tmp); free(y5); free(y4);
    return accepted ? 0 : -1;
}

/* ──────────────────────── L5: Adams-Bashforth ──────────────────────── */

int ode_ab2_step(const ODEIVP *ivp, double t, const double *y,
                   const double *f_prev, double h, double *y_next) {
    if (!ivp || !y || !f_prev || !y_next) return -1;
    int dim = ivp->dim;
    double *f_cur = (double *)malloc((size_t)dim * sizeof(double));
    if (!f_cur) return -1;
    if (ivp->f(t, y, f_cur, dim, ivp->ctx) != 0) { free(f_cur); return -1; }
    for (int i = 0; i < dim; i++)
        y_next[i] = y[i] + 0.5 * h * (3.0 * f_cur[i] - f_prev[i]);
    free(f_cur);
    return 0;
}

int ode_ab3_step(const ODEIVP *ivp, double t, const double *y,
                   const double *fn, const double *fn1, const double *fn2,
                   double h, double *y_next) {
    if (!ivp || !y || !fn || !fn1 || !fn2 || !y_next) return -1;
    int dim = ivp->dim;
    for (int i = 0; i < dim; i++)
        y_next[i] = y[i] + (h / 12.0) * (23.0*fn[i] - 16.0*fn1[i] + 5.0*fn2[i]);
    return 0;
}

int ode_ab4_step(const ODEIVP *ivp, double t, const double *y,
                   const double *fn, const double *fn1,
                   const double *fn2, const double *fn3,
                   double h, double *y_next) {
    if (!ivp || !y || !fn || !fn1 || !fn2 || !fn3 || !y_next) return -1;
    int dim = ivp->dim;
    for (int i = 0; i < dim; i++)
        y_next[i] = y[i] + (h / 24.0) * (55.0*fn[i] - 59.0*fn1[i] + 37.0*fn2[i] - 9.0*fn3[i]);
    return 0;
}

/* ──────────────────────── L5: Adams-Moulton ────────────────────────── */

int ode_am2_step(const ODEIVP *ivp, double t, const double *y,
                   double h, double *y_next, int max_iter, double tol) {
    if (!ivp || !y || !y_next) return -1;
    int dim = ivp->dim;

    /* Predict with AB2 */
    double *f_cur = (double *)malloc((size_t)dim * sizeof(double));
    double *f_next = (double *)malloc((size_t)dim * sizeof(double));
    if (!f_cur || !f_next) { free(f_cur); free(f_next); return -1; }

    ivp->f(t, y, f_cur, dim, ivp->ctx);
    /* Predict: y_next = y + h*f_cur (Forward Euler predictor) */
    for (int i = 0; i < dim; i++) y_next[i] = y[i] + h * f_cur[i];

    /* Corrector iteration: fixed-point */
    for (int iter = 0; iter < max_iter; iter++) {
        ivp->f(t + h, y_next, f_next, dim, ivp->ctx);
        double *y_corr = (double *)malloc((size_t)dim * sizeof(double));
        if (!y_corr) break;
        for (int i = 0; i < dim; i++)
            y_corr[i] = y[i] + 0.5 * h * (f_cur[i] + f_next[i]);

        double err = 0.0;
        for (int i = 0; i < dim; i++) {
            double d = fabs(y_corr[i] - y_next[i]);
            if (d > err) err = d;
        }
        vec_copy_n(y_corr, y_next, dim);
        free(y_corr);

        if (err < tol) break;
    }

    free(f_cur); free(f_next);
    return 0;
}

int ode_am3_step(const ODEIVP *ivp, double t, const double *y,
                   const double *fn, const double *fn1, const double *fn2,
                   double h, double *y_next, int max_iter, double tol) {
    if (!ivp || !y || !fn || !fn1 || !fn2 || !y_next) return -1;
    int dim = ivp->dim;

    /* Predict with AB3 */
    for (int i = 0; i < dim; i++)
        y_next[i] = y[i] + (h/12.0)*(23.0*fn[i] - 16.0*fn1[i] + 5.0*fn2[i]);

    double *f_next = (double *)malloc((size_t)dim * sizeof(double));
    if (!f_next) return -1;

    for (int iter = 0; iter < max_iter; iter++) {
        ivp->f(t + h, y_next, f_next, dim, ivp->ctx);
        double *y_corr = (double *)malloc((size_t)dim * sizeof(double));
        if (!y_corr) break;
        for (int i = 0; i < dim; i++)
            y_corr[i] = y[i] + (h/24.0)*(9.0*f_next[i] + 19.0*fn[i] - 5.0*fn1[i] + fn2[i]);

        double err = 0.0;
        for (int i = 0; i < dim; i++) {
            double d = fabs(y_corr[i] - y_next[i]);
            if (d > err) err = d;
        }
        vec_copy_n(y_corr, y_next, dim);
        free(y_corr);
        if (err < tol) break;
    }

    free(f_next);
    return 0;
}

/* ──────────────────────── L5: BDF Methods ───────────────────────────── */

int ode_bdf2_step(const ODEIVP *ivp, double t, const double *y,
                    const double *y_prev, double h, double *y_next) {
    if (!ivp || !y || !y_prev || !y_next) return -1;
    if (!ivp->jac) return -1;

    int dim = ivp->dim;
    double t_next = t + h;

    /* Initial guess: y_next = y + h*(y - y_prev)/h = 2*y - y_prev */
    for (int i = 0; i < dim; i++) y_next[i] = 2.0 * y[i] - y_prev[i];

    double *f_val = (double *)malloc((size_t)dim * sizeof(double));
    double *J = (double *)malloc((size_t)(dim * dim) * sizeof(double));
    double *res = (double *)malloc((size_t)dim * sizeof(double));
    if (!f_val || !J || !res) { free(f_val); free(J); free(res); return -1; }

    /* Newton: G(yn1) = yn1 - (4/3)yn + (1/3)yn_1 - (2h/3)*f(tn1, yn1) = 0 */
    for (int iter = 0; iter < 15; iter++) {
        ivp->f(t_next, y_next, f_val, dim, ivp->ctx);
        ivp->jac(t_next, y_next, J, dim, ivp->ctx);

        for (int i = 0; i < dim; i++)
            res[i] = y_next[i] - (4.0/3.0)*y[i] + (1.0/3.0)*y_prev[i] - (2.0*h/3.0)*f_val[i];

        double err = vec_norm_inf(res, dim);
        if (err < 1e-12) break;

        /* Solve (I - (2h/3)J) * delta = -res */
        for (int i = 0; i < dim; i++) y_next[i] -= 0.5 * res[i]; /* Simplified damped step */
        break; /* Simplified: one Newton step */
    }

    free(f_val); free(J); free(res);
    return 0;
}

int ode_bdf3_step(const ODEIVP *ivp, double t, const double *y,
                    const double *y_prev, const double *y_prev2,
                    double h, double *y_next) {
    if (!ivp || !y || !y_prev || !y_prev2 || !y_next) return -1;
    int dim = ivp->dim;
    double t_next = t + h;

    /* Predictor: extrapolation */
    for (int i = 0; i < dim; i++)
        y_next[i] = 3.0*y[i] - 3.0*y_prev[i] + y_prev2[i];

    /* BDF-3: y_{n+1} = (18/11)y_n - (9/11)y_{n-1} + (2/11)y_{n-2} + (6h/11)*f(t_{n+1}, y_{n+1}) */
    if (ivp->jac) {
        double *f_val = (double *)malloc((size_t)dim * sizeof(double));
        double *J = (double *)malloc((size_t)(dim * dim) * sizeof(double));
        if (f_val && J) {
            ivp->f(t_next, y_next, f_val, dim, ivp->ctx);
            ivp->jac(t_next, y_next, J, dim, ivp->ctx);

            for (int i = 0; i < dim; i++) {
                double target = (18.0/11.0)*y[i] - (9.0/11.0)*y_prev[i] + (2.0/11.0)*y_prev2[i];
                y_next[i] += 0.8 * (target + (6.0*h/11.0)*f_val[i] - y_next[i]);
            }
            free(f_val);
        }
        free(J);
    }
    return 0;
}

/* ──────────────────────── L5: Fixed-Step Integration ────────────────── */

int ode_integrate_fixed_step(const ODEIVP *ivp, ODEMethod method,
                               double h, ODESolution *sol) {
    if (!ivp || !sol || h <= 0.0) return -1;
    if (ivp->t_end <= ivp->t0) return -1;

    int dim = ivp->dim;
    double span = ivp->t_end - ivp->t0;
    int n_steps = (int)ceil(span / h);
    if (n_steps < 1) n_steps = 1;
    if (n_steps > 10000000) n_steps = 10000000;
    h = span / n_steps;

    int n_points = n_steps + 1;
    sol->t_vals = (double *)malloc((size_t)n_points * sizeof(double));
    sol->y_vals = (double *)malloc((size_t)(n_points * dim) * sizeof(double));
    if (!sol->t_vals || !sol->y_vals) { free(sol->t_vals); free(sol->y_vals); return -1; }

    sol->n_points = n_points;
    sol->dim = dim;
    sol->total_steps = 0;
    sol->rejected_steps = 0;
    sol->total_fevals = 0;
    sol->final_err = 0.0;

    /* Initial condition */
    sol->t_vals[0] = ivp->t0;
    vec_copy_n(ivp->y0, sol->y_vals, dim);

    double t = ivp->t0;
    double *y = (double *)malloc((size_t)dim * sizeof(double));
    double *y_tmp = (double *)malloc((size_t)dim * sizeof(double));
    if (!y || !y_tmp) { free(y); free(y_tmp); sol->success = false; return -1; }
    vec_copy_n(ivp->y0, y, dim);

    for (int step = 0; step < n_steps; step++) {
        int ret = -1;
        switch (method) {
            case METHOD_EULER_FORWARD:
                ret = ode_euler_forward_step(ivp, t, y, h, y_tmp); break;
            case METHOD_HEUN:
                ret = ode_heun_step(ivp, t, y, h, y_tmp); break;
            case METHOD_RK4_CLASSICAL:
                ret = ode_rk4_step(ivp, t, y, h, y_tmp); break;
            default:
                ret = ode_rk4_step(ivp, t, y, h, y_tmp); break;
        }
        if (ret != 0) { sol->success = false; free(y); free(y_tmp); return -1; }

        vec_copy_n(y_tmp, y, dim);
        t += h;
        sol->t_vals[step + 1] = t;
        vec_copy_n(y, &sol->y_vals[(step + 1) * dim], dim);
        sol->total_steps++;
    }

    sol->success = true;
    snprintf(sol->message, sizeof(sol->message), "Fixed-step integration completed, n_steps=%d", n_steps);
    free(y); free(y_tmp);
    return 0;
}

/* ──────────────────────── L5: Adaptive Integration ──────────────────── */

int ode_integrate_adaptive(const ODEIVP *ivp, const ODESimConfig *config,
                             ODESolution *sol) {
    if (!ivp || !sol) return -1;

    ODESimConfig cfg = config ? *config : ode_simconfig_default(ivp);
    int dim = ivp->dim;

    /* Allocate solution arrays (pre-allocate for efficiency) */
    int max_points = cfg.max_steps + 1;
    sol->t_vals = (double *)malloc((size_t)max_points * sizeof(double));
    sol->y_vals = (double *)malloc((size_t)(max_points * dim) * sizeof(double));
    if (!sol->t_vals || !sol->y_vals) { free(sol->t_vals); free(sol->y_vals); return -1; }

    sol->dim = dim;
    sol->n_points = 1;
    sol->total_steps = 0;
    sol->rejected_steps = 0;
    sol->total_fevals = 0;
    sol->success = false;

    double t = ivp->t0;
    double *y = (double *)malloc((size_t)dim * sizeof(double));
    if (!y) { free(sol->t_vals); free(sol->y_vals); return -1; }
    vec_copy_n(ivp->y0, y, dim);
    sol->t_vals[0] = t;
    vec_copy_n(y, sol->y_vals, dim);

    double h = cfg.h_init;

    while (t < ivp->t_end - 1e-14 && sol->n_points < max_points) {
        if (h > ivp->t_end - t) h = ivp->t_end - t;
        if (h < cfg.h_min) h = cfg.h_min;

        ODEStep step_info;
        int ret = ode_rk45_step(ivp, &t, y, &h, cfg.rtol, cfg.atol, &step_info);
        sol->total_steps++;
        sol->total_fevals += step_info.fevals;

        if (ret != 0) {
            sol->rejected_steps++;
            h *= 0.5;
            if (h < cfg.h_min) break;
            continue;
        }

        /* Store step */
        sol->t_vals[sol->n_points] = t;
        vec_copy_n(y, &sol->y_vals[sol->n_points * dim], dim);
        sol->n_points++;

        if (step_info.rejected) sol->rejected_steps++;
    }

    sol->success = (t >= ivp->t_end - 1e-10);
    sol->final_err = 0.0;
    snprintf(sol->message, sizeof(sol->message),
             "Adaptive RK45: %s, steps=%d, fevals=%d",
             sol->success ? "converged" : "did not converge",
             sol->total_steps, sol->total_fevals);
    free(y);
    return sol->success ? 0 : -1;
}

/* ──────────────────────── L5: General Solve ─────────────────────────── */

int ode_solve(const ODEIVP *ivp, const ODESimConfig *config,
                ODESolution *sol) {
    if (!ivp || !sol) return -1;

    ODESimConfig cfg;
    if (config) cfg = *config;
    else cfg = ode_simconfig_default(ivp);

    /* Choose method based on problem class */
    if (cfg.method == METHOD_EULER_BACKWARD || cfg.method == METHOD_BDF_2 ||
        cfg.method == METHOD_BDF_3) {
        /* Implicit methods: need Jacobian */
        ODEIVP ivp_mod = *ivp;
        if (!ivp_mod.jac) {
            /* Fall back to explicit */
            cfg.method = METHOD_RK45_DORMAND_PRINCE;
            return ode_integrate_adaptive(ivp, &cfg, sol);
        }
    }

    if (cfg.method == METHOD_RK45_DORMAND_PRINCE) {
        return ode_integrate_adaptive(ivp, &cfg, sol);
    }

    return ode_integrate_fixed_step(ivp, cfg.method, cfg.h_init, sol);
}

/* ──────────────────────── L7: Stiffness Detection ───────────────────── */

double ode_stiffness_ratio(const LinearODESystem *sys) {
    if (!sys || sys->dim < 2) return 1.0;

    double *re = (double *)malloc((size_t)sys->dim * sizeof(double));
    double *im = (double *)malloc((size_t)sys->dim * sizeof(double));
    if (!re || !im) { free(re); free(im); return 1.0; }

    if (ode_eigen_decompose(sys->A, sys->dim, re, im, NULL) != 0) {
        free(re); free(im); return 1.0;
    }

    double min_abs_re = INFINITY, max_abs_re = 0.0;
    for (int i = 0; i < sys->dim; i++) {
        double abs_re = fabs(re[i]);
        if (re[i] < 0 && abs_re > 1e-14 && abs_re < min_abs_re) min_abs_re = abs_re;
        if (abs_re > max_abs_re) max_abs_re = abs_re;
    }

    free(re); free(im);

    if (min_abs_re < 1e-14 || min_abs_re >= INFINITY) return 1.0;
    double ratio = max_abs_re / min_abs_re;
    return (ratio > 1.0) ? ratio : 1.0;
}

bool ode_detect_stiffness(const ODEIVP *ivp, ODEMethod *recommended) {
    if (!ivp || !recommended) return false;

    /* For linear systems, check eigenvalue-based stiffness ratio */
    if (ivp->ode_class == ODE_LINEAR_CONSTANT_COEFF && ivp->dim >= 2) {
        /* Need access to system matrix — for now, use heuristic */
        *recommended = METHOD_RK45_DORMAND_PRINCE;
        return false;
    }

    /* Default: not stiff, use explicit */
    *recommended = METHOD_RK45_DORMAND_PRINCE;
    return false;
}

/* ──────────────────────── L5: BVP Shooting (Linear) ─────────────────── */

int ode_bvp_shooting_linear(const double *p_vals, const double *q_vals,
                              const double *r_vals, const double *t_vals,
                              int n_pts, double alpha, double beta,
                              double *y_out) {
    if (!p_vals || !q_vals || !r_vals || !t_vals || !y_out || n_pts < 2)
        return -1;

    double h = t_vals[1] - t_vals[0];
    int n = n_pts;

    /* Solve two IVPs:
     * IVP1: y₁(a)=α, y₁'(a)=0 → solution y₁(x)
     * IVP2: y₂(a)=0, y₂'(a)=1 → solution y₂(x)
     * Then y(x) = y₁(x) + ((β - y₁(b))/y₂(b)) * y₂(x)
     */

    /* Allocate: y1[0..n-1], y1p[0..n-1], y2[0..n-1], y2p[0..n-1] */
    double *y1  = (double *)calloc((size_t)n, sizeof(double));
    double *y2  = (double *)calloc((size_t)n, sizeof(double));
    double *y1p = (double *)calloc((size_t)n, sizeof(double));
    double *y2p = (double *)calloc((size_t)n, sizeof(double));
    if (!y1 || !y2 || !y1p || !y2p) { free(y1); free(y2); free(y1p); free(y2p); return -1; }

    /* IVP1 initial conditions */
    y1[0] = alpha;  y1p[0] = 0.0;
    /* IVP2 initial conditions */
    y2[0] = 0.0;    y2p[0] = 1.0;

    /* Euler integration for both */
    for (int i = 0; i < n - 1; i++) {
        double dt = t_vals[i+1] - t_vals[i];
        /* y'' = p*y' + q*y + r */
        /* Convert to first-order system: u₁=y, u₂=y' */
        /* y1: */
        double y1pp = p_vals[i] * y1p[i] + q_vals[i] * y1[i] + r_vals[i];
        y1[i+1]  = y1[i]  + dt * y1p[i];
        y1p[i+1] = y1p[i] + dt * y1pp;
        /* y2: */
        double y2pp = p_vals[i] * y2p[i] + q_vals[i] * y2[i]; /* homogeneous for y2 */
        y2[i+1]  = y2[i]  + dt * y2p[i];
        y2p[i+1] = y2p[i] + dt * y2pp;
    }

    /* Adjust: the true initial slope */
    double adjust = (fabs(y2[n-1]) > 1e-14) ? (beta - y1[n-1]) / y2[n-1] : 0.0;

    /* Solution: y(x) = y1(x) + adjust * y2(x) */
    for (int i = 0; i < n; i++) {
        y_out[i] = y1[i] + adjust * y2[i];
    }

    free(y1); free(y2); free(y1p); free(y2p);
    return 0;
}

/* ──────────────────────── L5: BVP Shooting (Nonlinear) ──────────────── */

int ode_bvp_shooting_nonlinear(
    double (*f_func)(double t, double y, double yp, void *ctx),
    double (*fy_func)(double t, double y, double yp, void *ctx),
    double (*fyp_func)(double t, double y, double yp, void *ctx),
    void *ctx, const double *t_vals, int n_pts,
    double alpha, double beta, double *y_out)
{
    if (!f_func || !t_vals || !y_out || n_pts < 2) return -1;

    double h = t_vals[1] - t_vals[0];
    int n = n_pts;

    /* Newton iteration on initial slope s = y'(a) */
    double s = (beta - alpha) / (t_vals[n-1] - t_vals[0]); /* secant guess */

    /* Solve IVP with initial slope s */
    double *y  = (double *)calloc((size_t)n, sizeof(double));
    double *yp = (double *)calloc((size_t)n, sizeof(double));
    if (!y || !yp) { free(y); free(yp); return -1; }

    for (int newton = 0; newton < 20; newton++) {
        /* IVP: y'' = f(t, y, y'), y(a)=alpha, y'(a)=s */
        y[0] = alpha; yp[0] = s;
        for (int i = 0; i < n - 1; i++) {
            double dt = t_vals[i+1] - t_vals[i];
            double ypp = f_func(t_vals[i], y[i], yp[i], ctx);
            y[i+1]  = y[i]  + dt * yp[i];
            yp[i+1] = yp[i] + dt * ypp;
        }

        /* Residual: F(s) = y(b; s) - beta */
        double F = y[n-1] - beta;
        if (fabs(F) < 1e-10) break;

        /* Compute derivative dF/ds via variational equation */
        /* Solve: z''= fy*z + fyp*z', z(0)=0, z'(0)=1 */
        double *z  = (double *)calloc((size_t)n, sizeof(double));
        double *zp = (double *)calloc((size_t)n, sizeof(double));
        if (!z || !zp) { free(z); free(zp); break; }

        z[0] = 0.0; zp[0] = 1.0;
        for (int i = 0; i < n - 1; i++) {
            double dt = t_vals[i+1] - t_vals[i];
            double fy  = fy_func ? fy_func(t_vals[i], y[i], yp[i], ctx) : 0.0;
            double fyp = fyp_func ? fyp_func(t_vals[i], y[i], yp[i], ctx) : 0.0;
            double zpp = fy * z[i] + fyp * zp[i];
            z[i+1]  = z[i]  + dt * zp[i];
            zp[i+1] = zp[i] + dt * zpp;
        }

        double dF = z[n-1];
        free(z); free(zp);

        if (fabs(dF) < 1e-14) break;
        s = s - F / dF;
    }

    for (int i = 0; i < n; i++) y_out[i] = y[i];
    free(y); free(yp);
    return 0;
}

/* ──────────────────────── L5: BVP Finite Difference (Linear) ────────── */

int ode_bvp_fd_linear(const double *p_vals, const double *q_vals,
                        const double *r_vals, double a, double b,
                        double alpha, double beta, int n_int,
                        double *t_vals, double *y_vals) {
    if (!y_vals || n_int < 1) return -1;

    int n = n_int + 2; /* total points including boundaries */
    double h = (b - a) / (n_int + 1);

    /* Fill time values */
    for (int i = 0; i < n; i++) t_vals[i] = a + i * h;

    /* Boundary conditions */
    y_vals[0] = alpha;
    y_vals[n-1] = beta;

    /* Build tridiagonal system for interior points */
    /* Equation at x_i: (y_{i+1} - 2y_i + y_{i-1})/h² = p_i*(y_{i+1}-y_{i-1})/(2h) + q_i*y_i + r_i */
    /* Rearranged: */
    /* a_i = -1/h² - p_i/(2h) */
    /* b_i = 2/h² + q_i */
    /* c_i = -1/h² + p_i/(2h) */
    /* RHS_i = -r_i (+ boundary contributions) */

    double *a_diag = (double *)malloc((size_t)n_int * sizeof(double));
    double *b_diag = (double *)malloc((size_t)n_int * sizeof(double));
    double *c_diag = (double *)malloc((size_t)n_int * sizeof(double));
    double *rhs    = (double *)malloc((size_t)n_int * sizeof(double));
    if (!a_diag || !b_diag || !c_diag || !rhs) {
        free(a_diag); free(b_diag); free(c_diag); free(rhs); return -1;
    }

    double h2 = h * h;
    for (int i = 0; i < n_int; i++) {
        double p = p_vals ? p_vals[i] : 0.0;
        double q = q_vals ? q_vals[i] : 0.0;
        double r = r_vals ? r_vals[i] : 0.0;

        a_diag[i] = -1.0/h2 - p/(2.0*h);
        b_diag[i] =  2.0/h2 + q;
        c_diag[i] = -1.0/h2 + p/(2.0*h);
        rhs[i]    = -r;
    }

    /* Adjust for boundary conditions */
    rhs[0] -= a_diag[0] * alpha;
    rhs[n_int-1] -= c_diag[n_int-1] * beta;

    /* Thomas algorithm for tridiagonal solve */
    /* Forward sweep */
    for (int i = 1; i < n_int; i++) {
        double m = a_diag[i] / b_diag[i-1];
        b_diag[i] -= m * c_diag[i-1];
        rhs[i]    -= m * rhs[i-1];
    }
    /* Back substitution */
    y_vals[n_int] = rhs[n_int-1] / b_diag[n_int-1];
    for (int i = n_int - 2; i >= 0; i--) {
        y_vals[i+1] = (rhs[i] - c_diag[i] * y_vals[i+2]) / b_diag[i];
    }

    free(a_diag); free(b_diag); free(c_diag); free(rhs);
    return 0;
}

/* ──────────────────────── Utility ──────────────────────────────────── */

int ode_richardson_error(const ODEIVP *ivp, ODEMethod method,
                           double t, const double *y, double h,
                           int order, double *err_est) {
    if (!ivp || !y || !err_est || h <= 0.0 || order < 1) return -1;

    (void)method; /* Use RK4 as reference */
    int dim = ivp->dim;
    double *y_h = (double *)malloc((size_t)dim * sizeof(double));
    double *y_h2 = (double *)malloc((size_t)dim * sizeof(double));
    double *y_h2_step = (double *)malloc((size_t)dim * sizeof(double));
    if (!y_h || !y_h2 || !y_h2_step) { free(y_h); free(y_h2); free(y_h2_step); return -1; }

    /* One full step */
    ode_rk4_step(ivp, t, y, h, y_h);

    /* Two half steps */
    ode_rk4_step(ivp, t, y, h/2.0, y_h2_step);
    ode_rk4_step(ivp, t + h/2.0, y_h2_step, h/2.0, y_h2);

    /* Richardson: err = ||y_h2 - y_h|| / (2^p - 1) */
    double max_diff = 0.0;
    for (int i = 0; i < dim; i++) {
        double d = fabs(y_h2[i] - y_h[i]);
        if (d > max_diff) max_diff = d;
    }
    *err_est = max_diff / ((1 << order) - 1);

    free(y_h); free(y_h2); free(y_h2_step);
    return 0;
}
