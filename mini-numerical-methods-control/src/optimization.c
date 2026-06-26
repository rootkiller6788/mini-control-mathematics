/**
 * optimization.c — Numerical Optimization for Control Engineering
 *
 * Implements: L5 Methods (gradient descent, conjugate gradient, root finding),
 *             L8 Advanced (Newton, BFGS, line search)
 *
 * Optimization is the unifying mathematical framework for control:
 * LQR minimizes quadratic cost, MPC solves constrained optimization,
 * system identification minimizes prediction error, and robust control
 * solves min-max saddle-point problems.
 *
 * L8 Advanced: Stochastic gradient descent (SGD), Monte Carlo methods
 * for Bayesian inference in adaptive control, Lyapunov-based
 * optimization for time-varying systems, fuzzy logic tuning.
 */

#include "optimization.h"
#include "linear_systems.h"
#include <string.h>
#include <stdio.h>
#include <float.h>

/* ==========================================================================
 * L5: Root-Finding Methods
 * ========================================================================== */

double bisection_method(ScalarFunction f, void *params,
                        double a, double b, double tol, size_t max_iter,
                        int *converged) {
    *converged = 0;
    double fa = f(a, params), fb = f(b, params);
    if (fa * fb > 0.0) return a;  /* No sign change */
    if (fabs(fa) < tol) { *converged = 1; return a; }
    if (fabs(fb) < tol) { *converged = 1; return b; }

    double c;
    for (size_t i = 0; i < max_iter; i++) {
        c = (a + b) / 2.0;
        double fc = f(c, params);
        if (fabs(fc) < tol || (b - a) / 2.0 < tol) {
            *converged = 1;
            return c;
        }
        if (fa * fc < 0.0) { b = c; fb = fc; }
        else { a = c; fa = fc; }
    }
    return (a + b) / 2.0;
}

double newton_raphson(ScalarFunction f, ScalarFunction fprime, void *params,
                      double x0, double tol, size_t max_iter, int *converged) {
    *converged = 0;
    double x = x0;
    for (size_t i = 0; i < max_iter; i++) {
        double fx = f(x, params), fpx = fprime(x, params);
        if (fabs(fx) < tol) { *converged = 1; return x; }
        if (fabs(fpx) < MNC_TINY) return x;  /* Derivative too small */
        double dx = fx / fpx;
        x -= dx;
        if (fabs(dx) < tol) { *converged = 1; return x; }
    }
    return x;
}

double secant_method(ScalarFunction f, void *params,
                     double x0, double x1, double tol, size_t max_iter,
                     int *converged) {
    *converged = 0;
    double f0 = f(x0, params), f1 = f(x1, params);
    if (fabs(f1 - f0) < MNC_TINY) return x1;

    for (size_t i = 0; i < max_iter; i++) {
        if (fabs(f1) < tol) { *converged = 1; return x1; }
        double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        x0 = x1; f0 = f1;
        x1 = x2; f1 = f(x1, params);
        if (fabs(x1 - x0) < tol) { *converged = 1; return x1; }
    }
    return x1;
}

/* ==========================================================================
 * L5: Gradient Descent (First-Order Methods)
 * ========================================================================== */

OptResult* gradient_descent(MultivariateFunction f, MultivariateGradient grad,
                            void *params, size_t n, const double *x0,
                            double step_size, double tol, size_t max_iter) {
    OptResult *r = (OptResult*)malloc(sizeof(OptResult));
    if (!r) return NULL;
    r->x_opt = (double*)malloc(n * sizeof(double));
    if (!r->x_opt) { free(r); return NULL; }
    memcpy(r->x_opt, x0, n * sizeof(double));

    double *g = (double*)malloc(n * sizeof(double));
    if (!g) { free(r->x_opt); free(r); return NULL; }

    r->converged = 0;
    for (r->iterations = 0; r->iterations < max_iter; r->iterations++) {
        grad(r->x_opt, n, params, g);

        double gnorm = 0.0;
        for (size_t i = 0; i < n; i++) gnorm += g[i] * g[i];
        gnorm = sqrt(gnorm);
        r->gradient_norm = gnorm;
        if (gnorm < tol) { r->converged = 1; break; }

        for (size_t i = 0; i < n; i++)
            r->x_opt[i] -= step_size * g[i];
    }
    r->f_opt = f(r->x_opt, n, params);
    free(g);
    return r;
}

OptResult* gradient_descent_backtracking(MultivariateFunction f,
                                         MultivariateGradient grad,
                                         void *params, size_t n,
                                         const double *x0, double init_step,
                                         double tol, size_t max_iter) {
    OptResult *r = (OptResult*)malloc(sizeof(OptResult));
    if (!r) return NULL;
    r->x_opt = (double*)malloc(n * sizeof(double));
    if (!r->x_opt) { free(r); return NULL; }
    memcpy(r->x_opt, x0, n * sizeof(double));

    double *g = (double*)malloc(n * sizeof(double));
    double *x_try = (double*)malloc(n * sizeof(double));
    if (!g || !x_try) { free(g); free(x_try); free(r->x_opt); free(r); return NULL; }

    double c = 1e-4;  /* Armijo constant */
    r->converged = 0;

    for (r->iterations = 0; r->iterations < max_iter; r->iterations++) {
        double fk = f(r->x_opt, n, params);
        grad(r->x_opt, n, params, g);

        double gnorm = 0.0, gdotd = 0.0;
        for (size_t i = 0; i < n; i++) {
            gnorm += g[i] * g[i];
            gdotd -= g[i] * g[i];  /* d = -g(x) */
        }
        gnorm = sqrt(gnorm);
        r->gradient_norm = gnorm;
        if (gnorm < tol) { r->converged = 1; break; }

        /* Armijo backtracking line search */
        double alpha = init_step;
        while (alpha > MNC_TINY) {
            for (size_t i = 0; i < n; i++)
                x_try[i] = r->x_opt[i] - alpha * g[i];
            double f_try = f(x_try, n, params);
            if (f_try <= fk + c * alpha * gdotd) break;
            alpha *= 0.5;
        }

        for (size_t i = 0; i < n; i++)
            r->x_opt[i] -= alpha * g[i];
    }
    r->f_opt = f(r->x_opt, n, params);
    free(g); free(x_try);
    return r;
}

OptResult* gradient_descent_momentum(MultivariateFunction f,
                                     MultivariateGradient grad,
                                     void *params, size_t n, const double *x0,
                                     double step_size, double momentum,
                                     double tol, size_t max_iter) {
    OptResult *r = (OptResult*)malloc(sizeof(OptResult));
    if (!r) return NULL;
    r->x_opt = (double*)malloc(n * sizeof(double));
    if (!r->x_opt) { free(r); return NULL; }
    memcpy(r->x_opt, x0, n * sizeof(double));

    double *g = (double*)malloc(n * sizeof(double));
    double *v = (double*)calloc(n, sizeof(double));
    if (!g || !v) { free(g); free(v); free(r->x_opt); free(r); return NULL; }

    r->converged = 0;
    for (r->iterations = 0; r->iterations < max_iter; r->iterations++) {
        grad(r->x_opt, n, params, g);

        double gnorm = 0.0;
        for (size_t i = 0; i < n; i++) {
            v[i] = momentum * v[i] + step_size * g[i];
            r->x_opt[i] -= v[i];
            gnorm += g[i] * g[i];
        }
        r->gradient_norm = sqrt(gnorm);
        if (r->gradient_norm < tol) { r->converged = 1; break; }
    }
    r->f_opt = f(r->x_opt, n, params);
    free(g); free(v);
    return r;
}

/* ==========================================================================
 * L5: Conjugate Gradient (Linear)
 * ========================================================================== */

Vector* conjugate_gradient_linear(const Matrix *A, const Vector *b,
                                  const Vector *x0, double tol, size_t max_iter) {
    if (!A || !b || !x0 || A->rows != A->cols || A->rows != b->n) return NULL;
    size_t n = A->rows;

    Vector *x = vector_copy(x0);
    if (!x) return NULL;
    Vector *r = vector_copy(b);
    Vector *p = vector_copy(b);
    if (!r || !p) { vector_free(x); vector_free(r); vector_free(p); return NULL; }

    /* Initial residual r = b - Ax */
    for (size_t i = 0; i < n; i++) {
        double ax = 0.0;
        for (size_t j = 0; j < n; j++)
            ax += A->data[i * A->cols + j] * x->data[j];
        r->data[i] = b->data[i] - ax;
        p->data[i] = r->data[i];
    }

    double rsold = vector_dot(r, r);
    if (sqrt(rsold) < tol) { vector_free(r); vector_free(p); return x; }

    for (size_t k = 0; k < max_iter && k < n; k++) {
        /* Ap = A * p */
        double *Ap = (double*)malloc(n * sizeof(double));
        if (!Ap) break;
        for (size_t i = 0; i < n; i++) {
            Ap[i] = 0.0;
            for (size_t j = 0; j < n; j++)
                Ap[i] += A->data[i * A->cols + j] * p->data[j];
        }

        double pAp = 0.0;
        for (size_t i = 0; i < n; i++) pAp += p->data[i] * Ap[i];
        if (fabs(pAp) < MNC_TINY) { free(Ap); break; }

        double alpha = rsold / pAp;
        for (size_t i = 0; i < n; i++) {
            x->data[i] += alpha * p->data[i];
            r->data[i] -= alpha * Ap[i];
        }

        double rsnew = vector_dot(r, r);
        if (sqrt(rsnew) < tol) { free(Ap); break; }

        double beta = rsnew / rsold;
        for (size_t i = 0; i < n; i++)
            p->data[i] = r->data[i] + beta * p->data[i];
        rsold = rsnew;
        free(Ap);
    }

    vector_free(r); vector_free(p);
    return x;
}

/* ==========================================================================
 * L8: Nonlinear Conjugate Gradient (Fletcher-Reeves)
 * ========================================================================== */

OptResult* conjugate_gradient_nonlinear(MultivariateFunction f,
                                        MultivariateGradient grad,
                                        void *params, size_t n,
                                        const double *x0, double tol,
                                        size_t max_iter) {
    OptResult *r = (OptResult*)malloc(sizeof(OptResult));
    if (!r) return NULL;
    r->x_opt = (double*)malloc(n * sizeof(double));
    if (!r->x_opt) { free(r); return NULL; }
    memcpy(r->x_opt, x0, n * sizeof(double));

    double *g = (double*)malloc(n * sizeof(double));
    double *g_old = (double*)malloc(n * sizeof(double));
    double *d = (double*)malloc(n * sizeof(double));
    if (!g || !g_old || !d) {
        free(g); free(g_old); free(d); free(r->x_opt); free(r); return NULL;
    }

    grad(r->x_opt, n, params, g);
    double gnorm = 0.0;
    for (size_t i = 0; i < n; i++) { gnorm += g[i] * g[i]; d[i] = -g[i]; }
    r->gradient_norm = sqrt(gnorm);
    r->converged = 0;

    for (r->iterations = 0; r->iterations < max_iter; r->iterations++) {
        if (r->gradient_norm < tol) { r->converged = 1; break; }

        /* Wolfe line search along d */
        double alpha = line_search_wolfe(f, grad, params, n, r->x_opt, d, 1.0);

        for (size_t i = 0; i < n; i++) {
            r->x_opt[i] += alpha * d[i];
            g_old[i] = g[i];
        }

        grad(r->x_opt, n, params, g);
        gnorm = 0.0;
        double gnorm_old = 0.0;
        for (size_t i = 0; i < n; i++) {
            gnorm += g[i] * g[i];
            gnorm_old += g_old[i] * g_old[i];
        }
        r->gradient_norm = sqrt(gnorm);

        /* Fletcher-Reeves beta */
        double beta_fr = (gnorm_old > MNC_TINY) ? gnorm / gnorm_old : 0.0;
        for (size_t i = 0; i < n; i++)
            d[i] = -g[i] + beta_fr * d[i];
    }
    r->f_opt = f(r->x_opt, n, params);
    free(g); free(g_old); free(d);
    return r;
}

/* ==========================================================================
 * L8: Newton's Method for Optimization
 * ========================================================================== */

OptResult* newton_method_opt(MultivariateFunction f, MultivariateGradient grad,
                             void *params, size_t n, const double *x0,
                             double tol, size_t max_iter) {
    OptResult *r = (OptResult*)malloc(sizeof(OptResult));
    if (!r) return NULL;
    r->x_opt = (double*)malloc(n * sizeof(double));
    if (!r->x_opt) { free(r); return NULL; }
    memcpy(r->x_opt, x0, n * sizeof(double));

    double *g = (double*)malloc(n * sizeof(double));
    double *s = (double*)malloc(n * sizeof(double));
    double *h = (double*)malloc(n * n * sizeof(double));
    if (!g || !s || !h) {
        free(g); free(s); free(h); free(r->x_opt); free(r); return NULL;
    }

    r->converged = 0;
    for (r->iterations = 0; r->iterations < max_iter; r->iterations++) {
        grad(r->x_opt, n, params, g);

        double gnorm = 0.0;
        for (size_t i = 0; i < n; i++) gnorm += g[i] * g[i];
        r->gradient_norm = sqrt(gnorm);
        if (gnorm < tol) { r->converged = 1; break; }

        /* Finite-difference Hessian: H_{ij} ≈ (∂g_i/∂x_j) */
        double eps = 1e-6;
        for (size_t j = 0; j < n; j++) {
            double xj_save = r->x_opt[j];
            r->x_opt[j] += eps;
            grad(r->x_opt, n, params, s);
            r->x_opt[j] = xj_save;
            for (size_t i = 0; i < n; i++)
                h[i * n + j] = (s[i] - g[i]) / eps;
        }

        /* Solve H * dx = -g */
        Matrix *Hmat = matrix_create_from(n, n, h);
        Vector *gvec = vector_create_from(n, g);
        Vector *dx = vector_create(n);
        if (Hmat && gvec && dx && gaussian_elimination_solve(Hmat, gvec, dx) == 0) {
            for (size_t i = 0; i < n; i++)
                r->x_opt[i] -= dx->data[i];
        }
        matrix_free(Hmat); vector_free(gvec); vector_free(dx);
    }
    r->f_opt = f(r->x_opt, n, params);
    free(g); free(s); free(h);
    return r;
}

/* ==========================================================================
 * L8: BFGS Update and Wolfe Line Search
 * ========================================================================== */

void bfgs_update(Matrix *B, const double *s, const double *y, size_t n) {
    if (!B || !s || !y) return;
    double yTs = 0.0, sTBs = 0.0;
    double *Bs = (double*)malloc(n * sizeof(double));
    if (!Bs) return;

    for (size_t i = 0; i < n; i++) {
        Bs[i] = 0.0;
        for (size_t j = 0; j < n; j++)
            Bs[i] += B->data[i * B->cols + j] * s[j];
    }

    for (size_t i = 0; i < n; i++) {
        yTs += y[i] * s[i];
        sTBs += s[i] * Bs[i];
    }

    if (yTs <= MNC_TINY) { free(Bs); return; }

    /* B += y*y^T / y^T s - (B s)(B s)^T / s^T B s */
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            B->data[i * B->cols + j] += (y[i] * y[j]) / yTs
                                       - (Bs[i] * Bs[j]) / sTBs;
    free(Bs);
}

double line_search_wolfe(MultivariateFunction f, MultivariateGradient grad,
                         void *params, size_t n, const double *x,
                         const double *direction, double max_step) {
    double c1 = 1e-4, c2 = 0.9;
    double alpha = max_step;
    double f0 = f(x, n, params);

    double *g0 = (double*)malloc(n * sizeof(double));
    double *xt = (double*)malloc(n * sizeof(double));
    double *gt = (double*)malloc(n * sizeof(double));
    if (!g0 || !xt || !gt) { free(g0); free(xt); free(gt); return 1e-6; }

    grad(x, n, params, g0);
    double dg0 = 0.0;
    for (size_t i = 0; i < n; i++) dg0 += g0[i] * direction[i];

    for (size_t i = 0; i < 20; i++) {
        for (size_t j = 0; j < n; j++)
            xt[j] = x[j] + alpha * direction[j];
        double ft = f(xt, n, params);

        /* Armijo condition */
        if (ft > f0 + c1 * alpha * dg0) { alpha *= 0.5; continue; }

        /* Curvature condition */
        grad(xt, n, params, gt);
        double dgt = 0.0;
        for (size_t j = 0; j < n; j++) dgt += gt[j] * direction[j];
        if (dgt < c2 * dg0) { alpha *= 2.0; continue; }

        break;
    }

    free(g0); free(xt); free(gt);
    return alpha;
}

void opt_result_free(OptResult *r) {
    if (!r) return;
    free(r->x_opt);
    free(r);
}
