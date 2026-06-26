/**
 * newton_methods.c — Newton-Type Second-Order Methods Implementation
 *
 * Knowledge points:
 * L5: Newton's method with backtracking line search
 * L3: Newton decrement (affine-invariant optimality measure)
 * L5: Damped Newton for self-concordant functions
 * L5: BFGS quasi-Newton method (most popular in practice)
 * L5: DFP quasi-Newton method (historical first)
 * L5: SR1 symmetric rank-one update
 * L8: L-BFGS limited-memory BFGS for large-scale problems
 * L8: Trust-region Newton (Levenberg-Marquardt)
 * L5: Gauss-Newton for nonlinear least squares
 * L5: Levenberg-Marquardt for robust NLS
 *
 * Reference: Nocedal & Wright (2006), Dennis & Schnabel (1996)
 */

#include "newton_methods.h"
#include "gradient_methods.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>

/* ═══════════════════════════════════════════════════════════════════════
 * Configuration Defaults
 * ═══════════════════════════════════════════════════════════════════════ */

newton_config_t newton_config_default(void)
{
    newton_config_t cfg;
    cfg.max_iters = 200;
    cfg.grad_tol = GRAD_TOL;
    cfg.alpha_min = 1e-16;
    cfg.beta = 0.5;
    cfg.c1 = 1e-4;
    cfg.use_cholesky = 1;
    cfg.verbose = 0;
    return cfg;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Newton Decrement
 *
 * lambda(x) = (grad f(x)^T [Hess f(x)]^{-1} grad f(x))^{1/2}
 *
 * This is the norm of the Newton step in the Hessian metric.
 * It provides an affine-invariant measure of proximity to the optimum.
 *
 * Properties (for self-concordant functions):
 * - lambda(x) = 0 iff x is a stationary point
 * - lambda(x) < 1 => f(x) - f* <= lambda(x)^2
 * - lambda(x) < 3-sqrt(5) ≈ 0.26 => pure Newton converges quadratically
 *
 * The Newton decrement squared equals the predicted reduction:
 *     lambda(x)^2 = -grad f(x)^T * Delta_x_nt
 * ═══════════════════════════════════════════════════════════════════════ */

double newton_decrement(const convex_function_t *f, const vector_t *x,
                        const vector_t *grad, matrix_t *H, vector_t *work)
{
    assert(f && x && grad && H && work);
    size_t n = x->n;

    /* Compute Hessian at x */
    f->hess(f, x, H);

    /* Solve H * d_nt = -grad for Newton step d_nt in work */
    /* Copy to work vector: work = -grad */
    for (size_t i = 0; i < n; i++) {
        work->data[i] = -grad->data[i];
    }

    /* Try Cholesky solve */
    matrix_t *H_copy = matrix_alloc(n, n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            MAT_ELT(H_copy, i, j) = MAT_ELT(H, i, j);
        }
    }

    if (cholesky_decomp(H_copy) == 0) {
        /* H is PD; solve H * d = -grad */
        cholesky_solve(H_copy, work, work);
    } else {
        /* H is not PD; use gradient step direction scaled by ||grad|| */
        double gnorm = vector_norm_l2(grad);
        if (gnorm > 1e-12) {
            for (size_t i = 0; i < n; i++) {
                work->data[i] = -grad->data[i] / gnorm;
            }
        }
    }

    matrix_free(H_copy);

    /* lambda^2 = -grad^T * d_nt (= grad^T H^{-1} grad) */
    double lambda2 = 0.0;
    for (size_t i = 0; i < n; i++) {
        lambda2 += grad->data[i] * (-work->data[i]);
    }

    return (lambda2 > 0.0) ? sqrt(lambda2) : 0.0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Newton's Method with Backtracking Line Search
 *
 * Pure Newton step: d_k = -[Hess f(x_k)]^{-1} grad f(x_k)
 *
 * This is derived from the second-order Taylor expansion:
 *     f(x + d) ≈ f(x) + grad^T d + (1/2) d^T Hess d
 *
 * The Newton step minimizes this quadratic model exactly.
 *
 * Theorem (local quadratic convergence):
 * If grad f(x*) = 0, Hess f(x*) > 0, and Hess f is Lipschitz
 * continuous near x*, then for x_0 sufficiently close to x*:
 *     ||x_{k+1} - x*|| <= C ||x_k - x*||^2
 *
 * This is the hallmark of Newton's method: doubling the number
 * of correct digits at each iteration in the local convergence phase.
 *
 * Globalization: Use backtracking line search to ensure decrease
 * when far from the solution. The full Newton step (alpha=1) is
 * eventually accepted near the optimum.
 * ═══════════════════════════════════════════════════════════════════════ */

void newton_method(const unconstr_problem_t *problem,
                   const newton_config_t *config,
                   opt_result_t *result)
{
    assert(problem && config && result);
    assert(problem->f.eval && problem->f.grad && problem->f.hess);

    size_t n = problem->n;
    vector_t *x = vector_alloc(n);
    vector_t *g = vector_alloc(n);
    vector_t *d = vector_alloc(n);
    matrix_t *H = matrix_alloc(n, n);
    matrix_t *L = matrix_alloc(n, n);

    vector_copy(&problem->x0, x);

    backtracking_params_t ls_params;
    ls_params.alpha_init = 1.0;
    ls_params.rho = config->beta;
    ls_params.c1 = config->c1;
    ls_params.c2 = 0.9;
    ls_params.max_iter = 50;

    size_t iter;
    for (iter = 0; iter < config->max_iters; iter++) {
        problem->f.grad(&problem->f, x, g);
        double gnorm = vector_norm_l2(g);

        if (gnorm < config->grad_tol) break;

        /* Compute Hessian */
        problem->f.hess(&problem->f, x, H);

        /* Compute Newton step direction: d = -H^{-1} g */
        if (config->use_cholesky) {
            /* Copy H to L for Cholesky */
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    MAT_ELT(L, i, j) = MAT_ELT(H, i, j);
                }
            }

            if (cholesky_decomp(L) == 0) {
                /* H is PD; solve H * d = -g */
                for (size_t i = 0; i < n; i++) {
                    d->data[i] = -g->data[i];
                }
                cholesky_solve(L, d, d);
            } else {
                /* H not PD; fallback to negative gradient */
                for (size_t i = 0; i < n; i++) {
                    d->data[i] = -g->data[i];
                }
            }
        } else {
            /* No Cholesky; use gradient descent */
            for (size_t i = 0; i < n; i++) {
                d->data[i] = -g->data[i];
            }
        }

        /* Check Newton decrement for early termination */
        double lambda2 = -vector_dot(g, d);
        if (lambda2 / 2.0 < config->grad_tol) break;

        /* Line search */
        double alpha = backtracking_linesearch(&problem->f, x, d, g, &ls_params);
        if (alpha < config->alpha_min) {
            result->status = OPT_LINE_SEARCH_FAILED;
            break;
        }

        /* Update: x += alpha * d */
        for (size_t i = 0; i < n; i++) {
            x->data[i] += alpha * d->data[i];
        }
    }

    result->iterations = iter;
    result->grad_norm = vector_norm_l2(g);
    result->f_value = problem->f.eval(&problem->f, x);

    if (result->status != OPT_LINE_SEARCH_FAILED) {
        result->status = (iter < config->max_iters) ? OPT_SUCCESS
                                                     : OPT_MAXITER_REACHED;
    }

    if (!result->x.data) {
        result->x.data = (double*)malloc(n * sizeof(double));
        result->x.n = n;
    }
    memcpy(result->x.data, x->data, n * sizeof(double));

    vector_free(x); vector_free(g); vector_free(d);
    matrix_free(H); matrix_free(L);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Damped Newton's Method for Self-Concordant Functions
 *
 * For self-concordant f, the step size is chosen based on the
 * Newton decrement lambda(x):
 *     alpha = 1 if lambda(x) <= 1/4
 *     alpha = 1/(1+lambda(x)) otherwise
 *
 * Theorem (Nesterov & Nemirovski, 1994):
 * Damped Newton for self-concordant f converges globally with:
 * - Damped phase: f decreases by at least constant > 0
 * - Pure phase: quadratic convergence once lambda < 1/4
 * Total iterations: O(log log(1/epsilon))
 *
 * Self-concordance condition:
 * |f'''(x)| <= 2 * f''(x)^(3/2)
 *
 * Examples of self-concordant functions: log-barrier, log-det,
 * entropy, log-sum-exp.
 *
 * Reference: Nesterov & Nemirovski, "Interior-Point Polynomial
 *            Algorithms in Convex Programming" (1994)
 * ═══════════════════════════════════════════════════════════════════════ */

void damped_newton_method(const unconstr_problem_t *problem,
                          const newton_config_t *config,
                          opt_result_t *result)
{
    assert(problem && config && result && problem->f.hess);

    size_t n = problem->n;
    vector_t *x = vector_alloc(n);
    vector_t *g = vector_alloc(n);
    vector_t *d = vector_alloc(n);
    matrix_t *H = matrix_alloc(n, n);
    matrix_t *L = matrix_alloc(n, n);

    vector_copy(&problem->x0, x);

    size_t iter;
    for (iter = 0; iter < config->max_iters; iter++) {
        problem->f.grad(&problem->f, x, g);
        problem->f.hess(&problem->f, x, H);

        /* Compute Newton decrement */
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                MAT_ELT(L, i, j) = MAT_ELT(H, i, j);
            }
        }

        double lam;
        if (cholesky_decomp(L) == 0) {
            for (size_t i = 0; i < n; i++) d->data[i] = -g->data[i];
            cholesky_solve(L, d, d);
            /* Newton decrement = sqrt(g^T H^{-1} g) = sqrt(-g^T d) */
            lam = 0.0;
            for (size_t i = 0; i < n; i++) lam += g->data[i] * (-d->data[i]);
            lam = (lam > 0.0) ? sqrt(lam) : 0.0;
        } else {
            lam = vector_norm_l2(g);
            for (size_t i = 0; i < n; i++) d->data[i] = -g->data[i];
        }

        if (lam * lam / 2.0 < config->grad_tol) break;

        /* Damped step selection */
        double alpha = (lam <= 0.25) ? 1.0 : (1.0 / (1.0 + lam));

        /* Line search backup */
        double fx = problem->f.eval(&problem->f, x);
        double dir_deriv = vector_dot(g, d);
        while (alpha > config->alpha_min) {
            for (size_t i = 0; i < n; i++) {
                x->data[i] += alpha * d->data[i];
            }
            double f_new = problem->f.eval(&problem->f, x);
            if (f_new <= fx + config->c1 * alpha * dir_deriv) {
                break;
            }
            /* Undo step */
            for (size_t i = 0; i < n; i++) {
                x->data[i] -= alpha * d->data[i];
            }
            alpha *= config->beta;
        }
    }

    result->iterations = iter;
    result->f_value = problem->f.eval(&problem->f, x);
    result->grad_norm = vector_norm_l2(g);
    result->status = (iter < config->max_iters) ? OPT_SUCCESS
                                                 : OPT_MAXITER_REACHED;

    if (!result->x.data) {
        result->x.data = (double*)malloc(n * sizeof(double));
        result->x.n = n;
    }
    memcpy(result->x.data, x->data, n * sizeof(double));

    vector_free(x); vector_free(g); vector_free(d);
    matrix_free(H); matrix_free(L);
}

/* ═══════════════════════════════════════════════════════════════════════
 * BFGS Method (Broyden-Fletcher-Goldfarb-Shanno)
 *
 * Maintains an approximation H_k ≈ [Hess f(x_k)]^{-1} updated via:
 *     H_{k+1} = (I - rho_k s_k y_k^T) H_k (I - rho_k y_k s_k^T)
 *              + rho_k s_k s_k^T
 *
 * where s_k = x_{k+1} - x_k, y_k = grad f(x_{k+1}) - grad f(x_k),
 * rho_k = 1/(y_k^T s_k).
 *
 * BFGS is the most widely used quasi-Newton method because:
 * 1. Superlinear convergence under mild conditions
 * 2. Maintains positive definiteness of H_k (=> descent directions)
 * 3. Self-correcting: errors in H_k are damped over iterations
 * 4. Works well with inexact line searches (Wolfe conditions)
 *
 * Theorem (Powell, 1976; Byrd, Nocedal, Yuan, 1987):
 * For convex f with bounded level sets, BFGS with Wolfe line search
 * converges globally with liminf ||grad f(x_k)|| = 0.
 * Superlinear convergence rate under additional smoothness.
 *
 * Reference: Broyden (1970), Fletcher (1970), Goldfarb (1970), Shanno (1970)
 * ═══════════════════════════════════════════════════════════════════════ */

void bfgs_method(const unconstr_problem_t *problem,
                 size_t max_iters, opt_result_t *result)
{
    assert(problem && result);
    size_t n = problem->n;

    vector_t *x = vector_alloc(n);
    vector_t *x_new = vector_alloc(n);
    vector_t *g = vector_alloc(n);
    vector_t *g_new = vector_alloc(n);
    vector_t *s = vector_alloc(n);   /* Step: x_{k+1} - x_k */
    vector_t *y = vector_alloc(n);   /* Gradient difference */
    vector_t *Hs = vector_alloc(n);  /* H_k * s */
    /* H_k stored as dense matrix (initialized to I) */
    matrix_t *H = matrix_alloc(n, n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            MAT_ELT(H, i, j) = (i == j) ? 1.0 : 0.0;
        }
    }

    vector_copy(&problem->x0, x);
    problem->f.grad(&problem->f, x, g);

    backtracking_params_t ls_params;
    ls_params.alpha_init = 1.0;
    ls_params.rho = 0.5;
    ls_params.c1 = 1e-4;
    ls_params.c2 = 0.9;
    ls_params.max_iter = 30;

    size_t iter;
    for (iter = 0; iter < max_iters; iter++) {
        if (vector_norm_l2(g) < GRAD_TOL) break;

        /* Search direction: d = -H_k * g_k */
        vector_t *d = vector_alloc(n);
        matrix_vector_mul(H, g, d);
        for (size_t i = 0; i < n; i++) d->data[i] = -d->data[i];

        /* Line search */
        double alpha = backtracking_linesearch(&problem->f, x, d, g, &ls_params);

        /* Update: x_new = x + alpha * d */
        for (size_t i = 0; i < n; i++) {
            x_new->data[i] = x->data[i] + alpha * d->data[i];
        }

        /* Compute new gradient */
        problem->f.grad(&problem->f, x_new, g_new);

        /* Compute s and y */
        for (size_t i = 0; i < n; i++) {
            s->data[i] = x_new->data[i] - x->data[i];
            y->data[i] = g_new->data[i] - g->data[i];
        }

        double ys = vector_dot(y, s);

        if (ys > 1e-16) {
            double rho = 1.0 / ys;

            /* Compute Hs = H * y */
            matrix_vector_mul(H, y, Hs);

            /* BFGS update: H_new = H - rho*(s*Hs^T + Hs*s^T)
             *                       + rho*(1 + rho*y^T*Hs)*(s*s^T) */
            double yHs = vector_dot(y, Hs);

            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    double term1 = rho * (s->data[i] * Hs->data[j]
                                         + Hs->data[i] * s->data[j]);
                    double term2 = rho * (1.0 + rho * yHs)
                                  * s->data[i] * s->data[j];
                    MAT_ELT(H, i, j) += term2 - term1;
                }
            }
        }
        /* If ys <= 0, skip update to preserve positive definiteness */

        /* Shift for next iteration */
        vector_copy(x_new, x);
        vector_copy(g_new, g);
        vector_free(d);
    }

    result->iterations = iter;
    result->f_value = problem->f.eval(&problem->f, x);
    result->grad_norm = vector_norm_l2(g);
    result->status = (iter < max_iters) ? OPT_SUCCESS : OPT_MAXITER_REACHED;

    if (!result->x.data) {
        result->x.data = (double*)malloc(n * sizeof(double));
        result->x.n = n;
    }
    memcpy(result->x.data, x->data, n * sizeof(double));

    vector_free(x); vector_free(x_new);
    vector_free(g); vector_free(g_new);
    vector_free(s); vector_free(y); vector_free(Hs);
    matrix_free(H);
}

/* ═══════════════════════════════════════════════════════════════════════
 * DFP Method (Davidon-Fletcher-Powell)
 *
 * The dual of BFGS. The update for H_k (inverse Hessian approximation):
 *     H_{k+1} = H_k + (s_k s_k^T)/(y_k^T s_k)
 *                    - (H_k y_k y_k^T H_k)/(y_k^T H_k y_k)
 *
 * Historical significance: This was the first quasi-Newton method
 * (Davidon, 1959), originally developed for atomic energy research.
 * While BFGS generally performs better in practice, DFP is
 * theoretically interesting as it was derived from the conjugate
 * gradient perspective.
 *
 * DFP maintains positive definiteness under the same conditions as BFGS.
 *
 * Reference: Davidon (1959, published 1991), "Variable metric method
 *            for minimization"; Fletcher & Powell (1963)
 * ═══════════════════════════════════════════════════════════════════════ */

void dfp_method(const unconstr_problem_t *problem,
                size_t max_iters, opt_result_t *result)
{
    assert(problem && result);
    size_t n = problem->n;

    vector_t *x = vector_alloc(n), *x_new = vector_alloc(n);
    vector_t *g = vector_alloc(n), *g_new = vector_alloc(n);
    vector_t *s = vector_alloc(n), *y = vector_alloc(n);
    vector_t *Hy = vector_alloc(n);
    matrix_t *H = matrix_alloc(n, n);

    /* Initialize H = I */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            MAT_ELT(H, i, j) = (i == j) ? 1.0 : 0.0;
        }
    }

    vector_copy(&problem->x0, x);
    problem->f.grad(&problem->f, x, g);

    backtracking_params_t ls_params;
    ls_params.alpha_init = 1.0; ls_params.rho = 0.5;
    ls_params.c1 = 1e-4; ls_params.c2 = 0.9; ls_params.max_iter = 30;

    size_t iter;
    for (iter = 0; iter < max_iters; iter++) {
        if (vector_norm_l2(g) < GRAD_TOL) break;

        vector_t *d = vector_alloc(n);
        matrix_vector_mul(H, g, d);
        for (size_t i = 0; i < n; i++) d->data[i] = -d->data[i];

        double alpha = backtracking_linesearch(&problem->f, x, d, g, &ls_params);

        for (size_t i = 0; i < n; i++) {
            x_new->data[i] = x->data[i] + alpha * d->data[i];
        }

        problem->f.grad(&problem->f, x_new, g_new);

        /* s = step, y = grad difference */
        for (size_t i = 0; i < n; i++) {
            s->data[i] = x_new->data[i] - x->data[i];
            y->data[i] = g_new->data[i] - g->data[i];
        }

        double ys = vector_dot(y, s);
        matrix_vector_mul(H, y, Hy);
        double yHy = vector_dot(y, Hy);

        if (ys > 1e-16 && yHy > 1e-16) {
            /* DFP update */
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    /* H += s*s^T/ys - (Hy)*(Hy)^T/yHy */
                    MAT_ELT(H, i, j) += s->data[i] * s->data[j] / ys
                                      - Hy->data[i] * Hy->data[j] / yHy;
                }
            }
        }

        vector_copy(x_new, x);
        vector_copy(g_new, g);
        vector_free(d);
    }

    result->iterations = iter;
    result->f_value = problem->f.eval(&problem->f, x);
    result->grad_norm = vector_norm_l2(g);
    result->status = (iter < max_iters) ? OPT_SUCCESS : OPT_MAXITER_REACHED;

    if (!result->x.data) {
        result->x.data = (double*)malloc(n * sizeof(double));
        result->x.n = n;
    }
    memcpy(result->x.data, x->data, n * sizeof(double));

    vector_free(x); vector_free(x_new);
    vector_free(g); vector_free(g_new);
    vector_free(s); vector_free(y); vector_free(Hy);
    matrix_free(H);
}

/* ═══════════════════════════════════════════════════════════════════════
 * SR1 Method (Symmetric Rank-One Update)
 *
 * Update: H_{k+1} = H_k + (s_k - H_k y_k)(s_k - H_k y_k)^T
 *                         / ((s_k - H_k y_k)^T y_k)
 *
 * Unlike BFGS/DFP, SR1 does NOT guarantee positive definiteness.
 * This is both a weakness (can lose descent property) and a strength
 * (can better approximate indefinite Hessians, useful for non-convex
 * problems and constrained optimization where the reduced Hessian
 * may be indefinite).
 *
 * SR1 has the "hereditary" property: for quadratic functions with
 * exact line search, it generates the true Hessian after n steps
 * (same as BFGS/DFP).
 *
 * Reference: Broyden (1967), "Quasi-Newton methods and their
 *            application to function minimization"
 * ═══════════════════════════════════════════════════════════════════════ */

void sr1_method(const unconstr_problem_t *problem,
                size_t max_iters, opt_result_t *result)
{
    assert(problem && result);
    size_t n = problem->n;

    vector_t *x = vector_alloc(n), *x_new = vector_alloc(n);
    vector_t *g = vector_alloc(n), *g_new = vector_alloc(n);
    vector_t *s = vector_alloc(n), *y = vector_alloc(n);
    vector_t *r = vector_alloc(n);  /* Residual: s - H y */
    matrix_t *H = matrix_alloc(n, n);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            MAT_ELT(H, i, j) = (i == j) ? 1.0 : 0.0;
        }
    }

    vector_copy(&problem->x0, x);
    problem->f.grad(&problem->f, x, g);

    backtracking_params_t ls_params;
    ls_params.alpha_init = 1.0; ls_params.rho = 0.5;
    ls_params.c1 = 1e-4; ls_params.c2 = 0.9; ls_params.max_iter = 30;

    size_t iter;
    for (iter = 0; iter < max_iters; iter++) {
        if (vector_norm_l2(g) < GRAD_TOL) break;

        vector_t *d = vector_alloc(n);
        matrix_vector_mul(H, g, d);
        for (size_t i = 0; i < n; i++) d->data[i] = -d->data[i];

        double alpha = backtracking_linesearch(&problem->f, x, d, g, &ls_params);

        for (size_t i = 0; i < n; i++) {
            x_new->data[i] = x->data[i] + alpha * d->data[i];
        }

        problem->f.grad(&problem->f, x_new, g_new);

        for (size_t i = 0; i < n; i++) {
            s->data[i] = x_new->data[i] - x->data[i];
            y->data[i] = g_new->data[i] - g->data[i];
        }

        /* r = s - H*y */
        matrix_vector_mul(H, y, r);
        for (size_t i = 0; i < n; i++) {
            r->data[i] = s->data[i] - r->data[i];
        }

        double rTy = vector_dot(r, y);

        /* Only update if denominator is non-negligible */
        if (fabs(rTy) > 1e-12) {
            /* SR1 update: H += r * r^T / r^T y */
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    MAT_ELT(H, i, j) += r->data[i] * r->data[j] / rTy;
                }
            }
        }

        vector_copy(x_new, x);
        vector_copy(g_new, g);
        vector_free(d);
    }

    result->iterations = iter;
    result->f_value = problem->f.eval(&problem->f, x);
    result->grad_norm = vector_norm_l2(g);
    result->status = (iter < max_iters) ? OPT_SUCCESS : OPT_MAXITER_REACHED;

    if (!result->x.data) {
        result->x.data = (double*)malloc(n * sizeof(double));
        result->x.n = n;
    }
    memcpy(result->x.data, x->data, n * sizeof(double));

    vector_free(x); vector_free(x_new);
    vector_free(g); vector_free(g_new);
    vector_free(s); vector_free(y); vector_free(r);
    matrix_free(H);
}

/* ═══════════════════════════════════════════════════════════════════════
 * L-BFGS (Limited-Memory BFGS)
 *
 * Stores only the last m pairs (s_k, y_k) instead of the full n×n
 * Hessian approximation. The matrix-vector product H_k * v is
 * computed implicitly via the two-loop recursion:
 *
 *   q = v
 *   for i = k-1, ..., k-m:
 *     alpha_i = rho_i * s_i^T q
 *     q = q - alpha_i * y_i
 *   r = H_k^0 * q
 *   for i = k-m, ..., k-1:
 *     beta = rho_i * y_i^T r
 *     r = r + s_i * (alpha_i - beta)
 *   return r
 *
 * Memory: O(mn) vs O(n^2) for full BFGS.
 * This makes L-BFGS applicable to problems with millions of variables.
 *
 * Theorem (Liu & Nocedal, 1989):
 * L-BFGS with m >= 3 converges globally and R-linearly for convex
 * problems with Wolfe line search.
 *
 * L-BFGS is the de facto standard for large-scale unconstrained
 * optimization in machine learning (e.g., training logistic regression,
 * CRFs, neural network fine-tuning).
 *
 * Reference: Nocedal (1980), "Updating quasi-Newton matrices with
 *            limited storage"; Liu & Nocedal (1989)
 * ═══════════════════════════════════════════════════════════════════════ */

void lbfgs_method(const unconstr_problem_t *problem,
                  size_t m, size_t max_iters,
                  opt_result_t *result)
{
    assert(problem && result && m >= 1);
    size_t n = problem->n;

    vector_t *x = vector_alloc(n);
    vector_t *g = vector_alloc(n);
    vector_t *g_new = vector_alloc(n);
    vector_t *d = vector_alloc(n);

    /* Storage for the last m (s, y) pairs */
    vector_t **S = (vector_t**)malloc(m * sizeof(vector_t*));
    vector_t **Y = (vector_t**)malloc(m * sizeof(vector_t*));
    double *rho = (double*)malloc(m * sizeof(double));
    for (size_t i = 0; i < m; i++) {
        S[i] = vector_alloc(n);
        Y[i] = vector_alloc(n);
    }
    size_t n_stored = 0;  /* Number of pairs currently stored */
    size_t ptr = 0;       /* Circular buffer pointer */

    vector_copy(&problem->x0, x);
    problem->f.grad(&problem->f, x, g);

    backtracking_params_t ls_params;
    ls_params.alpha_init = 1.0; ls_params.rho = 0.5;
    ls_params.c1 = 1e-4; ls_params.c2 = 0.9; ls_params.max_iter = 30;

    size_t iter;
    for (iter = 0; iter < max_iters; iter++) {
        if (vector_norm_l2(g) < GRAD_TOL) break;

        /* Two-loop recursion to compute d = -H_k * g_k */
        double *alpha_arr = (double*)malloc(m * sizeof(double));
        vector_t *q = vector_alloc(n);
        vector_copy(g, q);

        /* First loop: unrolling */
        size_t cnt = (n_stored < m) ? n_stored : m;
        for (size_t i = 0; i < cnt; i++) {
            size_t idx = (ptr - 1 - i + m) % m;
            alpha_arr[i] = rho[idx] * vector_dot(S[idx], q);
            for (size_t j = 0; j < n; j++) {
                q->data[j] -= alpha_arr[i] * Y[idx]->data[j];
            }
        }

        /* Scale: H_k^0 = gamma * I, gamma = s_{k-1}^T y_{k-1} / y_{k-1}^T y_{k-1} */
        double gamma = 1.0;
        if (n_stored > 0) {
            size_t last = (ptr - 1 + m) % m;
            double sy = vector_dot(S[last], Y[last]);
            gamma = (fabs(sy) > 1e-12) ? sy / (vector_dot(Y[last], Y[last]) + 1e-12) : 1.0;
            if (gamma < 1e-4 || gamma > 1e4) gamma = 1.0;
        }

        /* Scale q */
        for (size_t j = 0; j < n; j++) d->data[j] = gamma * q->data[j];

        /* Second loop: rolling */
        for (size_t i = 0; i < cnt; i++) {
            size_t ii = cnt - 1 - i;
            size_t idx = (ptr - 1 - ii + m) % m;
            double beta = rho[idx] * vector_dot(Y[idx], d);
            for (size_t j = 0; j < n; j++) {
                d->data[j] += S[idx]->data[j] * (alpha_arr[ii] - beta);
            }
        }

        /* Make it a descent direction */
        for (size_t j = 0; j < n; j++) d->data[j] = -d->data[j];

        vector_free(q);
        free(alpha_arr);

        /* Line search */
        double alpha = backtracking_linesearch(&problem->f, x, d, g, &ls_params);

        /* Update x_new = x + alpha * d */
        vector_t *x_new = vector_alloc(n);
        for (size_t i = 0; i < n; i++) {
            x_new->data[i] = x->data[i] + alpha * d->data[i];
        }

        problem->f.grad(&problem->f, x_new, g_new);

        /* Store (s, y) pair */
        for (size_t i = 0; i < n; i++) {
            S[ptr]->data[i] = x_new->data[i] - x->data[i];
            Y[ptr]->data[i] = g_new->data[i] - g->data[i];
        }
        double ys = vector_dot(Y[ptr], S[ptr]);
        if (ys > 1e-16) {
            rho[ptr] = 1.0 / ys;
            ptr = (ptr + 1) % m;
            if (n_stored < m) n_stored++;
        }
        /* If ys <= 0, skip storing this pair */

        vector_copy(x_new, x);
        vector_copy(g_new, g);
        vector_free(x_new);
    }

    result->iterations = iter;
    result->f_value = problem->f.eval(&problem->f, x);
    result->grad_norm = vector_norm_l2(g);
    result->status = (iter < max_iters) ? OPT_SUCCESS : OPT_MAXITER_REACHED;

    if (!result->x.data) {
        result->x.data = (double*)malloc(n * sizeof(double));
        result->x.n = n;
    }
    memcpy(result->x.data, x->data, n * sizeof(double));

    vector_free(x); vector_free(g); vector_free(g_new); vector_free(d);
    for (size_t i = 0; i < m; i++) {
        vector_free(S[i]); vector_free(Y[i]);
    }
    free(S); free(Y); free(rho);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Trust-Region Newton Method (Levenberg-Marquardt style)
 *
 * At each iteration: minimize the quadratic model within a trust region:
 *     m_k(d) = f_k + g_k^T d + (1/2) d^T B_k d
 *     s.t. ||d|| <= Delta_k
 *
 * The trust-region radius Delta_k is adjusted based on the ratio:
 *     rho_k = (f(x_k) - f(x_k+d_k)) / (m_k(0) - m_k(d_k))
 *
 * If rho_k ≈ 1: model agrees well => increase Delta
 * If rho_k << 1: model disagrees => decrease Delta and reject step
 *
 * The subproblem is solved approximately:
 *     (B + lambda I) d = -g
 * with lambda chosen so ||d(lambda)|| ≈ Delta (if constraint active).
 *
 * This regularization (lambda I) is the Levenberg-Marquardt
 * stabilization, making the method robust when B is nearly
 * singular or indefinite.
 *
 * Reference: Conn, Gould, Toint (2000), "Trust-Region Methods"
 * ═══════════════════════════════════════════════════════════════════════ */

void trust_region_newton(const unconstr_problem_t *problem,
                         const newton_config_t *config,
                         opt_result_t *result)
{
    assert(problem && config && result && problem->f.hess);
    size_t n = problem->n;

    vector_t *x = vector_alloc(n);
    vector_t *g = vector_alloc(n);
    vector_t *d = vector_alloc(n);
    vector_t *x_trial = vector_alloc(n);
    matrix_t *B = matrix_alloc(n, n);
    matrix_t *L = matrix_alloc(n, n);

    vector_copy(&problem->x0, x);

    double Delta = 1.0;           /* Trust-region radius */
    double eta = 0.1;             /* Threshold for accepting step */
    size_t iter;

    for (iter = 0; iter < config->max_iters; iter++) {
        problem->f.grad(&problem->f, x, g);
        double gnorm = vector_norm_l2(g);
        if (gnorm < config->grad_tol) break;

        problem->f.hess(&problem->f, x, B);

        /* Approximate solution: (B + lambda I) d = -g */
        double lambda = 0.0;
        int step_ok = 0;

        for (int trial = 0; trial < 10 && !step_ok; trial++) {
            /* Copy B + lambda I to L for Cholesky */
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    MAT_ELT(L, i, j) = MAT_ELT(B, i, j);
                }
                MAT_ELT(L, i, i) += lambda;
            }

            if (cholesky_decomp(L) == 0) {
                for (size_t i = 0; i < n; i++) d->data[i] = -g->data[i];
                cholesky_solve(L, d, d);
                double dnorm = vector_norm_l2(d);

                if (dnorm <= 1.1 * Delta) {
                    step_ok = 1;
                } else {
                    /* Increase lambda to reduce step size */
                    lambda += Delta * (dnorm - Delta) / Delta;
                    if (lambda < 0.0) lambda = 1e-8;
                }
            } else {
                /* Not PD: increase lambda */
                lambda = (lambda == 0.0) ? 1e-4 : lambda * 10.0;
            }
        }

        if (!step_ok) {
            /* Fallback: Cauchy point within trust region */
            double gTg = vector_dot(g, g);
            double gTBg = 0.0;
            vector_t *Bg = vector_alloc(n);
            matrix_vector_mul(B, g, Bg);
            gTBg = vector_dot(g, Bg);
            vector_free(Bg);

            double tau = 1.0;
            if (gTBg > 1e-16) {
                tau = gTg / gTBg;
            }
            if (tau * vector_norm_l2(g) > Delta) {
                tau = Delta / gnorm;
            }
            for (size_t i = 0; i < n; i++) d->data[i] = -tau * g->data[i];
        }

        /* Evaluate actual vs predicted reduction */
        double f_curr = problem->f.eval(&problem->f, x);

        for (size_t i = 0; i < n; i++) {
            x_trial->data[i] = x->data[i] + d->data[i];
        }
        double f_trial = problem->f.eval(&problem->f, x_trial);

        double ared = f_curr - f_trial;  /* Actual reduction */
        double pred = -(vector_dot(g, d) + 0.5 * quadratic_form(B, d));

        double rho = (fabs(pred) > 1e-16) ? (ared / pred) : 0.0;

        /* Update trust-region radius */
        if (rho < 0.25) {
            Delta *= 0.25;
        } else if (rho > 0.75 && vector_norm_l2(d) > 0.9 * Delta) {
            Delta *= 2.0;
        }

        /* Accept or reject step */
        if (rho > eta) {
            vector_copy(x_trial, x);
        }
        /* If rho <= eta: reject step, keep x, retry with smaller Delta */
    }

    result->iterations = iter;
    result->f_value = problem->f.eval(&problem->f, x);
    result->grad_norm = vector_norm_l2(g);
    result->status = (iter < config->max_iters) ? OPT_SUCCESS
                                                 : OPT_MAXITER_REACHED;

    if (!result->x.data) {
        result->x.data = (double*)malloc(n * sizeof(double));
        result->x.n = n;
    }
    memcpy(result->x.data, x->data, n * sizeof(double));

    vector_free(x); vector_free(g); vector_free(d); vector_free(x_trial);
    matrix_free(B); matrix_free(L);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Gauss-Newton Method for Nonlinear Least Squares
 *
 * Problem: minimize f(x) = (1/2) ||r(x)||^2
 * where r: R^n -> R^m is the residual vector.
 *
 * Gauss-Newton approximation of the Hessian:
 *     H_GN = J^T J    (where J = Jacobian of r)
 *
 * This ignores the second-order term sum_i r_i * Hess(r_i), which
 * is small when residuals are small (good fit) or when the Hessians
 * of individual residuals tend to cancel.
 *
 * The Gauss-Newton step solves: J^T J d_GN = -J^T r
 *
 * Convergence: Linear for small-residual problems; can diverge for
 * large-residual problems. Levenberg-Marquardt adds regularization
 * for robustness.
 *
 * This is the workhorse for parameter estimation, curve fitting,
 * and many computer vision problems (bundle adjustment).
 * ═══════════════════════════════════════════════════════════════════════ */

void gauss_newton_method(size_t n_vars, size_t n_resid,
                         void (*r)(const vector_t *x, vector_t *r_out),
                         void (*J)(const vector_t *x, matrix_t *J_out),
                         const vector_t *x0, size_t max_iters,
                         opt_result_t *result)
{
    assert(r && J && x0 && result);
    size_t n = n_vars, m = n_resid;

    vector_t *x = vector_alloc(n);
    vector_t *resid = vector_alloc(m);
    vector_t *g = vector_alloc(n);    /* J^T r */
    matrix_t *Jac = matrix_alloc(m, n);
    matrix_t *JTJ = matrix_alloc(n, n);
    vector_t *d = vector_alloc(n);

    vector_copy(x0, x);

    size_t iter;
    for (iter = 0; iter < max_iters; iter++) {
        /* Compute residuals and Jacobian */
        r(x, resid);
        J(x, Jac);

        /* Compute gradient: g = J^T r */
        for (size_t i = 0; i < n; i++) {
            g->data[i] = 0.0;
            for (size_t k = 0; k < m; k++) {
                g->data[i] += MAT_ELT(Jac, k, i) * resid->data[k];
            }
        }

        if (vector_norm_l2(g) < GRAD_TOL) break;

        /* Compute J^T J */
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                double sum = 0.0;
                for (size_t k = 0; k < m; k++) {
                    sum += MAT_ELT(Jac, k, i) * MAT_ELT(Jac, k, j);
                }
                MAT_ELT(JTJ, i, j) = sum;
            }
        }

        /* Solve J^T J d = -g */
        matrix_t *L = matrix_alloc(n, n);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                MAT_ELT(L, i, j) = MAT_ELT(JTJ, i, j);
            }
        }

        for (size_t i = 0; i < n; i++) d->data[i] = -g->data[i];

        if (cholesky_decomp(L) == 0) {
            cholesky_solve(L, d, d);
        }
        /* If not PD, d remains as -g (gradient descent fallback) */

        matrix_free(L);

        /* Step with line search */
        double alpha = 1.0;
        double f_curr = 0.5 * vector_dot(resid, resid);

        for (int ls = 0; ls < 20; ls++) {
            vector_t *x_trial = vector_alloc(n);
            for (size_t i = 0; i < n; i++) {
                x_trial->data[i] = x->data[i] + alpha * d->data[i];
            }
            r(x_trial, resid);
            double f_trial = 0.5 * vector_dot(resid, resid);
            vector_free(x_trial);

            if (f_trial < f_curr) break;
            alpha *= 0.5;
        }

        for (size_t i = 0; i < n; i++) {
            x->data[i] += alpha * d->data[i];
        }
    }

    r(x, resid);
    result->f_value = 0.5 * vector_dot(resid, resid);
    result->grad_norm = vector_norm_l2(g);
    result->iterations = iter;
    result->status = (iter < max_iters) ? OPT_SUCCESS : OPT_MAXITER_REACHED;

    if (!result->x.data) {
        result->x.data = (double*)malloc(n * sizeof(double));
        result->x.n = n;
    }
    memcpy(result->x.data, x->data, n * sizeof(double));

    vector_free(x); vector_free(resid); vector_free(g);
    matrix_free(Jac); matrix_free(JTJ); vector_free(d);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Levenberg-Marquardt Method
 *
 * Interpolates between Gauss-Newton (lambda -> 0) and gradient
 * descent (lambda -> infinity):
 *     (J^T J + lambda I) d = -J^T r
 *
 * This provides robust global convergence while retaining fast
 * local convergence near the optimum for small-residual problems.
 *
 * Algorithm (Marquardt, 1963):
 * 1. Start with moderate lambda (e.g., 1e-3)
 * 2. If the step reduces the objective: accept, decrease lambda
 * 3. If the step increases the objective: reject, increase lambda
 *
 * The damping term lambda I adds positive definiteness to the
 * approximate Hessian, ensuring a descent direction. As lambda
 * increases, the step direction approaches negative gradient
 * (guaranteed descent); as lambda -> 0, the step approaches
 * Gauss-Newton (fast near optimum).
 *
 * This is the standard method for nonlinear least squares in
 * practice (used in Ceres Solver, scipy.optimize.least_squares,
 * MINPACK's lmder).
 *
 * Reference: Levenberg (1944), Marquardt (1963);
 *            More (1978), "The Levenberg-Marquardt algorithm:
 *            implementation and theory"
 * ═══════════════════════════════════════════════════════════════════════ */

void levenberg_marquardt(size_t n_vars, size_t n_resid,
                         void (*r)(const vector_t *x, vector_t *r_out),
                         void (*J)(const vector_t *x, matrix_t *J_out),
                         const vector_t *x0, double lambda0,
                         size_t max_iters, opt_result_t *result)
{
    assert(r && J && x0 && result);
    size_t n = n_vars, m = n_resid;

    vector_t *x = vector_alloc(n), *resid = vector_alloc(m);
    vector_t *g = vector_alloc(n), *d = vector_alloc(n);
    vector_t *x_trial = vector_alloc(n), *resid_trial = vector_alloc(m);
    matrix_t *Jac = matrix_alloc(m, n);
    matrix_t *JTJ = matrix_alloc(n, n), *L = matrix_alloc(n, n);

    vector_copy(x0, x);
    r(x, resid);
    double f_curr = 0.5 * vector_dot(resid, resid);
    double lambda = lambda0;

    size_t iter;
    for (iter = 0; iter < max_iters; iter++) {
        r(x, resid);
        J(x, Jac);

        /* Gradient: g = J^T r */
        for (size_t i = 0; i < n; i++) {
            g->data[i] = 0.0;
            for (size_t k = 0; k < m; k++) {
                g->data[i] += MAT_ELT(Jac, k, i) * resid->data[k];
            }
        }

        if (vector_norm_l2(g) < GRAD_TOL) break;

        /* J^T J */
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                double sum = 0.0;
                for (size_t k = 0; k < m; k++) {
                    sum += MAT_ELT(Jac, k, i) * MAT_ELT(Jac, k, j);
                }
                MAT_ELT(JTJ, i, j) = sum;
            }
        }

        /* Solve (J^T J + lambda I) d = -g */
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                MAT_ELT(L, i, j) = MAT_ELT(JTJ, i, j);
            }
            MAT_ELT(L, i, i) += lambda;
            d->data[i] = -g->data[i];
        }

        if (cholesky_decomp(L) == 0) {
            cholesky_solve(L, d, d);
        } else {
            /* Fallback: gradient descent */
            double gnorm = vector_norm_l2(g);
            if (gnorm > 1e-12) {
                for (size_t i = 0; i < n; i++) d->data[i] = -g->data[i] / gnorm;
            }
        }

        /* Trial step */
        for (size_t i = 0; i < n; i++) {
            x_trial->data[i] = x->data[i] + d->data[i];
        }
        r(x_trial, resid_trial);
        double f_trial = 0.5 * vector_dot(resid_trial, resid_trial);

        /* Predicted reduction from linear model */
        double pred = -(vector_dot(g, d) + 0.5 * vector_dot(d, g));
        /* Actually: pred = -(g^T d + 1/2 d^T J^T J d), and J^T J d ≈ -(g + lambda d)*/

        double rho = (f_curr - f_trial) / (pred + 1e-16);

        if (rho > 0.0) {
            /* Accept step */
            vector_copy(x_trial, x);
            f_curr = f_trial;
            lambda *= 0.1;  /* Decrease damping */
            if (lambda < 1e-8) lambda = 1e-8;
        } else {
            /* Reject step */
            lambda *= 10.0;  /* Increase damping */
            if (lambda > 1e8) lambda = 1e8;
        }
    }

    result->f_value = f_curr;
    result->grad_norm = vector_norm_l2(g);
    result->iterations = iter;
    result->status = (iter < max_iters) ? OPT_SUCCESS : OPT_MAXITER_REACHED;

    if (!result->x.data) {
        result->x.data = (double*)malloc(n * sizeof(double));
        result->x.n = n;
    }
    memcpy(result->x.data, x->data, n * sizeof(double));

    vector_free(x); vector_free(resid); vector_free(g); vector_free(d);
    vector_free(x_trial); vector_free(resid_trial);
    matrix_free(Jac); matrix_free(JTJ); matrix_free(L);
}
