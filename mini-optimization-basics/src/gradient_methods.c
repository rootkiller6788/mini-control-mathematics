/**
 * gradient_methods.c — First-Order Optimization Methods Implementation
 *
 * Knowledge points:
 * L5: Gradient descent with backtracking line search (Armijo condition)
 * L5: Conjugate gradient (linear CG + nonlinear Fletcher-Reeves)
 * L3: Lipschitz constant estimation via finite differences
 * L8: Nesterov accelerated gradient (optimal O(1/k^2) rate)
 * L5: Heavy ball method (Polyak momentum)
 * L5: Subgradient method for nondifferentiable convex functions
 * L5: Coordinate descent with exact coordinate line search
 * L8: Proximal gradient (ISTA/FISTA)
 *
 * Reference: Nocedal & Wright (2006), Boyd & Vandenberghe (2004)
 */

#include "gradient_methods.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>

/* ═══════════════════════════════════════════════════════════════════════
 * Configuration Defaults
 * ═══════════════════════════════════════════════════════════════════════ */

gradient_descent_config_t gd_config_default(void)
{
    gradient_descent_config_t cfg;
    cfg.step_size = 1.0;
    cfg.ls_type = LS_BACKTRACKING;
    cfg.ls_params.alpha_init = 1.0;
    cfg.ls_params.rho = 0.5;
    cfg.ls_params.c1 = 1e-4;
    cfg.ls_params.c2 = 0.9;
    cfg.ls_params.max_iter = 50;
    cfg.grad_tol = GRAD_TOL;
    cfg.max_iters = OPT_MAX_ITERS;
    cfg.use_nesterov = 0;
    cfg.momentum_beta = 0.0;
    cfg.verbose = 0;
    return cfg;
}

cg_config_t cg_config_default(void)
{
    cg_config_t cfg;
    cfg.max_iters = 1000;
    cfg.grad_tol = GRAD_TOL;
    cfg.variant = 0;
    cfg.verbose = 0;
    return cfg;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Backtracking Line Search (Armijo Condition)
 *
 * Armijo sufficient decrease condition:
 *     f(x + alpha d) <= f(x) + c1 * alpha * grad_f(x)^T d
 *
 * Theorem (Armijo 1966, Wolfe 1969):
 * For C^1 f with Lipschitz gradient, backtracking with 0 < c1 < 1
 * and 0 < rho < 1 terminates in finitely many steps with a step
 * satisfying sufficient decrease.
 *
 * Each iteration: one function evaluation. Complexity:
 * O(log(alpha_min)/log(rho)) evaluations worst case.
 * ═══════════════════════════════════════════════════════════════════════ */

double backtracking_linesearch(const convex_function_t *f,
                                const vector_t *x,
                                const vector_t *d,
                                const vector_t *grad_at_x,
                                const backtracking_params_t *params)
{
    assert(f && x && d && grad_at_x && params);
    assert(f->eval && "Function must have eval oracle");

    double dir_deriv = vector_dot(grad_at_x, d);

    /* d must be a descent direction */
    if (dir_deriv >= 0.0) {
        return 0.0;
    }

    double alpha = params->alpha_init;
    double fx = f->eval(f, x);
    size_t n = x->n;
    vector_t *x_trial = vector_alloc(n);

    for (size_t iter = 0; iter < params->max_iter; iter++) {
        for (size_t i = 0; i < n; i++) {
            x_trial->data[i] = x->data[i] + alpha * d->data[i];
        }

        double f_trial = f->eval(f, x_trial);

        /* Armijo condition */
        if (f_trial <= fx + params->c1 * alpha * dir_deriv) {
            vector_free(x_trial);
            return alpha;
        }

        alpha *= params->rho;
        if (alpha < 1e-16) break;
    }

    vector_free(x_trial);
    return alpha;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Gradient Descent with Backtracking Line Search
 *
 * Algorithm: Cauchy steepest descent (1847) with Armijo globalization.
 *
 * Iteration: x_{k+1} = x_k - alpha_k grad f(x_k)
 * where alpha_k satisfies the Armijo condition.
 *
 * Convergence (convex, L-smooth):
 *     f(x_k) - f* <= 2L ||x_0 - x*||^2 / k
 * This is the fundamental sublinear O(1/k) rate.
 *
 * Convergence (mu-strongly convex, L-smooth):
 *     f(x_k) - f* <= (1 - mu/L)^k * (f(x_0) - f*)
 * Linear convergence with condition number kappa = L/mu.
 *
 * Per-iteration cost: O(n) for gradient, O(n) for vector ops.
 * Total: O(n/epsilon) for convex, O(n * log(1/eps) * kappa) for strongly convex.
 * ═══════════════════════════════════════════════════════════════════════ */

void gradient_descent(const unconstr_problem_t *problem,
                      const gradient_descent_config_t *config,
                      opt_result_t *result)
{
    assert(problem && config && result);
    assert(problem->f.eval && problem->f.grad);

    size_t n = problem->n;
    vector_t *x = vector_alloc(n);
    vector_t *g = vector_alloc(n);
    vector_t *d = vector_alloc(n);
    vector_t *x_prev = NULL;

    vector_copy(&problem->x0, x);

    if (config->use_nesterov) {
        x_prev = vector_alloc(n);
        vector_copy(x, x_prev);
    }

    double t_prev = 1.0;
    size_t iter;

    for (iter = 0; iter < config->max_iters; iter++) {
        problem->f.grad(&problem->f, x, g);
        double gnorm = vector_norm_l2(g);

        if (gnorm < config->grad_tol) break;

        /* Nesterov extrapolation */
        if (config->use_nesterov && x_prev) {
            double t_new = 0.5 * (1.0 + sqrt(1.0 + 4.0 * t_prev * t_prev));
            double beta = (t_prev - 1.0) / t_new;
            t_prev = t_new;

            vector_copy(x, x_prev);
            for (size_t i = 0; i < n; i++) {
                x->data[i] = x->data[i] + beta * (x->data[i] - x_prev->data[i]);
            }
            problem->f.grad(&problem->f, x, g);
        }

        /* Search direction */
        for (size_t i = 0; i < n; i++) {
            d->data[i] = -g->data[i];
        }

        /* Momentum */
        if (config->momentum_beta > 0.0 && x_prev) {
            for (size_t i = 0; i < n; i++) {
                d->data[i] += config->momentum_beta * (x->data[i] - x_prev->data[i]);
            }
            vector_copy(x, x_prev);
        }

        /* Step size selection */
        double alpha;
        switch (config->ls_type) {
        case LS_FIXED:
            alpha = config->step_size;
            break;
        case LS_BACKTRACKING:
            alpha = backtracking_linesearch(&problem->f, x, d, g,
                                            &config->ls_params);
            break;
        case LS_DIMINISHING:
            alpha = config->step_size / (1.0 + (double)iter);
            break;
        default:
            alpha = backtracking_linesearch(&problem->f, x, d, g,
                                            &config->ls_params);
        }

        if (alpha <= 0.0) {
            result->status = OPT_LINE_SEARCH_FAILED;
            break;
        }

        /* x_{k+1} = x_k + alpha d_k */
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
    if (x_prev) vector_free(x_prev);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Exact Line Search for Quadratic Functions
 *
 * For f(x) = (1/2) x^T Q x - b^T x (Q SPD):
 * The exact minimizer along -grad is:
 *     alpha* = (grad^T grad) / (grad^T Q grad)
 *
 * This yields the Cauchy point. The convergence rate reveals the
 * fundamental limitation of first-order methods on ill-conditioned
 * problems:
 *     ||e_{k+1}||_Q <= ((kappa-1)/(kappa+1)) ||e_k||_Q
 * where kappa = lambda_max / lambda_min.
 *
 * For kappa >> 1, (kappa-1)/(kappa+1) ≈ 1 - 2/kappa → very slow.
 * ═══════════════════════════════════════════════════════════════════════ */

double steepest_descent_exact_step(const convex_function_t *f,
                                    vector_t *x, vector_t *g)
{
    assert(f && x && g && f->hess);

    size_t n = x->n;
    matrix_t *Q = matrix_alloc(n, n);
    vector_t *Qg = vector_alloc(n);

    f->hess(f, x, Q);
    matrix_vector_mul(Q, g, Qg);

    double gTg = vector_dot(g, g);
    double gTQg = vector_dot(g, Qg);
    double alpha = (gTQg > 1e-16) ? (gTg / gTQg) : 1e-4;

    for (size_t i = 0; i < n; i++) {
        x->data[i] -= alpha * g->data[i];
    }

    matrix_free(Q); vector_free(Qg);
    return alpha;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Linear Conjugate Gradient Method
 *
 * Solves Q x = b for SPD Q. Generates Q-orthogonal search directions
 * from the Krylov subspace K_k(Q, r_0).
 *
 * Theorem (Hestenes & Stiefel, 1952):
 * CG finds the exact solution in at most n iterations (exact arithmetic).
 * Error bound after k iterations:
 *     ||x_k - x*||_Q <= 2 * ((sqrt(kappa)-1)/(sqrt(kappa)+1))^k * ||x_0 - x*||_Q
 *
 * The sqrt(kappa) dependence makes CG vastly superior to steepest
 * descent for ill-conditioned linear systems. CG is the foundation
 * of iterative methods for large sparse SPD systems.
 * ═══════════════════════════════════════════════════════════════════════ */

void conjugate_gradient_linear(const matrix_t *Q, const vector_t *b,
                                vector_t *x, const cg_config_t *config)
{
    assert(Q && b && x && config);
    size_t n = Q->rows;
    assert(Q->cols == n && b->n == n && x->n == n);

    vector_t *r = vector_alloc(n);
    vector_t *p = vector_alloc(n);
    vector_t *Ap = vector_alloc(n);

    /* r_0 = b - Q x_0 */
    matrix_vector_mul(Q, x, Ap);
    for (size_t i = 0; i < n; i++) {
        r->data[i] = b->data[i] - Ap->data[i];
        p->data[i] = r->data[i];
    }

    double rsold = vector_dot(r, r);

    for (size_t k = 0; k < config->max_iters && k < n; k++) {
        if (sqrt(rsold) < config->grad_tol) break;

        matrix_vector_mul(Q, p, Ap);
        double pAp = vector_dot(p, Ap);
        if (pAp <= 1e-16) break;

        double alpha = rsold / pAp;

        for (size_t i = 0; i < n; i++) {
            x->data[i] += alpha * p->data[i];
        }
        for (size_t i = 0; i < n; i++) {
            r->data[i] -= alpha * Ap->data[i];
        }

        double rsnew = vector_dot(r, r);
        double beta = rsnew / rsold;

        for (size_t i = 0; i < n; i++) {
            p->data[i] = r->data[i] + beta * p->data[i];
        }
        rsold = rsnew;
    }

    vector_free(r); vector_free(p); vector_free(Ap);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Nonlinear Conjugate Gradient (Fletcher-Reeves variant)
 *
 * Extends CG to general convex functions. The Fletcher-Reeves
 * beta formula: beta_k = ||g_{k+1}||^2 / ||g_k||^2
 *
 * Unlike linear CG, the search directions are no longer exactly
 * conjugate, but the method preserves the low-memory advantage
 * (no Hessian storage). Combined with strong Wolfe line search,
 * global convergence is guaranteed.
 *
 * Reference: Fletcher & Reeves, "Function minimization by conjugate
 *            gradients", The Computer Journal (1964).
 * ═══════════════════════════════════════════════════════════════════════ */

void nonlinear_cg(const unconstr_problem_t *problem,
                  const cg_config_t *config, opt_result_t *result)
{
    assert(problem && config && result);
    size_t n = problem->n;

    vector_t *x = vector_alloc(n);
    vector_t *g = vector_alloc(n);
    vector_t *g_prev = vector_alloc(n);
    vector_t *d = vector_alloc(n);

    vector_copy(&problem->x0, x);
    problem->f.grad(&problem->f, x, g);
    double gnorm = vector_norm_l2(g);

    if (gnorm < config->grad_tol) {
        result->status = OPT_SUCCESS;
        result->iterations = 0;
        result->f_value = problem->f.eval(&problem->f, x);
        result->grad_norm = gnorm;
        goto ncg_done;
    }

    for (size_t i = 0; i < n; i++) {
        d->data[i] = -g->data[i];
        g_prev->data[i] = g->data[i];
    }

    backtracking_params_t ls_params;
    ls_params.alpha_init = 1.0;
    ls_params.rho = 0.5;
    ls_params.c1 = 1e-4;
    ls_params.c2 = 0.9;
    ls_params.max_iter = 50;

    size_t iter;
    for (iter = 0; iter < config->max_iters; iter++) {
        double alpha = backtracking_linesearch(&problem->f, x, d, g, &ls_params);
        if (alpha <= 0.0) break;

        for (size_t i = 0; i < n; i++) {
            x->data[i] += alpha * d->data[i];
        }

        problem->f.grad(&problem->f, x, g);
        gnorm = vector_norm_l2(g);

        if (gnorm < config->grad_tol) { iter++; break; }

        /* Fletcher-Reeves beta */
        double gprev_n2 = vector_dot(g_prev, g_prev);
        double g_n2 = vector_dot(g, g);
        double beta = g_n2 / gprev_n2;

        for (size_t i = 0; i < n; i++) {
            d->data[i] = -g->data[i] + beta * d->data[i];
            g_prev->data[i] = g->data[i];
        }
    }

    result->status = (gnorm < config->grad_tol) ? OPT_SUCCESS
                                                 : OPT_MAXITER_REACHED;
    result->iterations = iter;
    result->f_value = problem->f.eval(&problem->f, x);
    result->grad_norm = gnorm;

ncg_done:
    if (!result->x.data) {
        result->x.data = (double*)malloc(n * sizeof(double));
        result->x.n = n;
    }
    memcpy(result->x.data, x->data, n * sizeof(double));

    vector_free(x); vector_free(g);
    vector_free(g_prev); vector_free(d);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Lipschitz Constant Estimation
 *
 * For L-smooth f: ||grad f(x) - grad f(y)|| <= L ||x - y||
 *
 * Estimates L by sampling random pairs and computing the maximum
 * ratio of gradient difference to point difference.
 * The estimated L is a lower bound; the true L >= L_est.
 *
 * Optimal fixed step: alpha <= 2/L for gradient descent,
 * alpha <= 1/L for proximal gradient.
 *
 * Reference: Beck & Teboulle (2009), "Gradient-based algorithms
 *            with applications to signal recovery".
 * ═══════════════════════════════════════════════════════════════════════ */

double estimate_lipschitz(const convex_function_t *f, size_t n, size_t nsamples)
{
    assert(f && f->grad && n > 0);
    if (nsamples == 0) nsamples = 20;

    vector_t *x = vector_alloc(n);
    vector_t *y = vector_alloc(n);
    vector_t *gx = vector_alloc(n);
    vector_t *gy = vector_alloc(n);

    double L_est = 0.0;

    for (size_t s = 0; s < nsamples; s++) {
        /* Random points in [-1, 1]^n */
        for (size_t i = 0; i < n; i++) {
            x->data[i] = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
            y->data[i] = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
        }

        f->grad(f, x, gx);
        f->grad(f, y, gy);

        double diff_n2 = 0.0, grad_diff_n2 = 0.0;
        for (size_t i = 0; i < n; i++) {
            double dx = x->data[i] - y->data[i];
            double dg = gx->data[i] - gy->data[i];
            diff_n2 += dx * dx;
            grad_diff_n2 += dg * dg;
        }
        double diff_norm = sqrt(diff_n2);
        double grad_diff_norm = sqrt(grad_diff_n2);

        if (diff_norm > 1e-10) {
            double ratio = grad_diff_norm / diff_norm;
            if (ratio > L_est) L_est = ratio;
        }
    }

    vector_free(x); vector_free(y);
    vector_free(gx); vector_free(gy);

    return (L_est > 0.0) ? L_est : 1.0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Nesterov Accelerated Gradient (FISTA-style)
 *
 * Accelerated scheme for smooth convex minimization:
 *   y_k = x_k + ((t_k-1)/t_{k+1}) (x_k - x_{k-1})   [extrapolation]
 *   x_{k+1} = y_k - alpha grad f(y_k)                [gradient step]
 *   t_{k+1} = (1 + sqrt(1+4t_k^2))/2                [sequence update]
 *
 * Theorem (Nesterov, 1983):
 * For L-smooth convex f with bounded level sets:
 *     f(x_k) - f* <= 2L ||x_0 - x*||^2 / (k+1)^2
 *
 * This O(1/k^2) rate is optimal among first-order methods
 * for smooth convex optimization (matching the lower bound).
 *
 * Key insight: the extrapolation step creates "momentum" that
 * dampens oscillations, allowing larger effective step sizes.
 *
 * Reference: Nesterov, "A method for solving the convex programming
 *            problem with convergence rate O(1/k^2)", Dokl. Akad. Nauk (1983)
 * ═══════════════════════════════════════════════════════════════════════ */

void nesterov_accelerated_gradient(const unconstr_problem_t *problem,
                                    double step_size,
                                    size_t max_iters,
                                    opt_result_t *result)
{
    assert(problem && result);
    size_t n = problem->n;

    vector_t *x = vector_alloc(n);     /* Current iterate */
    vector_t *y = vector_alloc(n);     /* Extrapolated point */
    vector_t *g = vector_alloc(n);     /* Gradient at y */
    vector_t *x_prev = vector_alloc(n);/* Previous iterate */

    vector_copy(&problem->x0, x);
    vector_copy(x, y);
    vector_copy(x, x_prev);

    double t = 1.0;  /* Nesterov sequence parameter */

    size_t iter;
    for (iter = 0; iter < max_iters; iter++) {
        /* Compute gradient at extrapolated point y_k */
        problem->f.grad(&problem->f, y, g);

        double gnorm = vector_norm_l2(g);
        if (gnorm < GRAD_TOL) {
            /* Converged at extrapolated point; update x to y */
            vector_copy(y, x);
            break;
        }

        /* Save current x as previous */
        vector_copy(x, x_prev);

        /* Gradient step from y: x_{k+1} = y_k - alpha * grad f(y_k) */
        for (size_t i = 0; i < n; i++) {
            x->data[i] = y->data[i] - step_size * g->data[i];
        }

        /* Update Nesterov sequence */
        double t_new = 0.5 * (1.0 + sqrt(1.0 + 4.0 * t * t));

        /* Extrapolation: y_{k+1} = x_{k+1} + ((t_k-1)/t_{k+1}) * (x_{k+1} - x_k) */
        double beta = (t - 1.0) / t_new;
        for (size_t i = 0; i < n; i++) {
            y->data[i] = x->data[i] + beta * (x->data[i] - x_prev->data[i]);
        }

        t = t_new;
    }

    result->status = (iter < max_iters) ? OPT_SUCCESS : OPT_MAXITER_REACHED;
    result->iterations = iter;
    result->f_value = problem->f.eval(&problem->f, x);
    result->grad_norm = vector_norm_l2(g);

    if (!result->x.data) {
        result->x.data = (double*)malloc(n * sizeof(double));
        result->x.n = n;
    }
    memcpy(result->x.data, x->data, n * sizeof(double));

    vector_free(x); vector_free(y);
    vector_free(g); vector_free(x_prev);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Heavy Ball Method (Polyak Momentum)
 *
 * x_{k+1} = x_k - alpha * grad f(x_k) + beta * (x_k - x_{k-1})
 *
 * For strongly convex quadratic problems, optimal parameters:
 *     alpha* = 4 / (sqrt(L) + sqrt(mu))^2
 *     beta*  = (sqrt(L) - sqrt(mu))^2 / (sqrt(L) + sqrt(mu))^2
 *
 * Theorem (Polyak, 1964):
 * Heavy ball with optimal parameters achieves linear convergence
 * with rate sqrt(kappa)-dependence instead of kappa-dependence:
 *     ||x_k - x*|| <= C * ((sqrt(kappa)-1)/(sqrt(kappa)+1))^k
 *
 * This matches the accelerated rate. The momentum term provides
 * "inertia" that accelerates convergence along low-curvature
 * directions while damping oscillations along high-curvature ones.
 *
 * Reference: Polyak, "Some methods of speeding up the convergence
 *            of iteration methods", USSR Comp. Math. Math. Phys. (1964)
 * ═══════════════════════════════════════════════════════════════════════ */

void heavy_ball_method(const unconstr_problem_t *problem,
                       double alpha, double beta,
                       size_t max_iters, opt_result_t *result)
{
    assert(problem && result);
    size_t n = problem->n;

    vector_t *x = vector_alloc(n);
    vector_t *x_prev = vector_alloc(n);
    vector_t *g = vector_alloc(n);

    vector_copy(&problem->x0, x);
    vector_copy(x, x_prev);

    size_t iter;
    for (iter = 0; iter < max_iters; iter++) {
        problem->f.grad(&problem->f, x, g);

        double gnorm = vector_norm_l2(g);
        if (gnorm < GRAD_TOL) break;

        /* Heavy ball update */
        vector_t *x_new = vector_alloc(n);
        for (size_t i = 0; i < n; i++) {
            x_new->data[i] = x->data[i] - alpha * g->data[i]
                           + beta * (x->data[i] - x_prev->data[i]);
        }

        /* Shift: x_prev <- x, x <- x_new */
        vector_copy(x, x_prev);
        vector_copy(x_new, x);
        vector_free(x_new);
    }

    result->status = (iter < max_iters) ? OPT_SUCCESS : OPT_MAXITER_REACHED;
    result->iterations = iter;
    result->f_value = problem->f.eval(&problem->f, x);
    result->grad_norm = vector_norm_l2(g);

    if (!result->x.data) {
        result->x.data = (double*)malloc(n * sizeof(double));
        result->x.n = n;
    }
    memcpy(result->x.data, x->data, n * sizeof(double));

    vector_free(x); vector_free(x_prev); vector_free(g);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Subgradient Method
 *
 * For nondifferentiable convex f, a subgradient g in df(x) satisfies:
 *     f(y) >= f(x) + g^T (y - x)  for all y
 *
 * Iteration: x_{k+1} = x_k - alpha_k * g_k, g_k in df(x_k)
 *
 * Key properties:
 * - NOT a descent method: f(x_{k+1}) may increase
 * - Must track best iterate: f_best = min_k f(x_k)
 * - Step sizes must be diminishing: alpha_k -> 0, sum alpha_k = infinity
 *
 * Convergence (fixed step alpha):
 *     f(x_best) - f* <= (||x_0 - x*||^2 + alpha^2 * G^2 * K) / (2 * alpha * K)
 *                    -> alpha * G^2 / 2 as K -> infinity
 *
 * Convergence (diminishing alpha_k = alpha_0 / sqrt(k)):
 *     f(x_best) - f* <= O(log(K) / sqrt(K))
 *
 * The O(1/sqrt(k)) rate is optimal for nonsmooth convex optimization.
 *
 * Reference: Shor (1985), "Minimization Methods for Non-Differentiable Functions"
 *            Nesterov (2004), "Introductory Lectures on Convex Optimization"
 * ═══════════════════════════════════════════════════════════════════════ */

void subgradient_method(const unconstr_problem_t *problem,
                        void (*subgrad)(const convex_function_t *f,
                                        const vector_t *x, vector_t *g),
                        double alpha0, size_t max_iters,
                        opt_result_t *result)
{
    assert(problem && subgrad && result);
    size_t n = problem->n;

    vector_t *x = vector_alloc(n);
    vector_t *x_best = vector_alloc(n);
    vector_t *g = vector_alloc(n);

    vector_copy(&problem->x0, x);
    vector_copy(x, x_best);

    double f_best = problem->f.eval(&problem->f, x);

    size_t iter;
    for (iter = 0; iter < max_iters; iter++) {
        subgrad(&problem->f, x, g);

        /* Diminishing step size: alpha_k = alpha_0 / sqrt(k+1) */
        double alpha = alpha0 / sqrt((double)(iter + 1));

        /* Subgradient step */
        double gnorm = vector_norm_l2(g);
        if (gnorm < 1e-12) break;

        for (size_t i = 0; i < n; i++) {
            x->data[i] -= (alpha / gnorm) * g->data[i];
        }

        /* Track best iterate (subgradient is not monotone) */
        double fx = problem->f.eval(&problem->f, x);
        if (fx < f_best) {
            f_best = fx;
            vector_copy(x, x_best);
        }
    }

    result->status = (iter < max_iters) ? OPT_SUCCESS : OPT_MAXITER_REACHED;
    result->iterations = iter;
    result->f_value = f_best;

    if (!result->x.data) {
        result->x.data = (double*)malloc(n * sizeof(double));
        result->x.n = n;
    }
    memcpy(result->x.data, x_best->data, n * sizeof(double));

    vector_free(x); vector_free(x_best); vector_free(g);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Coordinate Descent
 *
 * Minimizes f along one coordinate at a time via exact line search:
 *     x^{(i)}_{k+1} = argmin_t f(x_k + t * e_i)
 *
 * For separable functions f(x) = sum_i f_i(x_i), coordinate descent
 * converges in one pass. For general convex f, it converges with
 * the same asymptotic rate as gradient descent.
 *
 * The method is particularly effective when:
 * - Coordinate-wise minimization has closed form (LASSO, box constraints)
 * - n is very large (no full gradient needed)
 * - Each coordinate update is cheap (O(1) instead of O(n))
 *
 * Reference: Tseng (2001), "Convergence of a Block Coordinate Descent
 *            Method for Nondifferentiable Minimization"
 * ═══════════════════════════════════════════════════════════════════════ */

void coordinate_descent(const unconstr_problem_t *problem,
                        size_t n_cycles, opt_result_t *result)
{
    assert(problem && result);
    size_t n = problem->n;

    vector_t *x = vector_alloc(n);
    vector_t *g = vector_alloc(n);
    vector_t *e_i = vector_alloc(n);  /* Coordinate direction */

    vector_copy(&problem->x0, x);

    for (size_t cycle = 0; cycle < n_cycles; cycle++) {
        for (size_t i = 0; i < n; i++) {
            /* Set e_i to i-th coordinate direction */
            for (size_t j = 0; j < n; j++) {
                e_i->data[j] = (j == i) ? 1.0 : 0.0;
            }

            /* Compute gradient (for direction) */
            problem->f.grad(&problem->f, x, g);
            double gi = g->data[i];

            /* Simple 1D Newton/line search along coordinate i */
            /* For quadratic f(x) = (1/2) x^T Q x + q^T x:
             *   f(x + t*e_i) = (1/2) Q_ii t^2 + (Q x + q)_i t + const
             *   minimizer: t = -(Q x + q)_i / Q_ii = -g_i / Q_ii
             */

            /* Use finite difference to estimate second derivative */
            double eps = 1e-6;
            double f0 = problem->f.eval(&problem->f, x);

            x->data[i] += eps;
            double fp = problem->f.eval(&problem->f, x);
            x->data[i] -= 2.0 * eps;
            double fm = problem->f.eval(&problem->f, x);
            x->data[i] += eps;  /* Restore */

            double d2 = (fp - 2.0*f0 + fm) / (eps * eps);
            double step = 0.0;

            if (d2 > 1e-12) {
                /* Newton step for coordinate i */
                step = -gi / d2;
            } else {
                /* Gradient descent step along coordinate i */
                step = -0.1 * gi;
            }

            /* Safe step: don't move too far */
            if (fabs(step) > 10.0) step = (step > 0) ? 10.0 : -10.0;

            x->data[i] += step;
        }

        /* Check convergence after each full cycle */
        problem->f.grad(&problem->f, x, g);
        double gnorm = vector_norm_l2(g);
        if (gnorm < GRAD_TOL) {
            result->status = OPT_SUCCESS;
            result->iterations = (cycle + 1) * n;
            break;
        }
    }

    if (result->status != OPT_SUCCESS) {
        result->status = OPT_MAXITER_REACHED;
        result->iterations = n_cycles * n;
    }

    result->f_value = problem->f.eval(&problem->f, x);
    result->grad_norm = vector_norm_l2(g);

    if (!result->x.data) {
        result->x.data = (double*)malloc(n * sizeof(double));
        result->x.n = n;
    }
    memcpy(result->x.data, x->data, n * sizeof(double));

    vector_free(x); vector_free(g); vector_free(e_i);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Proximal Gradient Method (ISTA)
 *
 * Solves: minimize f(x) + g(x)
 * where f is smooth convex, g is convex with easy proximal operator.
 *
 * Iteration: x_{k+1} = prox_{alpha * g}(x_k - alpha * grad f(x_k))
 *
 * Key instances:
 * - LASSO: f(x) = (1/2)||Ax-b||^2, g(x) = lambda||x||_1
 *          prox = soft-thresholding: sign(x_i) * max(|x_i|-alpha*lambda, 0)
 * - Elastic net, group lasso, nuclear norm minimization
 *
 * Convergence: O(1/k) for ISTA, O(1/k^2) for FISTA (Nesterov+prox).
 *
 * Theorem (Combettes & Wajs, 2005; Beck & Teboulle, 2009):
 * For convex f+L, the proximal gradient method converges to a
 * minimizer with rate O(1/k). Each iteration requires one gradient
 * evaluation and one proximal operation.
 *
 * Reference: Parikh & Boyd (2014), "Proximal Algorithms",
 *            Foundations and Trends in Optimization.
 * ═══════════════════════════════════════════════════════════════════════ */

void proximal_gradient(const convex_function_t *f,
                       const convex_function_t *g,
                       const vector_t *x0,
                       double step_size,
                       size_t max_iters,
                       opt_result_t *result)
{
    assert(f && f->grad && g && g->prox && x0 && result);
    size_t n = x0->n;

    vector_t *x = vector_alloc(n);
    vector_t *grad = vector_alloc(n);
    vector_t *temp = vector_alloc(n);

    vector_copy(x0, x);

    size_t iter;
    for (iter = 0; iter < max_iters; iter++) {
        /* Gradient step: temp = x - alpha * grad f(x) */
        f->grad(f, x, grad);
        for (size_t i = 0; i < n; i++) {
            temp->data[i] = x->data[i] - step_size * grad->data[i];
        }

        /* Proximal step: x_new = prox_{alpha * g}(temp) */
        g->prox(g, step_size, temp, x);

        /* Fixed-point residual for convergence check */
        for (size_t i = 0; i < n; i++) {
            temp->data[i] = x->data[i] - temp->data[i];
        }
        double resid = vector_norm_l2(temp) / step_size;
        if (resid < GRAD_TOL) break;

        /* Restore temp for next iteration */
        for (size_t i = 0; i < n; i++) {
            temp->data[i] = x->data[i] - step_size * grad->data[i];
        }
    }

    result->status = (iter < max_iters) ? OPT_SUCCESS : OPT_MAXITER_REACHED;
    result->iterations = iter;
    result->f_value = f->eval(f, x) + g->eval(g, x);
    result->grad_norm = vector_norm_l2(grad);

    if (!result->x.data) {
        result->x.data = (double*)malloc(n * sizeof(double));
        result->x.n = n;
    }
    memcpy(result->x.data, x->data, n * sizeof(double));

    vector_free(x); vector_free(grad); vector_free(temp);
}
