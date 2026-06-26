/**
 * example_newton_method.c — Newton and BFGS: Quadratic Minimization
 *
 * Demonstrates Newton's method and BFGS on convex quadratic problems.
 * Compares convergence rates between first-order and second-order methods.
 *
 * f(x) = (1/2) x^T P x + q^T x  with P = [[3, 1], [1, 2]], q = [-2, -3]
 * Minimum at x* = P^{-1}(-q) = [1/5, 7/5] = [0.2, 1.4]
 * f* = -1.3
 */

#include "optimization_core.h"
#include "newton_methods.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* Quadratic: f(x) = (1/2) x^T P x + q^T x */
static double quad_eval(const convex_function_t *f, const vector_t *x)
{
    const double *P = f->param;
    const double *q = f->param + 4;
    double val = 0.0;
    val += q[0]*x->data[0] + q[1]*x->data[1];
    val += 0.5 * (P[0]*x->data[0]*x->data[0] + P[1]*x->data[0]*x->data[1]
                + P[2]*x->data[1]*x->data[0] + P[3]*x->data[1]*x->data[1]);
    return val;
}

static void quad_grad(const convex_function_t *f, const vector_t *x, vector_t *g)
{
    const double *P = f->param;
    const double *q = f->param + 4;
    g->data[0] = P[0]*x->data[0] + P[1]*x->data[1] + q[0];
    g->data[1] = P[2]*x->data[0] + P[3]*x->data[1] + q[1];
}

static void quad_hess(const convex_function_t *f, const vector_t *x, matrix_t *H)
{
    (void)x;
    const double *P = f->param;
    MAT_ELT(H, 0, 0) = P[0]; MAT_ELT(H, 0, 1) = P[1];
    MAT_ELT(H, 1, 0) = P[2]; MAT_ELT(H, 1, 1) = P[3];
}

int main(void)
{
    printf("=== Newton and BFGS: Quadratic Minimization ===\n");

    /* P = [[3, 1], [1, 2]], q = [-2, -3] */
    double param[] = {3.0, 1.0, 1.0, 2.0, -2.0, -3.0};

    convex_function_t quad_fn;
    memset(&quad_fn, 0, sizeof(quad_fn));
    quad_fn.n = 2;
    quad_fn.param = param;
    quad_fn.param_len = 6;
    quad_fn.eval = quad_eval;
    quad_fn.grad = quad_grad;
    quad_fn.hess = quad_hess;

    unconstr_problem_t prob;
    memset(&prob, 0, sizeof(prob));
    prob.n = 2;
    prob.f = quad_fn;
    prob.has_analytic = 1;
    prob.opt_value = -1.3;

    /* x0 = [10, 10] — far from optimum */
    double x0d[] = {10.0, 10.0};
    prob.x0.data = x0d;
    prob.x0.n = 2;

    /* ─── Newton's Method ─── */
    printf("\n--- Newton's Method ---\n");
    newton_config_t ncfg = newton_config_default();
    ncfg.max_iters = 100;
    ncfg.grad_tol = 1e-10;

    opt_result_t n_result;
    memset(&n_result, 0, sizeof(n_result));
    newton_method(&prob, &ncfg, &n_result);

    printf("  Status: %d, Iters: %zu\n", n_result.status, n_result.iterations);
    printf("  Solution: (%.10f, %.10f)\n", n_result.x.data[0], n_result.x.data[1]);
    printf("  Objective: %.10f (opt = %.1f)\n", n_result.f_value, prob.opt_value);
    printf("  Gradient norm: %.2e\n", n_result.grad_norm);
    free(n_result.x.data);

    /* ─── BFGS ─── */
    printf("\n--- BFGS Method ---\n");
    prob.x0.data[0] = 10.0;
    prob.x0.data[1] = 10.0;

    opt_result_t b_result;
    memset(&b_result, 0, sizeof(b_result));
    bfgs_method(&prob, 1000, &b_result);

    printf("  Status: %d, Iters: %zu\n", b_result.status, b_result.iterations);
    printf("  Solution: (%.10f, %.10f)\n", b_result.x.data[0], b_result.x.data[1]);
    printf("  Objective: %.10f\n", b_result.f_value);
    printf("  Gradient norm: %.2e\n", b_result.grad_norm);
    free(b_result.x.data);

    /* ─── Damped Newton ─── */
    printf("\n--- Damped Newton ---\n");
    prob.x0.data[0] = 10.0;
    prob.x0.data[1] = 10.0;

    opt_result_t dn_result;
    memset(&dn_result, 0, sizeof(dn_result));
    damped_newton_method(&prob, &ncfg, &dn_result);

    printf("  Status: %d, Iters: %zu\n", dn_result.status, dn_result.iterations);
    printf("  Solution: (%.10f, %.10f)\n", dn_result.x.data[0], dn_result.x.data[1]);
    printf("  Objective: %.10f\n", dn_result.f_value);
    free(dn_result.x.data);

    /* ─── L-BFGS ─── */
    printf("\n--- L-BFGS (m=5) ---\n");
    prob.x0.data[0] = 10.0;
    prob.x0.data[1] = 10.0;

    opt_result_t lb_result;
    memset(&lb_result, 0, sizeof(lb_result));
    lbfgs_method(&prob, 5, 1000, &lb_result);

    printf("  Status: %d, Iters: %zu\n", lb_result.status, lb_result.iterations);
    printf("  Solution: (%.10f, %.10f)\n", lb_result.x.data[0], lb_result.x.data[1]);
    printf("  Objective: %.10f\n", lb_result.f_value);
    free(lb_result.x.data);

    printf("\n=== Newton Method Example Complete ===\n");
    return 0;
}
