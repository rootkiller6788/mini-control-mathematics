/**
 * example_gradient_descent.c — Gradient Descent: Rosenbrock Function
 *
 * Minimizes the 2D Rosenbrock function ("banana function"):
 *     f(x,y) = (1-x)^2 + 100*(y - x^2)^2
 *
 * Global minimum at (1, 1) with f* = 0.
 *
 * Demonstrates gradient descent with backtracking line search
 * and conjugate gradient methods.
 */

#include "optimization_core.h"
#include "gradient_methods.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* Rosenbrock function: f(x,y) = (1-x)^2 + 100*(y-x^2)^2 */
static double rosen_eval(const convex_function_t *f, const vector_t *x)
{
    (void)f;
    double a = x->data[1] - x->data[0]*x->data[0];
    double b = 1.0 - x->data[0];
    return 100.0 * a*a + b*b;
}

static void rosen_grad(const convex_function_t *f, const vector_t *x, vector_t *g)
{
    (void)f;
    double a = x->data[1] - x->data[0]*x->data[0];
    g->data[0] = -400.0 * x->data[0] * a - 2.0 * (1.0 - x->data[0]);
    g->data[1] = 200.0 * a;
}

static void rosen_hess(const convex_function_t *f, const vector_t *x, matrix_t *H)
{
    (void)f;
    /* Hessian of Rosenbrock */
    double a = x->data[1] - x->data[0]*x->data[0];
    MAT_ELT(H, 0, 0) = -400.0*a + 800.0*x->data[0]*x->data[0] + 2.0;
    MAT_ELT(H, 0, 1) = -400.0 * x->data[0];
    MAT_ELT(H, 1, 0) = -400.0 * x->data[0];
    MAT_ELT(H, 1, 1) = 200.0;
}

int main(void)
{
    printf("=== Gradient Descent: Rosenbrock Function ===\n\n");

    convex_function_t rosen;
    memset(&rosen, 0, sizeof(rosen));
    rosen.n = 2;
    rosen.eval = rosen_eval;
    rosen.grad = rosen_grad;
    rosen.hess = rosen_hess;

    unconstr_problem_t prob;
    memset(&prob, 0, sizeof(prob));
    prob.n = 2;
    prob.f = rosen;
    prob.has_analytic = 1;
    prob.opt_value = 0.0;

    /* Start from (-1.2, 1.0) — classic starting point */
    double x0d[] = {-1.2, 1.0};
    prob.x0.data = x0d;
    prob.x0.n = 2;

    /* ─── Gradient Descent ─── */
    printf("--- Gradient Descent (backtracking) ---\n");
    gradient_descent_config_t gd_cfg = gd_config_default();
    gd_cfg.ls_type = LS_BACKTRACKING;
    gd_cfg.grad_tol = 1e-6;
    gd_cfg.max_iters = 10000;
    gd_cfg.verbose = 0;

    opt_result_t gd_result;
    memset(&gd_result, 0, sizeof(gd_result));
    gradient_descent(&prob, &gd_cfg, &gd_result);

    printf("  Status: %d, Iters: %zu\n", gd_result.status, gd_result.iterations);
    printf("  Solution: (%.6f, %.6f)\n", gd_result.x.data[0], gd_result.x.data[1]);
    printf("  Objective: %.10f\n", gd_result.f_value);
    printf("  Gradient norm: %.8f\n", gd_result.grad_norm);
    free(gd_result.x.data);

    /* ─── Nonlinear CG ─── */
    printf("\n--- Nonlinear Conjugate Gradient (Fletcher-Reeves) ---\n");
    prob.x0.data[0] = -1.2;
    prob.x0.data[1] = 1.0;

    cg_config_t cg_cfg = cg_config_default();
    cg_cfg.max_iters = 5000;
    cg_cfg.grad_tol = 1e-6;

    opt_result_t cg_result;
    memset(&cg_result, 0, sizeof(cg_result));
    nonlinear_cg(&prob, &cg_cfg, &cg_result);

    printf("  Status: %d, Iters: %zu\n", cg_result.status, cg_result.iterations);
    printf("  Solution: (%.6f, %.6f)\n", cg_result.x.data[0], cg_result.x.data[1]);
    printf("  Objective: %.10f\n", cg_result.f_value);
    printf("  Gradient norm: %.8f\n", cg_result.grad_norm);
    free(cg_result.x.data);

    /* ─── Nesterov Accelerated Gradient ─── */
    printf("\n--- Nesterov Accelerated Gradient ---\n");
    prob.x0.data[0] = -1.2;
    prob.x0.data[1] = 1.0;

    opt_result_t nag_result;
    memset(&nag_result, 0, sizeof(nag_result));
    nesterov_accelerated_gradient(&prob, 0.001, 10000, &nag_result);

    printf("  Status: %d, Iters: %zu\n", nag_result.status, nag_result.iterations);
    printf("  Solution: (%.6f, %.6f)\n", nag_result.x.data[0], nag_result.x.data[1]);
    printf("  Objective: %.10f\n", nag_result.f_value);
    printf("  Gradient norm: %.8f\n", nag_result.grad_norm);
    free(nag_result.x.data);

    printf("\n=== Gradient Descent Example Complete ===\n");
    return 0;
}
