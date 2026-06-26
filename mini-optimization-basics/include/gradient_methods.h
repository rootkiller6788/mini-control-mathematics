/**
 * gradient_methods.h — First-Order Optimization Methods
 *
 * Implements gradient-based optimization algorithms for
 * unconstrained and constrained problems.
 *
 * L2 Concepts: Gradient descent, steepest descent, conjugate gradient
 * L5 Methods: Line search, momentum, Nesterov acceleration, prox-gradient
 * L3 Quantities: Learning rates, Lipschitz constants, convergence rates
 *
 * Reference: Nocedal & Wright, "Numerical Optimization" (2006), Ch.2-3,5
 *            Boyd & Vandenberghe, "Convex Optimization" (2004), Ch.9
 */

#ifndef GRADIENT_METHODS_H
#define GRADIENT_METHODS_H

#include "optimization_core.h"

typedef struct {
    double       step_size;
    linesearch_type_t ls_type;
    backtracking_params_t ls_params;
    double       grad_tol;
    size_t       max_iters;
    int          use_nesterov;
    double       momentum_beta;
    int          verbose;
} gradient_descent_config_t;

gradient_descent_config_t gd_config_default(void);

void gradient_descent(const unconstr_problem_t *problem,
                      const gradient_descent_config_t *config,
                      opt_result_t *result);

double backtracking_linesearch(const convex_function_t *f,
                                const vector_t *x,
                                const vector_t *d,
                                const vector_t *grad_at_x,
                                const backtracking_params_t *params);

double steepest_descent_exact_step(const convex_function_t *f,
                                    vector_t *x, vector_t *g);

typedef struct {
    size_t  max_iters;
    double  grad_tol;
    int     variant;
    int     verbose;
} cg_config_t;

cg_config_t cg_config_default(void);

void conjugate_gradient_linear(const matrix_t *Q, const vector_t *b,
                                vector_t *x, const cg_config_t *config);

void nonlinear_cg(const unconstr_problem_t *problem,
                  const cg_config_t *config, opt_result_t *result);

double estimate_lipschitz(const convex_function_t *f, size_t n, size_t nsamples);

void proximal_gradient(const convex_function_t *f,
                       const convex_function_t *g,
                       const vector_t *x0,
                       double step_size,
                       size_t max_iters,
                       opt_result_t *result);

void nesterov_accelerated_gradient(const unconstr_problem_t *problem,
                                    double step_size,
                                    size_t max_iters,
                                    opt_result_t *result);

void heavy_ball_method(const unconstr_problem_t *problem,
                       double alpha, double beta,
                       size_t max_iters, opt_result_t *result);

void subgradient_method(const unconstr_problem_t *problem,
                        void (*subgrad)(const convex_function_t *f,
                                        const vector_t *x, vector_t *g),
                        double alpha0, size_t max_iters,
                        opt_result_t *result);

void coordinate_descent(const unconstr_problem_t *problem,
                        size_t n_cycles, opt_result_t *result);

#endif /* GRADIENT_METHODS_H */
