/**
 * newton_methods.h — Newton-Type Second-Order Methods
 *
 * Implements Newton's method and quasi-Newton methods (BFGS, DFP)
 * for unconstrained optimization. Second-order methods exploit
 * Hessian curvature information for faster convergence.
 *
 * L2: Newton's method, quasi-Newton methods, secant condition
 * L5: Damped Newton, BFGS, DFP, SR1 methods
 * L8: L-BFGS (limited memory), trust-region Newton
 *
 * Reference: Nocedal & Wright, "Numerical Optimization" (2006), Ch.3,6,7
 */

#ifndef NEWTON_METHODS_H
#define NEWTON_METHODS_H

#include "optimization_core.h"

typedef struct {
    size_t  max_iters;
    double  grad_tol;
    double  alpha_min;
    double  beta;
    double  c1;
    int     use_cholesky;
    int     verbose;
} newton_config_t;

newton_config_t newton_config_default(void);

void newton_method(const unconstr_problem_t *problem,
                   const newton_config_t *config,
                   opt_result_t *result);

double newton_decrement(const convex_function_t *f, const vector_t *x,
                        const vector_t *grad, matrix_t *H, vector_t *work);

void damped_newton_method(const unconstr_problem_t *problem,
                          const newton_config_t *config,
                          opt_result_t *result);

void bfgs_method(const unconstr_problem_t *problem,
                 size_t max_iters, opt_result_t *result);

void dfp_method(const unconstr_problem_t *problem,
                size_t max_iters, opt_result_t *result);

void sr1_method(const unconstr_problem_t *problem,
                size_t max_iters, opt_result_t *result);

void lbfgs_method(const unconstr_problem_t *problem,
                  size_t m, size_t max_iters,
                  opt_result_t *result);

void trust_region_newton(const unconstr_problem_t *problem,
                         const newton_config_t *config,
                         opt_result_t *result);

void gauss_newton_method(size_t n_vars, size_t n_resid,
                         void (*r)(const vector_t *x, vector_t *r_out),
                         void (*J)(const vector_t *x, matrix_t *J_out),
                         const vector_t *x0, size_t max_iters,
                         opt_result_t *result);

void levenberg_marquardt(size_t n_vars, size_t n_resid,
                         void (*r)(const vector_t *x, vector_t *r_out),
                         void (*J)(const vector_t *x, matrix_t *J_out),
                         const vector_t *x0, double lambda0,
                         size_t max_iters, opt_result_t *result);

#endif /* NEWTON_METHODS_H */
