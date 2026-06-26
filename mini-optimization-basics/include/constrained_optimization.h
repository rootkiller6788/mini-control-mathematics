/**
 * constrained_optimization.h — Constrained Optimization Methods
 *
 * Covers Lagrange multipliers, KKT conditions, duality theory,
 * linear programming (simplex), and quadratic programming.
 *
 * L4: KKT Theorem, Weak/Strong Duality, Slater's condition
 * L5: Simplex method, Lagrange multiplier method, interior point
 * L6: Linear programs, quadratic programs, resource allocation
 *
 * Reference: Boyd & Vandenberghe, "Convex Optimization" (2004), Ch.5
 *            Bertsekas, "Constrained Optimization & Lagrange Multiplier Methods" (1982)
 *            Dantzig, "Linear Programming and Extensions" (1963)
 */

#ifndef CONSTRAINED_OPTIMIZATION_H
#define CONSTRAINED_OPTIMIZATION_H

#include "optimization_core.h"

double lagrangian_eval(const general_nlp_t *nlp, const vector_t *x,
                       double *result);

void lagrangian_gradient(const general_nlp_t *nlp, const vector_t *x,
                          vector_t *out);

int kkt_check(const general_nlp_t *nlp, const vector_t *x,
              const double *lambdas, const double *nus);

int slater_condition_check(const general_nlp_t *nlp);

int dual_function_eval(const general_nlp_t *nlp,
                       const double *lambdas, const double *nus,
                       double *dual_val);

void duality_gap(const general_nlp_t *nlp, const vector_t *x,
                 const double *lambdas, const double *nus,
                 double *dual_gap);

void lagrange_multiplier_solve(const eq_constrained_problem_t *problem,
                                opt_result_t *result);

void quadratic_penalty_method(const eq_constrained_problem_t *problem,
                               double rho_init, double rho_factor,
                               size_t max_outer, opt_result_t *result);

void log_barrier_method(const ineq_constrained_problem_t *problem,
                         double t_init, double mu,
                         size_t max_outer, opt_result_t *result);

typedef struct {
    size_t    max_iters;
    int       phase_one_only;
    int       verbose;
    int       pivot_rule;
} simplex_config_t;

simplex_config_t simplex_config_default(void);

void simplex_solve(linear_program_t *lp,
                   const simplex_config_t *config,
                   opt_result_t *result);

void lp_dual_solve(const linear_program_t *lp, vector_t *y_dual,
                    double *dual_val);

void qp_active_set_solve(quadratic_program_t *qp,
                          size_t max_iters, opt_result_t *result);

void projected_gradient(const unconstr_problem_t *problem,
                        const convex_set_t *C,
                        double step_size,
                        size_t max_iters,
                        opt_result_t *result);

int kkt_system_solve(const matrix_t *P, const vector_t *q,
                     const matrix_t *A, const vector_t *b,
                     vector_t *x, vector_t *nu);

int transportation_lp_solve(size_t n_suppliers, size_t n_demands,
                             const matrix_t *costs,
                             const vector_t *supply,
                             const vector_t *demand,
                             matrix_t *flow, double *opt_cost);

#endif /* CONSTRAINED_OPTIMIZATION_H */
