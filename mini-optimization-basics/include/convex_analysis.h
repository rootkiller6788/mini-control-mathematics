/**
 * convex_analysis.h — Convex Analysis and Duality Theory
 *
 * Foundational tools for convex optimization:
 * convexity verification, conjugate functions, Fenchel duality,
 * monotone operators, separating hyperplanes.
 *
 * L1: Convex sets, convex functions, epigraphs, perspective
 * L2: Operations preserving convexity, conjugate functions
 * L4: Separating hyperplane theorem, Fenchel duality
 * L8: Subdifferential calculus, Moreau decomposition
 *
 * Reference: Rockafellar, "Convex Analysis" (1970)
 *            Boyd & Vandenberghe, "Convex Optimization" (2004), Ch.2,3
 */

#ifndef CONVEX_ANALYSIS_H
#define CONVEX_ANALYSIS_H

#include "optimization_core.h"

int convex_set_verify(const convex_set_t *C, size_t n, size_t nsamples);

int convex_function_verify(const convex_function_t *f, size_t nsamples);

double convex_conjugate_eval(const convex_function_t *f, const vector_t *y);

int fenchel_conjugate_compute(const convex_function_t *f,
                               const vector_t *y,
                               const vector_t *x0,
                               size_t max_iters,
                               double *conj_val,
                               vector_t *x_opt);

convex_function_t convex_sum(const convex_function_t *funcs,
                              const double *weights, size_t k);

convex_function_t convex_max(const convex_function_t *funcs, size_t k);

convex_function_t convex_affine_compose(const convex_function_t *f,
                                          const matrix_t *A,
                                          const vector_t *b);

int separating_hyperplane(const convex_set_t *C1, const convex_set_t *C2,
                          vector_t *a, double *b);

int supporting_hyperplane(const convex_set_t *C, const vector_t *x,
                          vector_t *a, double *b);

void project_l2_ball(const vector_t *v, const vector_t *center,
                     double radius, vector_t *proj);
void project_l1_ball(const vector_t *v, double radius, vector_t *proj);
void project_linf_ball(const vector_t *v, double radius, vector_t *proj);
void project_simplex(const vector_t *v, vector_t *proj);
void project_box(const vector_t *v, const vector_t *lower,
                 const vector_t *upper, vector_t *proj);
void project_halfspace(const vector_t *v, const vector_t *a,
                       double b, vector_t *proj);

void moreau_decomposition(const convex_function_t *f,
                          const vector_t *v,
                          vector_t *prox_val,
                          vector_t *comp_val);

void fenchel_dual_solve(const convex_function_t *f,
                        const convex_function_t *g,
                        const matrix_t *A,
                        const vector_t *y0,
                        size_t max_iters,
                        opt_result_t *result,
                        vector_t *y_opt);

void subgradient_compute(const convex_function_t *f,
                         const vector_t *x, vector_t *g);

int alternating_projections(int m,
                             void (**proj_funcs)(const vector_t*, vector_t*),
                             const vector_t *x0, size_t max_iters,
                             vector_t *x_out);

double perspective_eval(const convex_function_t *f,
                        const vector_t *x, double t);

double infimal_convolution(const convex_function_t *f,
                           const convex_function_t *g,
                           const vector_t *x);

void moreau_envelope(const convex_function_t *f, double lambda,
                     const vector_t *x, double *m_env, vector_t *m_grad);

#endif /* CONVEX_ANALYSIS_H */
