/**
 * optimization.h — Numerical Optimization for Control
 *
 * Covers: L5 Methods (optimization algorithms), L8 Advanced Methods
 *
 * Optimization is central to control: LQR minimizes quadratic cost, MPC solves
 * constrained optimization online, system identification minimizes prediction
 * error, and robust control solves min-max problems.
 *
 * Reference: Nocedal & Wright, Numerical Optimization, 2nd ed. (2006)
 *            Boyd & Vandenberghe, Convex Optimization (2004)
 *            Luenberger & Ye, Linear and Nonlinear Programming (2016)
 */

#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include "numerical_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Optimization Problem Definitions
 * ========================================================================== */

/**
 * ScalarFunction — 1D function f: ℝ → ℝ
 *
 * L1 Definition: Univariate real-valued function used in root-finding
 * and line search. The function signature includes a params pointer for
 * passing problem-specific data through the optimization algorithms.
 */
typedef double (*ScalarFunction)(double x, void *params);

/**
 * MultivariateFunction — f: ℝⁿ → ℝ
 *
 * L1 Definition: Cost function in optimal control. For LQR:
 * J = ∫ (xᵀQx + uᵀRu) dt. For parameter estimation: J = Σ (yᵢ - ŷᵢ)².
 */
typedef double (*MultivariateFunction)(const double *x, size_t n, void *params);

/**
 * MultivariateGradient — ∇f: ℝⁿ → ℝⁿ
 *
 * L1 Definition: Gradient vector of partial derivatives.
 * Required for gradient-based optimization methods.
 */
typedef void (*MultivariateGradient)(const double *x, size_t n, void *params,
                                      double *grad);

/**
 * OptResult — optimization result
 *
 * L1 Definition: Stores the optimal point, function value, and convergence info.
 */
typedef struct {
    double *x_opt;        /**< Optimal point */
    double  f_opt;        /**< Function value at optimum */
    double  gradient_norm;/**< ||∇f(x_opt)|| */
    size_t  iterations;   /**< Number of iterations */
    int     converged;    /**< 1 if convergence criterion met */
} OptResult;

/* ==========================================================================
 * L5: Root-Finding Methods (1D)
 * ========================================================================== */

/**
 * bisection_method — find root of f(x) = 0 on [a, b]
 *
 * L5 Method: Binary search. Guaranteed convergence if f(a)·f(b) < 0.
 * Convergence: linear, rate 1/2. Robust but slow.
 * Complexity: O(log((b-a)/tol)) function evaluations.
 */
double bisection_method(ScalarFunction f, void *params,
                        double a, double b, double tol, size_t max_iter,
                        int *converged);

/**
 * newton_raphson — Newton-Raphson root finding
 *
 * L5 Method: x_{k+1} = x_k - f(x_k)/f'(x_k).
 * Convergence: quadratic near root (|e_{k+1}| ≤ C·|e_k|²).
 * Requires derivative. May diverge far from root.
 */
double newton_raphson(ScalarFunction f, ScalarFunction fprime, void *params,
                      double x0, double tol, size_t max_iter, int *converged);

/**
 * secant_method — Secant method (derivative-free Newton)
 *
 * L5 Method: x_{k+1} = x_k - f(x_k)·(x_k - x_{k-1})/(f(x_k) - f(x_{k-1})).
 * Convergence: superlinear, order φ ≈ 1.618 (golden ratio).
 * Two initial points required.
 */
double secant_method(ScalarFunction f, void *params,
                     double x0, double x1, double tol, size_t max_iter,
                     int *converged);

/* ==========================================================================
 * L5: First-Order Optimization Methods
 * ========================================================================== */

/**
 * gradient_descent — standard gradient descent with fixed step size
 *
 * L5 Method: x_{k+1} = x_k - α ∇f(x_k).
 * Convergence: linear for strongly convex functions.
 * Used in iterative learning control and online parameter estimation.
 */
OptResult* gradient_descent(MultivariateFunction f, MultivariateGradient grad,
                            void *params, size_t n, const double *x0,
                            double step_size, double tol, size_t max_iter);

/**
 * gradient_descent_backtracking — gradient descent with Armijo backtracking
 *
 * L5 Method: Adaptive step size via Armijo condition:
 * f(x_k + α·d_k) ≤ f(x_k) + c·α·∇f(x_k)ᵀd_k  with d_k = -∇f(x_k).
 * Ensures sufficient decrease. More robust than fixed step.
 */
OptResult* gradient_descent_backtracking(MultivariateFunction f,
                                         MultivariateGradient grad,
                                         void *params, size_t n,
                                         const double *x0, double init_step,
                                         double tol, size_t max_iter);

/**
 * gradient_descent_momentum — gradient descent with momentum (Polyak heavy ball)
 *
 * L8 Method: v_{k+1} = β·v_k + α·∇f(x_k); x_{k+1} = x_k - v_{k+1}.
 * Accelerates convergence, especially for ill-conditioned problems.
 * Used in Nesterov accelerated gradient methods for optimal control.
 */
OptResult* gradient_descent_momentum(MultivariateFunction f,
                                     MultivariateGradient grad,
                                     void *params, size_t n, const double *x0,
                                     double step_size, double momentum,
                                     double tol, size_t max_iter);

/* ==========================================================================
 * L5/L8: Conjugate Gradient Methods
 * ========================================================================== */

/**
 * conjugate_gradient_linear — solve Ax = b for SPD A using CG
 *
 * L5 Method: Iterative Krylov subspace method, optimal over Krylov subspaces.
 * Convergence: O(√κ) iterations for well-conditioned problems.
 * Only requires matrix-vector products (no explicit A storage needed).
 * Essential for large-scale MPC (interior point methods use CG).
 */
Vector* conjugate_gradient_linear(const Matrix *A, const Vector *b,
                                  const Vector *x0, double tol, size_t max_iter);

/**
 * conjugate_gradient_nonlinear — Fletcher-Reeves nonlinear CG
 *
 * L8 Method: Generalization of linear CG to nonlinear optimization.
 * β_k^{FR} = ||∇f_{k+1}||² / ||∇f_k||².
 * Low storage: only requires gradient evaluations, no Hessian.
 */
OptResult* conjugate_gradient_nonlinear(MultivariateFunction f,
                                        MultivariateGradient grad,
                                        void *params, size_t n,
                                        const double *x0, double tol,
                                        size_t max_iter);

/* ==========================================================================
 * L8: Second-Order Methods
 * ========================================================================== */

/**
 * newton_method_opt — Newton's method for unconstrained optimization
 *
 * L8 Method: x_{k+1} = x_k - [∇²f(x_k)]⁻¹ ∇f(x_k).
 * Quadratic convergence near optimum. Requires Hessian.
 *
 * For optimal control: Newton methods solve the Riccati equation iteratively
 * in Differential Dynamic Programming (DDP).
 */
OptResult* newton_method_opt(MultivariateFunction f, MultivariateGradient grad,
                             void *params, size_t n, const double *x0,
                             double tol, size_t max_iter);

/**
 * bfgs_update — BFGS update for approximate Hessian
 *
 * L8 Method: B_{k+1} = B_k + (y_k y_kᵀ)/(y_kᵀ s_k) - (B_k s_k s_kᵀ B_k)/(s_kᵀ B_k s_k)
 * where s_k = x_{k+1} - x_k, y_k = ∇f_{k+1} - ∇f_k.
 * Maintains SPD approximation with O(n²) storage.
 */
void bfgs_update(Matrix *B, const double *s, const double *y, size_t n);

/**
 * line_search_wolfe — Wolfe condition line search
 *
 * L8 Method: Finds step size α satisfying both Armijo (sufficient decrease)
 * and curvature condition: ∇f(x_k + α·d)ᵀd ≥ c₂·∇f(x_k)ᵀd.
 * Ensures convergence of BFGS and nonlinear CG.
 */
double line_search_wolfe(MultivariateFunction f, MultivariateGradient grad,
                         void *params, size_t n, const double *x,
                         const double *direction, double max_step);

/**
 * opt_result_free — deallocate optimization result
 */
void opt_result_free(OptResult *r);

#ifdef __cplusplus
}
#endif

#endif /* OPTIMIZATION_H */
