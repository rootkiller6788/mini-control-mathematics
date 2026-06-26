/**
 * nonlinear_methods.h — Nonlinear Numerical Methods for Control
 *
 * Covers: L8 Advanced Methods, L9 Frontiers
 *
 * Nonlinear control systems ẋ = f(x, u) require specialized numerical methods
 * beyond linear algebra. This module covers Newton-Kantorovich for nonlinear
 * systems, continuation/homotopy, shooting methods for BVPs, and
 * computational methods for nonlinear model predictive control.
 *
 * Reference: Deuflhard, Newton Methods for Nonlinear Problems (2004)
 *            Allgower & Georg, Numerical Continuation Methods (1990)
 *            Ascher, Mattheij, Russell, Numerical Solution of BVPs for ODEs (1995)
 */

#ifndef NONLINEAR_METHODS_H
#define NONLINEAR_METHODS_H

#include "numerical_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Nonlinear System Types
 * ========================================================================== */

/**
 * NLEquation — nonlinear equation F(x) = 0, F: ℝⁿ → ℝⁿ
 *
 * L1 Definition: System of n nonlinear equations in n unknowns.
 * Examples: trim conditions for aircraft (forces and moments balance),
 * equilibrium points of nonlinear dynamics.
 */
typedef struct {
    /** Evaluate F at point x, store result in Fx */
    void (*eval)(const double *x, double *F, void *params, size_t n);
    /** Evaluate Jacobian J_{ij} = ∂F_i/∂x_j at point x */
    void (*jacobian)(const double *x, Matrix *J, void *params, size_t n);
    void  *params;
    size_t n;
} NLEquation;

/**
 * NLSolution — solution to nonlinear system
 */
typedef struct {
    double *x;           /**< Solution point */
    double  residual;    /**< ||F(x)|| at solution */
    size_t  iterations;  /**< Newton iterations */
    int     converged;
} NLSolution;

/**
 * BVPProblem — two-point boundary value problem
 *
 * L8 Definition: ẋ = f(t, x) with boundary conditions g(x(t₀), x(t_f)) = 0.
 * Arises in optimal control (Pontryagin's principle yields a BVP).
 */
typedef struct {
    void (*dynamics)(double t, const double *x, double *dx, void *params, size_t n);
    void (*boundary)(const double *x0, const double *xf, double *bc, void *params, size_t n);
    void  *params;
    size_t n;       /**< State dimension */
    double t0;      /**< Initial time */
    double tf;      /**< Final time */
} BVPProblem;

/* ==========================================================================
 * L8: Newton-Kantorovich and Quasi-Newton Methods
 * ========================================================================== */

/**
 * newton_kantorovich — Newton's method for nonlinear systems
 *
 * L8 Method: x_{k+1} = x_k - J(x_k)⁻¹ F(x_k).
 * Kantorovich theorem: guarantees convergence if initial guess is close enough
 * and Jacobian is Lipschitz continuous. Quadratic convergence locally.
 *
 * L7 Application: Solving steady-state/equilibrium conditions for nonlinear
 * power systems, chemical reactors (Arrhenius kinetics), and aircraft trim.
 */
NLSolution* newton_kantorovich(const NLEquation *eq, const double *x0,
                                double tol, size_t max_iter);

/**
 * newton_kantorovich_global — globally convergent Newton with line search
 *
 * L8 Method: Damped Newton: x_{k+1} = x_k - λ_k J(x_k)⁻¹ F(x_k)
 * where λ_k is chosen via Armijo backtracking on ||F(x)||.
 * Converges from poor initial guesses.
 */
NLSolution* newton_kantorovich_global(const NLEquation *eq, const double *x0,
                                       double tol, size_t max_iter);

/**
 * broyden_method — Broyden's quasi-Newton method (no Jacobian required)
 *
 * L8 Method: Maintains approximate Jacobian B_k via rank-1 updates:
 * B_{k+1} = B_k + (y_k - B_k s_k) s_kᵀ / (s_kᵀ s_k)
 * Superlinear convergence. Only O(n²) per iteration vs O(n³) for Newton.
 *
 * L7 Application: Real-time nonlinear MPC where Jacobian computation
 * is expensive (e.g., CFD-based models).
 */
NLSolution* broyden_method(const NLEquation *eq, const double *x0,
                            double tol, size_t max_iter);

/* ==========================================================================
 * L8: Continuation / Homotopy Methods
 * ========================================================================== */

/**
 * continuation_natural — natural parameter continuation
 *
 * L8 Method: Solve F(x, λ) = 0 for a sequence of λ values.
 * Uses predictor-corrector: tangent predictor, Newton corrector.
 * Traces solution branches through parameter space including turning points.
 *
 * L9 Applications: Bifurcation analysis in nonlinear dynamics,
 * optimal control continuation with respect to cost function weight.
 */
NLSolution* continuation_natural(NLEquation *eq, const double *x0,
                                  double lambda0, double lambda_end,
                                  double delta_lambda, double tol);

/**
 * continuation_pseudo_arclength — pseudo-arclength continuation
 *
 * L8 Method: Parameterization by arclength s along solution curve.
 * Handles turning points (folds) where natural continuation fails.
 * The standard method for bifurcation analysis (AUTO, MATCONT).
 *
 * L9: Used in computational bifurcation theory to trace limit cycles
 * and study Hopf bifurcations in nonlinear control systems.
 */
NLSolution* continuation_pseudo_arclength(NLEquation *eq, const double *x0,
                                           double lambda0, double lambda_end,
                                           double ds, double tol);

/* ==========================================================================
 * L8: Shooting Methods for BVPs
 * ========================================================================== */

/**
 * shooting_simple — simple (single) shooting for BVPs
 *
 * L8 Method: Parameterize unknown initial condition, integrate forward,
 * match terminal condition via Newton iteration. Works for mildly nonlinear BVPs.
 *
 * L7 Application: Solving optimal control BVPs from Pontryagin's minimum
 * principle. Example: minimum-time spacecraft rendezvous (Apollo guidance).
 */
NLSolution* shooting_simple(BVPProblem *bvp, const double *x0_guess,
                             size_t n_steps, double tol, size_t max_iter);

/**
 * shooting_multiple — multiple shooting for BVPs
 *
 * L8 Method: Divide [t₀, t_f] into subintervals, shoot on each, enforce
 * continuity. More stable than simple shooting for sensitive problems.
 *
 * L7 Application: Re-entry trajectory optimization (SpaceX Falcon 9),
 * long-duration optimal control with sensitive dynamics.
 */
NLSolution* shooting_multiple(BVPProblem *bvp, const double *x_guess,
                               size_t n_subintervals, size_t n_steps_per,
                               double tol, size_t max_iter);

/* ==========================================================================
 * L9: Advanced Nonlinear Methods
 * ========================================================================== */

/**
 * compute_equilibrium — find equilibrium of ẋ = f(x) by solving f(x) = 0
 *
 * L9: For nonlinear systems, equilibrium computation is the first step
 * in linearization-based control design. Uses Newton-Kantorovich internally.
 */
NLSolution* compute_equilibrium(NLEquation *eq, const double *x0,
                                 double tol, size_t max_iter);

/**
 * compute_limit_cycle_poincare — detect limit cycle via Poincaré section
 *
 * L9 Frontier: Nonlinear oscillations in control (e.g., relay feedback,
 * saturating actuators). Detects periodic orbits by solving for
 * fixed point of Poincaré map using shooting.
 */
NLSolution* compute_limit_cycle_poincare(NLEquation *eq, const double *x0,
                                          double period_guess, double tol);

/**
 * nl_solution_free — deallocate nonlinear solution
 */
void nl_solution_free(NLSolution *sol);

#ifdef __cplusplus
}
#endif

#endif /* NONLINEAR_METHODS_H */
