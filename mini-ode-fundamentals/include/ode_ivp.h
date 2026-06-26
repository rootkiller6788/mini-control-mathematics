/**
 * @file    ode_ivp.h
 * @brief   Initial value problem formulation and fundamental theory.
 *
 * L1-L4: Definitions, existence/uniqueness, Lipschitz conditions,
 *        Picard iteration, Gronwall inequality.
 *
 * Reference: Coddington & Levinson "Theory of ODEs" (1955)
 *            Hairer, Nørsett, Wanner "Solving ODEs I" (1993)
 *
 * Course Mapping:
 *   MIT 6.241 – Dynamic systems IVP formulation
 *   Stanford ENGR 207B – Existence and uniqueness
 *   Cambridge 3F2 – Well-posedness of IVPs
 *   Caltech CDS 101 – Solution fundamentals
 */

#ifndef ODE_IVP_H
#define ODE_IVP_H

#include "ode_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ──────────────────────── L1: IVP Construction ──────────────────────── */

/**
 * @brief Create an IVP specification.
 *
 * @param f       RHS function
 * @param jac     Jacobian (may be NULL)
 * @param ctx     User context (may be NULL)
 * @param t0      Initial time
 * @param t_end   Final time
 * @param y0      Initial condition array (length dim)
 * @param dim     System dimension
 * @param ode_class  Classification
 * @return        Initialized ODEIVP struct
 */
ODEIVP ode_ivp_create(ODERHSFunc f, ODEJacobianFunc jac, void *ctx,
                       double t0, double t_end, const double *y0,
                       int dim, ODEClass ode_class);

/**
 * @brief Validate IVP for well-posedness (check dimension consistency,
 *        finite values, time ordering).
 *
 * @param ivp  IVP to validate
 * @return     true if IVP is well-formed
 */
bool ode_ivp_validate(const ODEIVP *ivp);

/**
 * @brief Free resources associated with a solution.
 */
void ode_solution_free(ODESolution *sol);

/**
 * @brief Clone an IVP (deep copy of y0).
 *
 * @param ivp  Original IVP
 * @return     Cloned IVP (caller must free y0 with free())
 */
ODEIVP ode_ivp_clone(const ODEIVP *ivp);

/**
 * @brief Destroy IVP resources (frees allocated y0).
 */
void ode_ivp_destroy(ODEIVP *ivp);

/* ──────────────────────── L2: Lipschitz Condition ───────────────────── */

/**
 * @brief Estimate the Lipschitz constant L of f(t,y) w.r.t. y
 *        over a domain using finite difference sampling.
 *
 * The Lipschitz condition: ||f(t,y₁) - f(t,y₂)|| ≤ L·||y₁ - y₂||
 *
 * This is essential for the Picard-Lindelöf theorem (L4).
 *
 * @param ivp         IVP specification
 * @param y_center    Center of sampling domain
 * @param radius      Sampling radius around y_center
 * @param n_samples   Number of random samples in ball
 * @return            Estimated Lipschitz constant (0 if cannot estimate)
 */
double ode_lipschitz_estimate(const ODEIVP *ivp, const double *y_center,
                                double radius, int n_samples);

/**
 * @brief Check if the RHS satisfies a one-sided Lipschitz condition
 *        (useful for contractivity analysis).
 *
 * <f(t,y₁) - f(t,y₂), y₁ - y₂> ≤ ν·||y₁ - y₂||²
 *
 * @param ivp         IVP specification
 * @param y_center    Center of sampling domain
 * @param radius      Sampling radius
 * @param n_samples   Number of samples
 * @param nu_out      Output: one-sided Lipschitz constant
 * @return            true if condition holds
 */
bool ode_one_sided_lipschitz(const ODEIVP *ivp, const double *y_center,
                               double radius, int n_samples, double *nu_out);

/* ──────────────────────── L4: Picard Iteration ──────────────────────── */

/**
 * @brief Perform one step of Picard iteration:
 *        y_{k+1}(t) = y₀ + ∫_{t₀}^{t} f(s, y_k(s)) ds
 *
 * Reference: Theorem of Picard-Lindelöf (1890)
 *   If f is continuous in t and Lipschitz in y on a rectangle
 *   R = { (t,y) : |t-t₀| ≤ a, ||y-y₀|| ≤ b }, then there exists
 *   a unique solution on [t₀-α, t₀+α] where α = min(a, b/M),
 *   M = max_{R} ||f(t,y)||.
 *
 * @param ivp          IVP specification (used for RHS)
 * @param y_prev       Previous iterate y_k(t) as function values
 * @param t_vals       Time discretization (length n_points)
 * @param n_points     Number of time points
 * @param y_next       Output: next iterate y_{k+1}(t) (length: n_points*dim)
 * @return             Maximum change ||y_{k+1} - y_k||_∞
 */
double ode_picard_iteration(const ODEIVP *ivp,
                             const double *y_prev,
                             const double *t_vals,
                             int n_points,
                             double *y_next);

/**
 * @brief Solve IVP via Picard iteration until convergence.
 *
 * @param ivp          IVP specification
 * @param t_vals       Time points (length n_points)
 * @param n_points     Number of points
 * @param max_iter     Maximum iterations
 * @param rtol         Relative convergence tolerance
 * @param y_out        Output solution (length: n_points*dim)
 * @return             Number of iterations (0 if failed)
 */
int ode_picard_solve(const ODEIVP *ivp, const double *t_vals, int n_points,
                      int max_iter, double rtol, double *y_out);

/* ──────────────────────── L4: Gronwall's Inequality ──────────────────── */

/**
 * @brief Gronwall's inequality (differential form):
 *        If u'(t) ≤ β(t)·u(t), then u(t) ≤ u(0)·exp(∫₀ᵗ β(s) ds).
 *
 * Used to bound solution growth and prove continuous dependence on initial data.
 *
 * @param u0         Initial value u(0)
 * @param beta_vals  β(t) values (length n_points)
 * @param t_vals     Time points (length n_points)
 * @param n_points   Number of points
 * @param u_bound    Output: upper bound for u(t) at each t
 * @return           0 on success
 */
int ode_gronwall_bound(double u0, const double *beta_vals,
                        const double *t_vals, int n_points,
                        double *u_bound);

/**
 * @brief Verify continuous dependence on initial conditions:
 *        ||y(t; y₀₁) - y(t; y₀₂)|| ≤ ||y₀₁ - y₀₂||·exp(L·|t-t₀|)
 *
 * where L is the Lipschitz constant. This follows from Gronwall.
 *
 * @param ivp          IVP with Lipschitz constant L
 * @param L            Lipschitz constant
 * @param y0_a         First initial condition
 * @param y0_b         Second initial condition
 * @param t            Time to evaluate bound
 * @return             Upper bound on ||y(t; y₀₁) - y(t; y₀₂)||
 */
double ode_continuous_dependence_bound(const ODEIVP *ivp, double L,
                                         const double *y0_a,
                                         const double *y0_b, double t);

/* ──────────────────────── Default Simulation Config ──────────────────── */

/**
 * @brief Get a sensible default simulation configuration for the given IVP.
 *
 * Chooses method based on ode_class, sets reasonable tolerances.
 */
ODESimConfig ode_simconfig_default(const ODEIVP *ivp);

/* ──────────────────────── Solution Utility ──────────────────────────── */

/**
 * @brief Print a summary of the ODE solution to stdout.
 */
void ode_solution_print(const ODESolution *sol);

/**
 * @brief Interpolate solution at an arbitrary time using linear interpolation.
 *
 * @param sol    Solution
 * @param t      Query time (must be within [t0, t_end])
 * @param y_out  Interpolated state (length: sol->dim)
 * @return       0 on success, -1 if t out of range
 */
int ode_solution_interpolate(const ODESolution *sol, double t, double *y_out);

/**
 * @brief Compute L2 norm of solution error (if reference solution known).
 */
double ode_solution_error_norm(const ODESolution *computed,
                                 const ODESolution *reference);

#ifdef __cplusplus
}
#endif

#endif /* ODE_IVP_H */
