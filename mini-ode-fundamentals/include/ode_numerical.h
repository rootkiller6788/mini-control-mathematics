/**
 * @file    ode_numerical.h
 * @brief   Numerical integration methods for ODEs.
 *
 * L5: Computational Methods – single-step, multistep, adaptive, implicit.
 *
 * Reference: Butcher "Numerical Methods for Ordinary Differential Equations" (2016)
 *            Hairer, Nørsett, Wanner "Solving Ordinary Differential Equations I & II"
 *            Ascher & Petzold "Computer Methods for ODEs" (1998)
 *
 * Course Mapping:
 *   MIT 6.241 – Numerical simulation of dynamic systems
 *   Stanford ENGR 207B – Computational methods
 *   ETH 151-0591 – Numerical integration
 *   Caltech CDS 101 – Time stepping
 *   Tsinghua 数值分析 – Numerical ODE methods
 */

#ifndef ODE_NUMERICAL_H
#define ODE_NUMERICAL_H

#include "ode_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ────────────────────────── L5: Single-Step Methods ──────────────────── */

/**
 * @brief Forward Euler method (explicit, first order):
 *        y_{n+1} = y_n + h·f(t_n, y_n)
 *
 * Complexity: O(dim) per step.  Order: 1.  A-stable: No.
 *
 * @param ivp     IVP specification
 * @param t       Current time t_n
 * @param y       Current state y_n (length dim)
 * @param h       Step size
 * @param y_next  Output: y_{n+1} (length dim)
 * @return        0 on success
 */
int ode_euler_forward_step(const ODEIVP *ivp, double t, const double *y,
                             double h, double *y_next);

/**
 * @brief Backward (implicit) Euler method:
 *        y_{n+1} = y_n + h·f(t_{n+1}, y_{n+1})
 *
 * Solved via Newton iteration. A-stable, suitable for stiff problems.
 * Complexity: O(dim³·iter) per step due to linear solve.
 *
 * @param ivp     IVP (requires Jacobian for Newton iteration)
 * @param t       Current time t_n
 * @param y       Current state y_n
 * @param h       Step size
 * @param y_next  Output: y_{n+1}
 * @return        0 on success
 */
int ode_euler_backward_step(const ODEIVP *ivp, double t, const double *y,
                              double h, double *y_next);

/**
 * @brief Heun's method (Improved Euler, explicit, second order):
 *        k₁ = f(t_n, y_n)
 *        k₂ = f(t_n + h, y_n + h·k₁)
 *        y_{n+1} = y_n + (h/2)·(k₁ + k₂)
 *
 * This is the simplest two-stage explicit Runge-Kutta method.
 *
 * @param ivp     IVP specification
 * @param t       Current time
 * @param y       Current state
 * @param h       Step size
 * @param y_next  Output: next state
 * @return        0 on success
 */
int ode_heun_step(const ODEIVP *ivp, double t, const double *y,
                    double h, double *y_next);

/**
 * @brief Classical RK4 method (explicit, fourth order):
 *        k₁ = f(t_n, y_n)
 *        k₂ = f(t_n + h/2, y_n + h/2·k₁)
 *        k₃ = f(t_n + h/2, y_n + h/2·k₂)
 *        k₄ = f(t_n + h, y_n + h·k₃)
 *        y_{n+1} = y_n + (h/6)·(k₁ + 2k₂ + 2k₃ + k₄)
 *
 * Workhorse of ODE solving. Order 4, 4 stages, 4 function evaluations/step.
 *
 * @param ivp     IVP specification
 * @param t       Current time
 * @param y       Current state
 * @param h       Step size
 * @param y_next  Output: next state
 * @return        0 on success
 */
int ode_rk4_step(const ODEIVP *ivp, double t, const double *y,
                   double h, double *y_next);

/**
 * @brief RK45 Dormand-Prince adaptive step (explicit, 5th order).
 *
 * Uses 6 stages to produce 5th-order solution with embedded 4th-order
 * error estimate. Step size is adapted to meet tolerance.
 *
 * Reference: Dormand & Prince "A family of embedded Runge-Kutta formulae"
 *            (1980). This is the method behind MATLAB's ode45.
 *
 * @param ivp     IVP specification
 * @param t       Current time (updated on output)
 * @param y       Current state (updated on output)
 * @param h       Step size (updated on output with suggested next step)
 * @param rtol    Relative tolerance
 * @param atol    Absolute tolerance
 * @param step_info Output: step information
 * @return        0 on success (step may have been accepted or rejected)
 */
int ode_rk45_step(const ODEIVP *ivp, double *t, double *y, double *h,
                    double rtol, double atol, ODEStep *step_info);

/* ────────────────────────── L5: Multistep Methods ──────────────────── */

/**
 * @brief Adams-Bashforth 2-step method (explicit, O(h²)):
 *        y_{n+1} = y_n + (h/2)·(3·f(t_n, y_n) - f(t_{n-1}, y_{n-1}))
 *
 * Uses 1 previous derivative value (2-step). Requires starting method.
 * Complexity: O(dim) per step (1 function evaluation).
 *
 * @param ivp     IVP specification
 * @param t       Current time t_n
 * @param y       Current state y_n
 * @param f_prev  Previous derivative f(t_{n-1}, y_{n-1})
 * @param h       Step size
 * @param y_next  Output: y_{n+1}
 * @return        0 on success
 */
int ode_ab2_step(const ODEIVP *ivp, double t, const double *y,
                   const double *f_prev, double h, double *y_next);

/**
 * @brief Adams-Bashforth 3-step method (explicit, O(h³)):
 *        y_{n+1} = y_n + (h/12)·(23f_n - 16f_{n-1} + 5f_{n-2})
 */
int ode_ab3_step(const ODEIVP *ivp, double t, const double *y,
                   const double *fn, const double *fn1, const double *fn2,
                   double h, double *y_next);

/**
 * @brief Adams-Bashforth 4-step method (explicit, O(h⁴)):
 *        y_{n+1} = y_n + (h/24)·(55f_n - 59f_{n-1} + 37f_{n-2} - 9f_{n-3})
 */
int ode_ab4_step(const ODEIVP *ivp, double t, const double *y,
                   const double *fn, const double *fn1,
                   const double *fn2, const double *fn3,
                   double h, double *y_next);

/**
 * @brief Adams-Moulton 2-step (implicit trapezoidal rule, O(h²)):
 *        y_{n+1} = y_n + (h/2)·(f(t_{n+1}, y_{n+1}) + f(t_n, y_n))
 *
 * Solved via simple fixed-point iteration (or Newton for stiff).
 */
int ode_am2_step(const ODEIVP *ivp, double t, const double *y,
                   double h, double *y_next, int max_iter, double tol);

/**
 * @brief Adams-Moulton 3-step (implicit, O(h³)):
 *        y_{n+1} = y_n + (h/24)·(9f_{n+1} + 19f_n - 5f_{n-1} + f_{n-2})
 */
int ode_am3_step(const ODEIVP *ivp, double t, const double *y,
                   const double *fn, const double *fn1, const double *fn2,
                   double h, double *y_next, int max_iter, double tol);

/* ────────────────────────── L5: BDF Methods (Stiff) ────────────────── */

/**
 * @brief BDF-2 (Backward Differentiation Formula, order 2, A-stable):
 *        y_{n+1} = (4/3)y_n - (1/3)y_{n-1} + (2h/3)·f(t_{n+1}, y_{n+1})
 *
 * The workhorse for stiff ODEs. A-stable up to order 2.
 * Solved via Newton iteration.
 *
 * @param ivp     IVP (requires Jacobian)
 * @param t       Current time t_n
 * @param y       Current state
 * @param y_prev  Previous state y_{n-1}
 * @param h       Step size
 * @param y_next  Output: y_{n+1}
 * @return        0 on success
 */
int ode_bdf2_step(const ODEIVP *ivp, double t, const double *y,
                    const double *y_prev, double h, double *y_next);

/**
 * @brief BDF-3 (order 3, A(α)-stable with α≈86°).
 */
int ode_bdf3_step(const ODEIVP *ivp, double t, const double *y,
                    const double *y_prev, const double *y_prev2,
                    double h, double *y_next);

/* ────────────────────────── L5: Full Integration ────────────────────── */

/**
 * @brief Integrate an IVP using a fixed-step method.
 *
 * Wraps any single-step method and produces a full solution trajectory.
 *
 * @param ivp     IVP specification
 * @param method  Integration method
 * @param h       Fixed step size
 * @param sol     Output: ODESolution (caller must free with ode_solution_free)
 * @return        0 on success
 */
int ode_integrate_fixed_step(const ODEIVP *ivp, ODEMethod method,
                               double h, ODESolution *sol);

/**
 * @brief Integrate an IVP using adaptive step size control.
 *
 * Uses RK45 with error estimation to adapt step size.
 *
 * @param ivp     IVP specification
 * @param config  Simulation configuration (rtol, atol, h_init, etc.)
 * @param sol     Output: ODESolution
 * @return        0 on success
 */
int ode_integrate_adaptive(const ODEIVP *ivp, const ODESimConfig *config,
                             ODESolution *sol);

/**
 * @brief General-purpose ODE integrator that selects method automatically.
 *
 * Chooses between explicit/implicit, fixed/adaptive based on problem.
 *
 * @param ivp     IVP specification
 * @param config  Configuration (may be NULL for defaults)
 * @param sol     Output: ODESolution
 * @return        0 on success
 */
int ode_solve(const ODEIVP *ivp, const ODESimConfig *config,
                ODESolution *sol);

/* ────────────────────────── L7: Stiffness Detection ────────────────── */

/**
 * @brief Estimate the stiffness ratio of a linear ODE system.
 *
 * Stiffness ratio = max|Re(λ)| / min|Re(λ)| for stable λ.
 * If > 1000, the system is considered stiff and implicit methods
 * should be preferred.
 *
 * Reference: Hairer & Wanner "Solving ODEs II" (1996), Stiff Problems
 *
 * @param sys      Linear system
 * @return         Stiffness ratio (≥ 1; 0 on error)
 */
double ode_stiffness_ratio(const LinearODESystem *sys);

/**
 * @brief Determine whether an IVP is stiff and recommend an integration
 *        method. Considers eigenvalue spectrum, nonlinearity, and tolerances.
 *
 * @param ivp            IVP
 * @param recommended    Output: recommended method
 * @return               true if system appears stiff
 */
bool ode_detect_stiffness(const ODEIVP *ivp, ODEMethod *recommended);

/* ────────────────────────── BVP Solvers ────────────────────────────── */

/**
 * @brief Linear shooting method for two-point boundary value problem:
 *        y'' = p(t)y' + q(t)y + r(t),  y(a) = α,  y(b) = β
 *
 * Converts BVP to two IVPs and uses superposition.
 *
 * Reference: Burden & Faires "Numerical Analysis" (2011)
 *
 * @param p_vals   Coefficient p(t) at grid points (length n_pts)
 * @param q_vals   Coefficient q(t) at grid points
 * @param r_vals   Forcing r(t) at grid points
 * @param t_vals   Grid points (length n_pts)
 * @param n_pts    Number of grid points
 * @param alpha    Left boundary condition y(a)
 * @param beta     Right boundary condition y(b)
 * @param y_out    Solution y(t) (length n_pts)
 * @return         0 on success
 */
int ode_bvp_shooting_linear(const double *p_vals, const double *q_vals,
                              const double *r_vals, const double *t_vals,
                              int n_pts, double alpha, double beta,
                              double *y_out);

/**
 * @brief Nonlinear shooting method using Newton iteration:
 *        y'' = f(t, y, y'),  y(a) = α, y(b) = β
 *
 * @param f_func       RHS f(t, y, y')
 * @param fy_func      ∂f/∂y
 * @param fyp_func     ∂f/∂y'
 * @param ctx          User context
 * @param t_vals       Grid (length n_pts)
 * @param n_pts        Number of points
 * @param alpha        y(a)
 * @param beta         y(b)
 * @param y_out        Solution (length n_pts)
 * @return             0 on success
 */
int ode_bvp_shooting_nonlinear(
    double (*f_func)(double t, double y, double yp, void *ctx),
    double (*fy_func)(double t, double y, double yp, void *ctx),
    double (*fyp_func)(double t, double y, double yp, void *ctx),
    void *ctx, const double *t_vals, int n_pts,
    double alpha, double beta, double *y_out);

/**
 * @brief Linear finite difference method for BVP:
 *        y'' = p(t)y' + q(t)y + r(t),  y(a)=α, y(b)=β
 *
 * Solves a tridiagonal system. O(h²) accurate.
 *
 * @param p_vals   p(t) at interior grid points (length n_int)
 * @param q_vals   q(t) at interior grid points
 * @param r_vals   r(t) at interior grid points
 * @param a        Left boundary
 * @param b        Right boundary
 * @param alpha    y(a)
 * @param beta     y(b)
 * @param n_int    Number of interior points
 * @param t_vals   Full grid (output, length n_int+2)
 * @param y_vals   Solution (output, length n_int+2)
 * @return         0 on success
 */
int ode_bvp_fd_linear(const double *p_vals, const double *q_vals,
                        const double *r_vals, double a, double b,
                        double alpha, double beta, int n_int,
                        double *t_vals, double *y_vals);

/* ────────────────────────── Utility Functions ──────────────────────── */

/**
 * @brief Estimate the local truncation error using Richardson extrapolation
 *        (step-halving):  err ≈ (y_h - y_{h/2}) / (2^p - 1)
 *
 * @param ivp       IVP specification
 * @param method    Integration method step function
 * @param t         Current time
 * @param y         Current state
 * @param h         Step size
 * @param order     Method order p
 * @param err_est   Output: estimated error
 * @return          0 on success
 */
int ode_richardson_error(const ODEIVP *ivp, ODEMethod method,
                           double t, const double *y, double h,
                           int order, double *err_est);

#ifdef __cplusplus
}
#endif

#endif /* ODE_NUMERICAL_H */
