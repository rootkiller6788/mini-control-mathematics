/**
 * integration.h — ODE Integration Methods for Control Systems
 *
 * Covers: L4 Conservation Laws (ODE forms), L5 Engineering Methods (integrators)
 *
 * State-space control systems are governed by ODEs: dx/dt = f(t, x, u).
 * Numerical integration is the computational backbone of simulation,
 * trajectory prediction, and model predictive control (MPC).
 *
 * Reference: Hairer, Norsett, Wanner, Solving Ordinary Differential Equations I (1993)
 *            Shampine & Gordon, Computer Solution of Ordinary Differential Equations (1975)
 *            Ascher & Petzold, Computer Methods for ODEs and DAEs (1998)
 */

#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "numerical_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: ODE System Definition
 * ========================================================================== */

/**
 * ODESystem — defines an initial value problem (IVP)
 *
 * L1 Definition: An ODE system is dx/dt = f(t, x) with x(t₀) = x₀.
 * Here f: ℝ × ℝⁿ → ℝⁿ is the vector field (right-hand side).
 */
typedef struct {
    /** Right-hand side function: dx/dt = rhs(t, state, params, n) */
    void (*rhs)(double t, const double *state, double *deriv, void *params, size_t n);
    void   *params;    /**< User-defined parameters passed to rhs */
    size_t  n;         /**< State dimension */
} ODESystem;

/**
 * ODESolution — stores the trajectory of an ODE solution
 *
 * L1 Definition: A trajectory is the set of points {(tₖ, xₖ)} for k = 0,...,N.
 */
typedef struct {
    double  *time;     /**< Array of N+1 time points */
    double  *states;   /**< Array of (N+1) × n state vectors, row-major */
    size_t   n_states; /**< State dimension */
    size_t   n_steps;  /**< Number of time steps (= N) */
} ODESolution;

/* ==========================================================================
 * L4: ODE as Conservation Form
 * ========================================================================== */

/**
 * ODEStatistics — integration quality metrics
 *
 * L3 Engineering Quantities: number of function evaluations, step rejections,
 * local error estimates — these quantify computational cost and accuracy.
 */
typedef struct {
    size_t n_f_eval;           /**< Number of RHS evaluations */
    size_t n_steps_total;      /**< Total steps attempted */
    size_t n_steps_rejected;   /**< Rejected steps (adaptive methods) */
    double final_local_error;  /**< Estimated local error at final step */
} ODEStatistics;

/* ==========================================================================
 * L5: Single-Step Methods (Runge-Kutta Family)
 * ========================================================================== */

/**
 * euler_forward_step — single step of Forward Euler: x_{n+1} = x_n + h·f(t_n, x_n)
 *
 * L5 Method: Explicit Euler, 1st-order, conditionally stable.
 * Complexity: O(n) per step, 1 RHS eval per step.
 * Stability: |1 + hλ| ≤ 1 where λ are eigenvalues of ∂f/∂x.
 */
int euler_forward_step(ODESystem *sys, double t, const double *x,
                       double h, double *x_next);

/**
 * euler_backward_step — single step of Implicit Euler: x_{n+1} = x_n + h·f(t_{n+1}, x_{n+1})
 *
 * L5 Method: Implicit Euler, 1st-order, A-stable (unconditionally stable for linear).
 * Uses Newton iteration to solve implicit equation.
 * Complexity: O(k·n) per step (k = Newton iterations).
 */
int euler_backward_step(ODESystem *sys, double t, const double *x,
                        double h, double *x_next, size_t max_newton, double tol);

/**
 * heun_step — Improved Euler (Heun's method, explicit trapezoidal)
 * k1 = f(t, x), k2 = f(t+h, x+h·k1), x_{n+1} = x_n + (h/2)·(k1 + k2)
 *
 * L5 Method: 2nd-order, 2-stage explicit RK. Local error O(h³).
 */
int heun_step(ODESystem *sys, double t, const double *x,
              double h, double *x_next);

/**
 * rk4_step — Classical 4th-order Runge-Kutta
 *
 * L5 Method: The workhorse of ODE integration. 4 stages, O(h⁵) local error.
 * k1 = f(t, x)
 * k2 = f(t + h/2, x + h·k1/2)
 * k3 = f(t + h/2, x + h·k2/2)
 * k4 = f(t + h, x + h·k3)
 * x_{n+1} = x_n + (h/6)·(k1 + 2k2 + 2k3 + k4)
 */
int rk4_step(ODESystem *sys, double t, const double *x,
             double h, double *x_next);

/**
 * rkf45_step — Runge-Kutta-Fehlberg 5(4) adaptive step
 *
 * L5/L8 Method: Embedded RK pair of orders 5 and 4. Produces both 5th-order
 * solution and 4th-order error estimate. Used for adaptive step size control.
 * Returns suggested next step size in *h_next.
 */
int rkf45_step(ODESystem *sys, double t, const double *x,
               double h, double *x_next, double *error_est,
               double *h_next, double rtol, double atol);

/* ==========================================================================
 * L5/L8: Multi-Step Methods
 * ========================================================================== */

/**
 * adams_bashforth2_step — 2-step Adams-Bashforth
 *
 * L5 Method: x_{n+2} = x_{n+1} + h·(3/2·f_{n+1} - 1/2·f_n)
 * Explicit linear multistep, 2nd order. Requires previous derivative.
 */
int adams_bashforth2_step(ODESystem *sys, double t, double h,
                          const double *x_prev, const double *f_prev,
                          const double *x_curr, double *x_next);

/**
 * adams_bashforth3_step — 3-step Adams-Bashforth
 *
 * L5 Method: x_{n+3} = x_{n+2} + h·(23/12·f_{n+2} - 16/12·f_{n+1} + 5/12·f_n)
 * 3rd order, requires 3 previous derivatives.
 */
int adams_bashforth3_step(ODESystem *sys, double t, double h,
                          const double *f_n, const double *f_n1,
                          const double *f_n2, const double *x_curr,
                          double *x_next);

/**
 * adams_bashforth4_step — 4-step Adams-Bashforth
 *
 * L5 Method: 4th order explicit multistep.
 * Coefficients: 55/24, -59/24, 37/24, -9/24
 */
int adams_bashforth4_step(ODESystem *sys, double t, double h,
                          const double *f_n, const double *f_n1,
                          const double *f_n2, const double *f_n3,
                          const double *x_curr, double *x_next);

/* ==========================================================================
 * L5: Full Integration Drivers
 * ========================================================================== */

/**
 * solve_ode_fixed_step — integrate ODE over [t0, tf] with fixed step size h
 *
 * L5 Method: Uses RK4 internally. Returns full trajectory.
 * Complexity: O(N·n) where N = (tf - t0) / h.
 */
ODESolution* solve_ode_fixed_step(ODESystem *sys, double t0, double tf,
                                  const double *x0, double h);

/**
 * solve_ode_adaptive — integrate ODE with adaptive step size control (RKF45)
 *
 * L5/L8 Method: Adaptive step size based on local error estimation.
 * Complexity: variable, depends on stiffness and tolerance.
 */
ODESolution* solve_ode_adaptive(ODESystem *sys, double t0, double tf,
                                const double *x0, double h_init,
                                double rtol, double atol,
                                ODEStatistics *stats);

/**
 * odesolution_free — deallocate ODE solution
 */
void odesolution_free(ODESolution *sol);

/**
 * odesolution_interpolate — evaluate solution at arbitrary time via linear interpolation
 */
int odesolution_interpolate(const ODESolution *sol, double t, double *x_out);

#ifdef __cplusplus
}
#endif

#endif /* INTEGRATION_H */
