/**
 * @file    ode_stability.h
 * @brief   Stability analysis for ODEs: Lyapunov methods, stability of
 *          equilibria, instability criteria, absolute stability of methods.
 *
 * L1-L2, L4, L8: Stability definitions, Lyapunov theorems, advanced topics.
 *
 * Reference: Khalil "Nonlinear Systems" (2002)
 *            Lyapunov "The General Problem of the Stability of Motion" (1892)
 *            LaSalle & Lefschetz "Stability by Liapunov's Direct Method" (1961)
 *            Chetaev "The Stability of Motion" (1961)
 *
 * Course Mapping:
 *   MIT 6.243 – Lyapunov stability theory
 *   Stanford ENGR 209A – Nonlinear stability
 *   Berkeley ME 234 – Stability analysis
 *   ETH 151-0563 – Robust stability
 *   Cambridge 4F3 – Nonlinear stability
 */

#ifndef ODE_STABILITY_H
#define ODE_STABILITY_H

#include "ode_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ──────────────────────── L1: Stability Definitions ─────────────────── */

/**
 * @brief Check Lyapunov stability of an equilibrium:
 *        ∀ ε > 0, ∃ δ > 0 such that ||x(0)|| < δ ⇒ ||x(t)|| < ε ∀t ≥ 0.
 *
 * This function checks numerically for a given equilibrium.
 *
 * @param ivp       IVP specification (equilibrium at origin assumed)
 * @param epsilon   Perturbation radius ε
 * @param n_trials  Number of random initial conditions within δ-ball
 * @param T_max     Simulation time horizon
 * @param delta_out Output: feasible δ (0 if unstable)
 * @return          true if appears Lyapunov stable
 */
bool ode_check_lyapunov_stability(const ODEIVP *ivp, double epsilon,
                                    int n_trials, double T_max,
                                    double *delta_out);

/**
 * @brief Check asymptotic stability: the equilibrium is Lyapunov stable
 *        AND x(t) → 0 as t → ∞ for all ||x(0)|| < δ.
 *
 * @param ivp       IVP (equilibrium at origin)
 * @param delta     Initial perturbation bound
 * @param n_trials  Number of trials
 * @param T_max     Time horizon
 * @param decay_tol Tolerance for detecting convergence to origin
 * @return          true if appears asymptotically stable
 */
bool ode_check_asymptotic_stability(const ODEIVP *ivp, double delta,
                                      int n_trials, double T_max,
                                      double decay_tol);

/**
 * @brief Check exponential stability:
 *        ||x(t)|| ≤ M·exp(-α·t)·||x(0)|| for some M ≥ 1, α > 0.
 *
 * @param ivp       IVP (equilibrium at origin)
 * @param delta     Initial perturbation bound
 * @param n_trials  Number of trials
 * @param T_max     Time horizon
 * @param alpha_out Output: estimated decay rate α
 * @param M_out     Output: estimated overshoot bound M
 * @return          true if appears exponentially stable
 */
bool ode_check_exponential_stability(const ODEIVP *ivp, double delta,
                                       int n_trials, double T_max,
                                       double *alpha_out, double *M_out);

/* ──────────────────────── L4: Lyapunov's Direct Method ──────────────── */

/**
 * @brief Lyapunov's direct (second) method: test if a candidate function
 *        V(x) is a valid Lyapunov function by checking:
 *        1) V(0) = 0
 *        2) V(x) > 0 for all x ≠ 0 in some neighborhood (positive definite)
 *        3) V̇(x) = ∇V·f(x) < 0 for all x ≠ 0 (negative definite derivative)
 *
 * Reference: Lyapunov (1892), trans. Fuller (1992)
 *
 * @param vf      Vector field (autonomous assumed)
 * @param V       Candidate Lyapunov function
 * @param radius  Neighborhood radius to check
 * @param n_pts   Number of random sample points
 * @param V_pos   Output: true if V confirmed positive definite in ball
 * @param Vdot_neg Output: true if V̇ confirmed negative definite in ball
 * @return        true if V is a valid Lyapunov function (both conditions hold)
 */
bool ode_lyapunov_validate(const VectorField *vf, const LyapunovFunction *V,
                             double radius, int n_pts,
                             bool *V_pos, bool *Vdot_neg);

/**
 * @brief Compute V̇ = ∇V(x)·f(x), the orbital derivative of V along
 *        trajectories of the vector field f.
 *
 * @param vf       Vector field
 * @param V        Lyapunov function
 * @param x        Point in state space (length dim)
 * @param Vdot_out Output: V̇(x)
 * @return         0 on success
 */
int ode_lyapunov_derivative(const VectorField *vf, const LyapunovFunction *V,
                              const double *x, double *Vdot_out);

/* ──────────────────────── L4: Lyapunov Equation (Linear) ────────────── */

/**
 * @brief Solve the continuous-time Lyapunov equation: A·P + P·Aᵀ = -Q
 *
 * For a linear system ẋ = Ax, if P ≻ 0 exists for any Q ≻ 0, then
 * the system is asymptotically stable and V(x) = xᵀPx is a Lyapunov function.
 *
 * Uses Bartels-Stewart algorithm (Schur decomposition).
 *
 * Complexity: O(dim³)
 *
 * Reference: Bartels & Stewart "Solution of the Matrix Equation AX + XB = C"
 *            (1972). This is a fundamental result in linear control theory.
 *
 * Course Mapping: MIT 6.241 (Lyapunov equation), EE221A (stability analysis)
 *
 * @param A       System matrix (dim×dim)
 * @param Q       Right-hand side matrix (must be symmetric positive definite)
 * @param dim     Dimension
 * @param P       Output: solution matrix P (symmetric if Q is symmetric)
 * @return        0 on success (P positive definite = system stable)
 */
int ode_lyapunov_equation(const double *A, const double *Q, int dim,
                            double *P);

/**
 * @brief Construct a quadratic Lyapunov function V(x) = xᵀPx for
 *        linear system ẋ = Ax by solving the Lyapunov equation AᵀP + PA = -I.
 *
 * @param A       System matrix
 * @param dim     Dimension
 * @param V       Output: Lyapunov function (P stored in ctx as dim×dim matrix)
 * @return        true if a positive definite P was found (system stable)
 */
bool ode_construct_quadratic_lyapunov(const double *A, int dim,
                                        LyapunovFunction *V);

/* ──────────────────────── L4: LaSalle's Invariance Principle ────────── */

/**
 * @brief LaSalle's invariance principle: For an autonomous system with
 *        V(x) > 0 and V̇(x) ≤ 0, every trajectory approaches the largest
 *        invariant set in {x : V̇(x) = 0}.
 *
 * This function identifies the set where V̇(x) = 0 and checks if it
 * contains only the equilibrium (global asymptotic stability).
 *
 * Reference: LaSalle "Some Extensions of Liapunov's Second Method" (1960)
 *
 * @param vf        Vector field
 * @param V         Lyapunov function with V̇ ≤ 0
 * @param radius    Search radius
 * @param n_pts     Number of check points
 * @param only_zero Output: true if V̇=0 only at origin
 * @return          0 on success
 */
int ode_lasalle_check(const VectorField *vf, const LyapunovFunction *V,
                        double radius, int n_pts, bool *only_zero);

/* ──────────────────────── L4: Chetaev Instability Theorem ───────────── */

/**
 * @brief Chetaev's instability theorem: If there exists a function V(x)
 *        and a region Ω where V > 0 and V̇ > 0 for x ≠ 0, with the origin
 *        on the boundary of Ω, then the origin is unstable.
 *
 * @param vf       Vector field
 * @param V        Candidate Chetaev function
 * @param n_pts    Number of sample points
 * @param unstable Output: true if instability is indicated
 * @return         0 on success
 */
int ode_chetaev_check(const VectorField *vf, const LyapunovFunction *V,
                        int n_pts, bool *unstable);

/* ──────────────────────── L8: Region of Attraction ──────────────────── */

/**
 * @brief Estimate the region of attraction (ROA) for an asymptotically
 *        stable equilibrium using Lyapunov sublevel sets.
 *
 * The ROA contains the largest sublevel set {x: V(x) ≤ c} that is
 * entirely contained in the region where V̇(x) < 0.
 *
 * Reference: Khalil (2002), Chap. 8
 *
 * @param vf       Vector field
 * @param V        Lyapunov function
 * @param radius   Maximum search radius
 * @param n_pts    Sampling points for level set search
 * @param c_max    Output: largest verified sublevel set constant
 * @return         0 on success
 */
int ode_estimate_roa(const VectorField *vf, const LyapunovFunction *V,
                       double radius, int n_pts, double *c_max);

/* ──────────────────────── L1: Stability of Numerical Methods ────────── */

/**
 * @brief Determine the absolute stability region boundary for an explicit
 *        Runge-Kutta method applied to the test equation ẏ = λy.
 *
 * The stability function R(z) where z = hλ must satisfy |R(z)| ≤ 1.
 * For RK4: R(z) = 1 + z + z²/2 + z³/6 + z⁴/24.
 *
 * @param method     ODE method
 * @param z_real     Array of real parts (output)
 * @param z_imag     Array of imaginary parts (output)
 * @param n_pts      Number of boundary points
 * @param min_real   Output: minimum real z where |R(z)| ≤ 1
 * @return           0 on success
 */
int ode_absolute_stability_region(ODEMethod method,
                                    double *z_real, double *z_imag,
                                    int n_pts, double *min_real);

/**
 * @brief Check time step satisfies linear stability condition:
 *        h < h_max where h_max = stability_boundary / max|λ|.
 *
 * @param A       System matrix (for eigenvalue computation)
 * @param dim     Dimension
 * @param method  Integration method
 * @param h_max   Output: maximum stable step size
 * @return        0 on success
 */
int ode_stable_step_size(const double *A, int dim, ODEMethod method,
                           double *h_max);

/* ──────────────────────── L8: ISS (Input-to-State Stability) ────────── */

/**
 * @brief Check input-to-state stability (ISS) for a perturbed system
 *        ẋ = f(x, u) where u(t) is a bounded disturbance.
 *
 * ISS: There exist β∈KL, γ∈K∞ such that
 *      ||x(t)|| ≤ β(||x(0)||, t) + γ(||u||_∞).
 *
 * This is a fundamental concept in nonlinear control (Sontag, 1989).
 *
 * @param ivp        IVP with disturbance
 * @param u_bound    Bound on disturbance ||u||_∞
 * @param n_trials   Number of trial disturbances
 * @param T_max      Time horizon
 * @param iss_gain   Output: estimate of ISS gain γ
 * @return           true if ISS property appears to hold
 */
bool ode_check_iss(const ODEIVP *ivp, double u_bound,
                     int n_trials, double T_max, double *iss_gain);

#ifdef __cplusplus
}
#endif

#endif /* ODE_STABILITY_H */
