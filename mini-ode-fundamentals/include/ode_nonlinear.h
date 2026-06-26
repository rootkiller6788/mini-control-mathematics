/**
 * @file    ode_nonlinear.h
 * @brief   Nonlinear ODE analysis: phase plane, nullclines, equilibrium
 *          classification, limit cycles, bifurcations.
 *
 * L2, L4, L6, L8: Nonlinear concepts, canonical systems, advanced topics.
 *
 * Reference: Strogatz "Nonlinear Dynamics and Chaos" (2015)
 *            Wiggins "Introduction to Applied Nonlinear Dynamical Systems" (2003)
 *            Kuznetsov "Elements of Applied Bifurcation Theory" (2004)
 *
 * Course Mapping:
 *   MIT 6.243 – Nonlinear control systems
 *   Stanford ENGR 209A – Nonlinear dynamics
 *   Berkeley ME 234 – Nonlinear systems analysis
 *   Caltech CDS 240 – Nonlinear dynamics and bifurcations
 *   Cambridge 4F3 – Nonlinear and predictive control
 */

#ifndef ODE_NONLINEAR_H
#define ODE_NONLINEAR_H

#include "ode_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ──────────────────────── L2: Equilibrium Analysis ──────────────────── */

/**
 * @brief Find equilibrium points of an autonomous system ẏ = f(y)
 *        by solving f(y) = 0 using Newton's method with multiple
 *        random initial guesses.
 *
 * @param vf         Vector field
 * @param init_guess Initial guesses for Newton (array of guess*dim doubles)
 * @param n_guesses  Number of initial guesses
 * @param tol        Convergence tolerance
 * @param max_iter   Maximum Newton iterations per guess
 * @param eq_points  Output: equilibrium points found (length ≤ n_guesses*dim)
 * @param n_found    Output: number of distinct equilibria found
 * @return           0 on success
 */
int ode_find_equilibria(const VectorField *vf, const double *init_guess,
                          int n_guesses, double tol, int max_iter,
                          double *eq_points, int *n_found);

/**
 * @brief Determine the type of a 2D equilibrium via linearization.
 *
 * Computes the Jacobian J = ∂f/∂y at the equilibrium, its eigenvalues,
 * and classifies the equilibrium type.
 *
 * @param vf        Vector field (must be 2D)
 * @param eq_point  The equilibrium point (length 2)
 * @param class_out Output: equilibrium classification
 * @return          0 on success
 */
int ode_classify_equilibrium_2d_vf(const VectorField *vf,
                                     const double *eq_point,
                                     EquilibriaClassification *class_out);

/* ──────────────────────── L2: Phase Plane Analysis ──────────────────── */

/**
 * @brief Compute nullclines for a 2D autonomous system:
 *        x-nullcline: f₁(x, y) = 0
 *        y-nullcline: f₂(x, y) = 0
 *
 * Nullclines partition the phase plane into regions where direction
 * of flow is monotonic. Intersections are equilibrium points.
 *
 * @param vf          Vector field (2D)
 * @param portrait    Output: PhasePortrait with nullcline data
 * @return            0 on success
 */
int ode_compute_nullclines(const VectorField *vf, PhasePortrait *portrait);

/**
 * @brief Generate phase portrait trajectories for a 2D system
 *        by integrating from a grid of initial conditions.
 *
 * @param vf          2D vector field
 * @param x_range     [x_min, x_max]
 * @param y_range     [y_min, y_max]
 * @param n_samples  Number of ICs per dimension
 * @param T           Integration time for each trajectory
 * @param portrait    Output: PhasePortrait
 * @return            0 on success
 */
int ode_generate_phase_portrait(const VectorField *vf,
                                  const double *x_range,
                                  const double *y_range,
                                  int n_samples, double T,
                                  PhasePortrait *portrait);

/**
 * @brief Free resources in a PhasePortrait.
 */
void ode_phase_portrait_free(PhasePortrait *portrait);

/* ──────────────────────── L2: Limit Cycle Detection ─────────────────── */

/**
 * @brief Detect limit cycles by computing Poincaré map on a section.
 *
 * A limit cycle manifests as a fixed point of the Poincaré map.
 * This function integrates from a Poincaré section and checks if
 * the trajectory returns to the section.
 *
 * @param vf            Vector field
 * @param section_tol   Tolerance for section crossing detection
 * @param T_max         Maximum integration time
 * @param cycle_period  Output: estimated period (0 if none found)
 * @param cycle_points  Output: trajectory on limit cycle
 * @param n_points      Output: number of points
 * @return              true if a limit cycle was detected
 */
bool ode_detect_limit_cycle(const VectorField *vf, double section_tol,
                              double T_max, double *cycle_period,
                              double *cycle_points, int *n_points);

/**
 * @brief Poincaré section crossing detection for a trajectory.
 *
 * Detects when a trajectory crosses the hyperplane n·(y - y₀) = 0
 * (used for period mapping, impact detection).
 *
 * @param vf       Vector field
 * @param y0       Initial state
 * @param n_vec    Normal vector of section plane
 * @param y_sec    Point on section plane
 * @param T_max    Maximum integration time
 * @param y_cross  Output: crossing state
 * @param t_cross  Output: crossing time
 * @return         true if crossing detected
 */
bool ode_poincare_section(const VectorField *vf, const double *y0,
                            const double *n_vec, const double *y_sec,
                            double T_max, double *y_cross, double *t_cross);

/* ──────────────────────── L6: Canonical Nonlinear Systems ────────────── */

/**
 * @brief Right-hand side for the harmonic oscillator:
 *        ẍ + (c/m)·ẋ + (k/m)·x = F(t)/m
 *
 * Converted to first-order system: y₁ = x, y₂ = ẋ
 *   ẏ₁ = y₂
 *   ẏ₂ = -(c/m)·y₂ - (k/m)·y₁ + F(t)/m
 *
 * @param t      Time
 * @param y      [x, ẋ] (length 2)
 * @param f_out  [ẋ, ẍ] (length 2)
 * @param dim    2
 * @param ctx    Cast to OscillatorParams*
 * @return       0
 */
int ode_rhs_harmonic_oscillator(double t, const double *y, double *f_out,
                                  int dim, void *ctx);

/**
 * @brief Right-hand side for the pendulum:
 *        θ̈ + (b/m)·θ̇ + (g/L)·sin(θ) = A·cos(ω_d·t)
 *
 * y₁ = θ, y₂ = θ̇
 *   ẏ₁ = y₂
 *   ẏ₂ = -(b/m)·y₂ - (g/L)·sin(y₁) + A·cos(ω_d·t)
 */
int ode_rhs_pendulum(double t, const double *y, double *f_out,
                       int dim, void *ctx);

/**
 * @brief Right-hand side for the Van der Pol oscillator:
 *        ẍ - μ·(1 - x²)·ẋ + x = A·cos(ω·t)
 *
 * y₁ = x, y₂ = ẋ
 *   ẏ₁ = y₂
 *   ẏ₂ = μ·(1 - y₁²)·y₂ - y₁ + A·cos(ω·t)
 */
int ode_rhs_vanderpol(double t, const double *y, double *f_out,
                        int dim, void *ctx);

/**
 * @brief Right-hand side for the Lotka-Volterra predator-prey system:
 *        ẋ = α·x - β·x·y
 *        ẏ = δ·x·y - γ·y
 */
int ode_rhs_lotka_volterra(double t, const double *y, double *f_out,
                             int dim, void *ctx);

/**
 * @brief Right-hand side for the Lorenz system:
 *        ẋ = σ·(y - x)
 *        ẏ = x·(ρ - z) - y
 *        ż = x·y - β·z
 *
 * Classic chaotic system for ρ > ρ_crit ≈ 24.74 (σ=10, β=8/3).
 */
int ode_rhs_lorenz(double t, const double *y, double *f_out,
                     int dim, void *ctx);

/* ──────────────────────── L2: Lyapunov Exponents ─────────────────────── */

/**
 * @brief Estimate the maximal Lyapunov exponent for a dynamical system
 *        using the two-particle (Benettin) method.
 *
 * λ_max > 0 indicates chaos. λ_max < 0 indicates stability.
 *
 * Reference: Benettin et al. "Lyapunov Characteristic Exponents for
 *            Smooth Dynamical Systems" (1980)
 *
 * @param vf          Vector field
 * @param y0          Initial condition
 * @param delta0      Initial separation (small perturbation)
 * @param T_total     Total integration time
 * @param n_steps     Number of renormalization steps
 * @param lambda_max  Output: largest Lyapunov exponent
 * @return            0 on success
 */
int ode_lyapunov_exponent(const VectorField *vf, const double *y0,
                            double delta0, double T_total, int n_steps,
                            double *lambda_max);

/**
 * @brief Compute the full Lyapunov spectrum using QR-based Benettin method.
 *
 * @param vf         Vector field
 * @param y0         Initial condition
 * @param T_total    Total integration time
 * @param n_steps    Number of reorthogonalization steps
 * @param spectrum   Output: Lyapunov spectrum (length dim, sorted descending)
 * @return           0 on success
 */
int ode_lyapunov_spectrum(const VectorField *vf, const double *y0,
                            double T_total, int n_steps,
                            double *spectrum);

/* ──────────────────────── L8: Bifurcation Analysis ───────────────────── */

/**
 * @brief Detect bifurcation in a parameter-dependent 1D system ẋ = f(x; μ)
 *        by checking sign changes of ∂f/∂x at equilibrium.
 *
 * @param f_mu_x      Function f(x; μ) for given μ
 * @param x_range     Range [x_min, x_max] to search
 * @param mu_min, mu_max Parameter range
 * @param n_mu        Number of μ samples
 * @param result      Output: bifurcation info
 * @return            0 on success
 */
int ode_detect_bifurcation_1d(double (*f_mu_x)(double x, double mu),
                                const double *x_range,
                                double mu_min, double mu_max, int n_mu,
                                BifurcationResult *result);

/**
 * @brief Check the Hopf bifurcation conditions for a 2D system:
 *        1. Equilibrium exists for parameter μ
 *        2. Jacobian has purely imaginary eigenvalues at μ = μ_crit
 *        3. Eigenvalues cross imaginary axis transversally
 *        4. First Lyapunov coefficient ≠ 0
 *
 * @param vf           Vector field parametrized by μ
 * @param mu_range     [mu_min, mu_max]
 * @param n_mu         Number of parameter samples
 * @param result       Output: bifurcation info
 * @return             true if Hopf bifurcation detected
 */
bool ode_detect_hopf_2d(const VectorField *vf, const double *mu_range,
                           int n_mu, BifurcationResult *result);

/* ──────────────────────── L4: Index Theory ──────────────────────────── */

/**
 * @brief Compute the Poincaré index of a 2D equilibrium.
 *
 * The index is the winding number of the vector field around a small
 * closed curve encircling the equilibrium. Indices: node/focus/center = +1,
 * saddle = -1.
 *
 * Reference: Strogatz (2015), Chap. 6.8
 *
 * @param vf          2D Vector field
 * @param eq_point    Equilibrium point (length 2)
 * @param radius      Radius of encircling curve
 * @param index_out   Output: Poincaré index
 * @return            0 on success
 */
int ode_poincare_index(const VectorField *vf, const double *eq_point,
                         double radius, int *index_out);

#ifdef __cplusplus
}
#endif

#endif /* ODE_NONLINEAR_H */
