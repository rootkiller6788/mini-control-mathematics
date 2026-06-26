/**
 * @file    ode_linear.h
 * @brief   Linear ODE systems: fundamental matrix, state transition matrix,
 *          eigenvalue analysis, matrix exponential.
 *
 * L1-L4: Linear system definitions, superposition, modal decomposition,
 *        Wronskian, Liouville's formula, Floquet theory.
 *
 * Reference: Chen "Linear System Theory and Design" (1999)
 *            Hirsch & Smale "Differential Equations, Dynamical Systems" (2013)
 *
 * Course Mapping:
 *   MIT 6.241 – State transition matrix, matrix exponential
 *   Stanford ENGR 207B – Modal decomposition
 *   Berkeley ME 232 – Linear system analysis
 *   ETH 151-0555 – Linear control systems
 *   Tsinghua 现代控制理论 – State space
 */

#ifndef ODE_LINEAR_H
#define ODE_LINEAR_H

#include "ode_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ──────────────────────── L1: Linear System Creation ────────────────── */

/**
 * @brief Create a linear ODE system: ẏ = A y + b(t).
 *
 * @param A             System matrix (dim×dim row-major, copied)
 * @param b             Forcing vector (may be NULL for homogeneous)
 * @param dim           Dimension
 * @param time_varying  Is A time-dependent?
 * @return              LinearODESystem (caller owns A and b arrays)
 */
LinearODESystem ode_linear_create(const double *A, const double *b,
                                    int dim, bool time_varying);

/**
 * @brief Free resources associated with a linear ODE system.
 */
void ode_linear_destroy(LinearODESystem *sys);

/* ──────────────────────── L4: Matrix Exponential ─────────────────────── */

/**
 * @brief Compute the matrix exponential exp(A·t) using scaling-and-squaring
 *        with Padé approximation (Higham 2009).
 *
 * The state transition matrix for ẏ = A y is Φ(t) = exp(A·t).
 *
 * Reference: Higham "The Scaling and Squaring Method for the Matrix Exponential
 *            Revisited" (2009) – algorithm for expm.
 *
 * Complexity: O(dim³)
 *
 * @param A        Input matrix (dim×dim)
 * @param dim      Dimension
 * @param t        Time scaling factor
 * @param expAt    Output: exp(A·t) (dim×dim, row-major)
 * @return         0 on success
 */
int ode_matrix_exp(const double *A, int dim, double t, double *expAt);

/**
 * @brief Compute the matrix exponential via eigenvalue decomposition
 *        when A is diagonalizable. Fallback if not diagonalizable.
 */
int ode_matrix_exp_eigen(const double *A, int dim, double t, double *expAt);

/* ──────────────────────── L4: Modal Decomposition ─────────────────────── */

/**
 * @brief Compute eigenvalues and eigenvectors of a real matrix
 *        using the QR algorithm with double shifts.
 *
 * Critical for stability analysis: eigenvalues determine
 * solution behavior of linear ODE systems.
 *
 * @param A          Input matrix (dim×dim, row-major)
 * @param dim        Dimension
 * @param real_parts Output: real parts of eigenvalues (length dim)
 * @param imag_parts Output: imaginary parts of eigenvalues (length dim)
 * @param vecs       Output: eigenvectors as columns (dim×dim, may be NULL)
 * @return           0 on success
 */
int ode_eigen_decompose(const double *A, int dim,
                          double *real_parts, double *imag_parts,
                          double *vecs);

/**
 * @brief QR decomposition of a matrix (Householder reflections).
 *
 * Used internally for eigenvalue computation.
 */
int ode_qr_decompose(double *A, int rows, int cols, double *Q, double *R);

/**
 * @brief Solve 2×2 eigenvalue problem analytically for
 *        stability classification of planar systems.
 *
 * For A = [[a, b], [c, d]]:
 *   λ₁,₂ = (tr(A) ± √(tr(A)² - 4·det(A))) / 2
 *
 * Classification maps directly to phase portrait types (L2).
 *
 * @param a,b,c,d       Matrix entries
 * @param lambda1_real  Output: real part of λ₁
 * @param lambda1_imag  Output: imag part of λ₁
 * @param lambda2_real  Output: real part of λ₂
 * @param lambda2_imag  Output: imag part of λ₂
 * @return              EQUILIBRIUM classification
 */
EquilibriaClassification ode_classify_equilibrium_2d(
    double a, double b, double c, double d,
    double *lambda1_real, double *lambda1_imag,
    double *lambda2_real, double *lambda2_imag);

/* ──────────────────────── L4: Wronskian & Liouville ──────────────────── */

/**
 * @brief Compute the Wronskian W(t) = det(Φ(t)) for a linear system.
 *
 * Liouville's formula: W(t) = W(0)·exp(∫₀ᵗ tr(A(s)) ds).
 *
 * For constant A: W(t) = W(0)·exp(tr(A)·t).
 *
 * Reference: Coddington & Levinson (1955), Chap. 3
 *
 * @param sys     Linear ODE system
 * @param t       Time at which to evaluate
 * @param info    Output: WronskianInfo
 * @return        0 on success
 */
int ode_wronskian_compute(const LinearODESystem *sys, double t,
                           WronskianInfo *info);

/**
 * @brief Verify Liouville's formula for a constant-coefficient system.
 *
 * @param sys   Linear system with constant A
 * @param W0    Initial Wronskian (at t=0)
 * @param t     Time
 * @return      Computed W(t) via Liouville's formula
 */
double ode_liouville_formula(const LinearODESystem *sys, double W0, double t);

/* ──────────────────────── L4: Floquet Theory (Periodic) ───────────────── */

/**
 * @brief Compute the monodromy matrix Φ(T) for a periodic linear system
 *        ẏ = A(t)y where A(t+T) = A(t).
 *
 * The stability of periodic systems is governed by the Floquet multipliers
 * (eigenvalues of the monodromy matrix). If all |μ| < 1, the system is
 * asymptotically stable.
 *
 * Reference: Yakubovich & Starzhinskii "Linear Differential Equations
 *            with Periodic Coefficients" (1975)
 *
 * @param ivp           IVP with periodic RHS (period T)
 * @param T             Period
 * @param n_substeps    Subdivisions per period for integration
 * @param monodromy     Output: monodromy matrix M (dim×dim)
 * @param multipliers   Output: Floquet multipliers |μᵢ| (length dim)
 * @return              StabilityType from Floquet analysis
 */
StabilityType ode_floquet_analysis(const ODEIVP *ivp, double T,
                                     int n_substeps, double *monodromy,
                                     double *multipliers);

/* ──────────────────────── L2: Superposition Principle ────────────────── */

/**
 * @brief Verify the superposition principle for a linear homogeneous ODE:
 *        If y₁(t) and y₂(t) are solutions of ẏ = A(t)y, then
 *        c₁·y₁(t) + c₂·y₂(t) is also a solution for any c₁, c₂.
 *
 * This function constructs two solutions numerically and verifies that
 * the linear combination satisfies the ODE to within tolerance.
 *
 * @param sys     Linear system
 * @param r      Sampling radius for initial conditions
 * @param n_tests Number of random (c₁, c₂) pairs to test
 * @param max_err Output: maximum residual ||ẏ_comb - A·y_comb||
 * @return        true if superposition holds within tolerance
 */
bool ode_verify_superposition(const LinearODESystem *sys, double r,
                                int n_tests, double *max_err);

/* ──────────────────────── L2: Eigenvalue Stability ───────────────────── */

/**
 * @brief Determine stability type of an equilibrium at origin for ẏ = A y.
 *
 * All Re(λ) < 0  → asymptotically stable
 * Some Re(λ) = 0, none > 0 → marginally stable
 * Any Re(λ) > 0  → unstable
 *
 * @param A       System matrix (dim×dim)
 * @param dim     Dimension
 * @return        StabilityType
 */
StabilityType ode_stability_from_eigenvalues(const double *A, int dim);

/**
 * @brief Compute the dominant eigenvalue (largest real part) and the
 *        associated time constant τ = -1/Re(λ_dom) for stable systems.
 *
 * @param A         System matrix
 * @param dim       Dimension
 * @param dom_real  Output: real part of dominant eigenvalue
 * @param tau       Output: dominant time constant (-1/dom_real)
 * @return          StabilityType
 */
StabilityType ode_dominant_eigenvalue(const double *A, int dim,
                                        double *dom_real, double *tau);

/* ──────────────────────── L3: State Transition Matrix ────────────────── */

/**
 * @brief Compute the state transition matrix Φ(t, t₀) for a linear ODE.
 *
 * Properties:
 *   Φ(t₀, t₀) = I
 *   ∂Φ(t, t₀)/∂t = A(t)·Φ(t, t₀)
 *   x(t) = Φ(t, t₀)·x(t₀) + ∫_{t₀}^{t} Φ(t, τ)·b(τ) dτ
 *
 * For constant A: Φ(t, t₀) = exp(A·(t-t₀)).
 *
 * @param A       System matrix (constant or function value at t)
 * @param dim     Dimension
 * @param t       Current time
 * @param t0      Initial time
 * @param Phi     Output: transition matrix (dim×dim)
 * @return        0 on success
 */
int ode_state_transition(const double *A, int dim, double t, double t0,
                           double *Phi);

/* ──────────────────────── Sturm Comparison ──────────────────────────── */

/**
 * @brief Sturm separation theorem: Between any two consecutive zeros of
 *        one nontrivial solution of y'' + p(x)y = 0, there is exactly
 *        one zero of any other linearly independent solution.
 *
 * This function counts zeros of solutions for comparison.
 *
 * @param p_func    Function p(x) (caller provides values)
 * @param sol_a     First solution values (length n_pts)
 * @param sol_b     Second solution values (length n_pts)
 * @param x_vals    x evaluation points (length n_pts)
 * @param n_pts     Number of points
 * @param zeros_a   Output: number of zeros of sol_a
 * @param zeros_b   Output: number of zeros of sol_b
 * @return          0 on success
 */
int ode_sturm_separation(const double *p_func, const double *sol_a,
                           const double *sol_b, const double *x_vals,
                           int n_pts, int *zeros_a, int *zeros_b);

#ifdef __cplusplus
}
#endif

#endif /* ODE_LINEAR_H */
