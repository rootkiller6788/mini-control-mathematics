/**
 * eigen_methods.h — Eigenvalue Computation for Stability Analysis
 *
 * Covers: L3 Math Structures, L4 Conservation (spectral properties), L5 Methods
 *
 * Eigenvalues determine stability, natural frequencies, and time constants
 * of control systems. Solutions of ẋ = Ax behave as Σ cᵢ e^{λᵢt} vᵢ.
 *
 * Reference: Wilkinson, The Algebraic Eigenvalue Problem (1965)
 *            Parlett, The Symmetric Eigenvalue Problem (1998)
 *            Saad, Numerical Methods for Large Eigenvalue Problems (2011)
 */

#ifndef EIGEN_METHODS_H
#define EIGEN_METHODS_H

#include "numerical_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Eigenvalue Result Types
 * ========================================================================== */

/**
 * EigenResult — stores eigenvalues and (optionally) eigenvectors
 *
 * L1 Definition: λ is an eigenvalue of A if ∃ v ≠ 0 such that Av = λv.
 * For complex eigenvalues, real_parts and imag_parts hold Re(λ) and Im(λ).
 */
typedef struct {
    double *real_parts;   /**< Real components of eigenvalues */
    double *imag_parts;   /**< Imaginary components (0 for real eigenvalues) */
    size_t  n;            /**< Number of eigenvalues (= matrix dimension) */
    Matrix *eigenvectors; /**< Real-form eigenvectors (n×2k for k complex pairs), NULL if not computed */
    int     converged;    /**< 1 if QR iteration converged */
} EigenResult;

/**
 * SpectralRadius — spectral information summary
 *
 * L2 Concept: ρ(A) = maxᵢ |λᵢ|. For discrete-time systems, stability
 * requires ρ(A) < 1. For continuous-time, Re(λᵢ) < 0 for all i.
 */
typedef struct {
    double spectral_radius;     /**< ρ(A) = max |λ| */
    double spectral_abscissa;   /**< α(A) = max Re(λ) — continuous-time stability */
    double min_real_part;       /**< Minimum real part */
    double max_imag_part;       /**< Maximum |Im(λ)| — damping indicator */
    int    is_hurwitz;          /**< 1 if all Re(λ) < 0 (continuous-time stable) */
    int    is_schur;            /**< 1 if all |λ| < 1 (discrete-time stable) */
} SpectralRadius;

/* ==========================================================================
 * L5: Eigenvalue Algorithms
 * ========================================================================== */

/**
 * power_iteration — compute dominant eigenvalue/eigenvector via power method
 *
 * L5 Method: x_{k+1} = A x_k / ||A x_k||, converges to eigenvector of |λ_max|.
 * Rayleigh quotient: λ ≈ xᵀAx / xᵀx.
 * Complexity: O(k·n²) for k iterations.
 * Convergence: linear, rate |λ₂/λ₁|.
 */
int power_iteration(const Matrix *A, double *eigenvalue, Vector *eigenvector,
                    size_t max_iter, double tol);

/**
 * inverse_iteration — compute eigenvalue nearest to shift σ
 *
 * L5 Method: Solve (A - σI)x_{k+1} = x_k repeatedly.
 * Convergence: cubic with Rayleigh quotient update.
 * Used for computing specific eigenvalues (e.g., those near imaginary axis).
 */
int inverse_iteration(const Matrix *A, double shift, double *eigenvalue,
                      Vector *eigenvector, size_t max_iter, double tol);

/**
 * rayleigh_quotient_iteration — Rayleigh quotient iteration
 *
 * L5/L8 Method: Shift σ_k = x_kᵀ A x_k / x_kᵀ x_k, then inverse iteration
 * with this shift. Cubic convergence! The fastest method for a single
 * eigenvalue, but requires solving a linear system each iteration.
 */
int rayleigh_quotient_iteration(const Matrix *A, double *eigenvalue,
                                Vector *eigenvector, size_t max_iter, double tol);

/**
 * qr_algorithm_eigenvalues — compute all eigenvalues via QR algorithm
 *
 * L5 Method: A₀ = A, then for k = 0,1,...: A_k = Q_k R_k, A_{k+1} = R_k Q_k.
 * A_k → quasi-upper-triangular (real Schur form) with eigenvalues on diagonal.
 * Hessenberg reduction first, then QR with shifts (Francis algorithm).
 *
 * Complexity: O(n³) for all eigenvalues.
 * This is the standard method used in MATLAB (eig) and LAPACK (dgeev).
 */
EigenResult* qr_algorithm_eigenvalues(const Matrix *A);

/**
 * eigen_symmetric — compute eigenvalues of symmetric A via QR (faster)
 *
 * L5 Method: For symmetric A, tridiagonalize first (O(n³)), then
 * symmetric QR (implicit QL with Wilkinson shifts). O(n²) per eigenvalue,
 * total O(n³) but with smaller constant than general case.
 */
EigenResult* eigen_symmetric(const Matrix *A);

/* ==========================================================================
 * L3/L5: Spectral Analysis for Control
 * ========================================================================== */

/**
 * spectral_radius_compute — compute ρ(A) and related quantities
 *
 * L2 Concept: Spectral radius determines asymptotic behavior of x_{k+1} = A x_k
 * and bounds ||A^k|| ≤ C ρ(A)^k for large k.
 *
 * L6 Application: For discrete LTI ẋ = Ax, convergence rate of iterative
 * learning control depends on spectral radius of the error propagation matrix.
 */
SpectralRadius* spectral_radius_compute(const EigenResult *eig);

/**
 * spectral_radius_free — deallocate SpectralRadius
 */
void spectral_radius_free(SpectralRadius *sr);

/**
 * eigen_result_free — deallocate EigenResult
 */
void eigen_result_free(EigenResult *eig);

/**
 * eigen_is_stable_continuous — check if all eigenvalues have negative real parts
 *
 * L6 Problem: Routh-Hurwitz criterion is algebraic; numerical verification
 * confirms via computed eigenvalues.
 */
int eigen_is_stable_continuous(const EigenResult *eig, double tol);

/**
 * eigen_is_stable_discrete — check if all eigenvalues have |λ| < 1
 *
 * L6 Problem: Schur stability for discrete systems.
 */
int eigen_is_stable_discrete(const EigenResult *eig, double tol);

/**
 * eigen_damping_ratios — compute damping ratios ζᵢ = -Re(λᵢ)/|λᵢ|
 *
 * L3 Engineering Quantity: Damping ratio ζ characterizes the oscillatory
 * behavior. ζ < 0: unstable; 0 ≤ ζ < 1: underdamped; ζ = 1: critically
 * damped; ζ > 1: overdamped.
 */
Vector* eigen_damping_ratios(const EigenResult *eig);

/**
 * eigen_natural_frequencies — compute ωₙᵢ = |λᵢ| for each mode
 *
 * L3 Engineering Quantity: Natural frequencies in rad/s. Relevant for
 * vibration control and flexible structure control.
 */
Vector* eigen_natural_frequencies(const EigenResult *eig);

/* ==========================================================================
 * L8: Large-Scale Methods
 * ========================================================================== */

/**
 * arnoldi_iteration — Arnoldi iteration for Krylov subspace basis
 *
 * L8 Method: Builds orthonormal basis V_m of K_m(A, v₁) = span{v₁, Av₁, ..., A^{m-1}v₁}.
 * H_m = V_mᵀ A V_m is upper Hessenberg. Eigenvalues of H_m (Ritz values)
 * approximate extreme eigenvalues of A.
 *
 * Used in model reduction (moment matching) and large sparse eigenvalue problems.
 */
int arnoldi_iteration(const Matrix *A, const Vector *v1, size_t m,
                      Matrix **V_out, Matrix **H_out);

/**
 * arnoldi_eigenvalues — compute Ritz eigenvalues from Arnoldi decomposition
 *
 * L8 Method: Approximates eigenvalues of large A using much smaller H_m.
 * Good for finding eigenvalues near the boundary of the spectrum.
 */
EigenResult* arnoldi_eigenvalues(const Matrix *A, const Vector *v1,
                                 size_t m, size_t n_eigen);

#ifdef __cplusplus
}
#endif

#endif /* EIGEN_METHODS_H */
