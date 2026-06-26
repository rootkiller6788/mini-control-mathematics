/**
 * @file eigen.h
 * @brief Eigenvalue and eigenvector computation methods
 * @module mini-linear-algebra-control
 *
 * Knowledge coverage:
 *   L1: Eigenvalues, eigenvectors, characteristic polynomial, spectrum
 *   L2: Diagonalization, Jordan form, invariant subspaces
 *   L3: Spectral radius, spectral abscissa, eigenvalue sensitivity
 *   L4: Spectral theorem for symmetric matrices, Perron-Frobenius theorem
 *   L5: Power iteration, inverse iteration, Rayleigh quotient iteration, QR algorithm
 *   L6: Stability analysis via eigenvalues, modal decomposition
 *
 * Reference: Golub & Van Loan, Matrix Computations (4th Ed, 2013)
 *            Wilkinson, The Algebraic Eigenvalue Problem (1965)
 *            Parlett, The Symmetric Eigenvalue Problem (1998)
 */

#ifndef MINI_LINALG_EIGEN_H
#define MINI_LINALG_EIGEN_H

#include "matrix.h"
#include "vector.h"

/* =========================================================================
 * L1: Characteristic Polynomial
 * ========================================================================= */

/**
 * @brief Compute characteristic polynomial coefficients: det(lambda*I - A).
 * For nxn matrix, returns n+1 coefficients of polynomial:
 * p(lambda) = c[0] + c[1]*lambda + ... + c[n]*lambda^n
 * where c[n] = 1 (monic).
 *
 * Algorithm: Faddeev-Leverrier method. Complexity: O(n^4).
 * (L1: The characteristic polynomial is a fundamental definition of eigenvalues)
 */
int charpoly_faddeev_leverrier(const Matrix *A, double *coeffs);

/**
 * @brief Evaluate characteristic polynomial at a point lambda using
 * matrix determinant: det(lambda*I - A).
 */
double charpoly_evaluate(const Matrix *A, double lambda);

/* =========================================================================
 * L5: Power Iteration Methods
 * ========================================================================= */

/**
 * @brief Power iteration: find dominant eigenvalue and eigenvector.
 * v_{k+1} = A v_k / ||A v_k||. Converges to eigenvector of largest |lambda|.
 *
 * Theorem (L4): If A has a unique dominant eigenvalue, power iteration
 * converges linearly with rate |lambda_2 / lambda_1|.
 *
 * (L5: Simplest iterative eigenvalue algorithm)
 */
int eigen_power_iteration(const Matrix *A, Vector *v, double *lambda,
                           double tol, int max_iter);

/**
 * @brief Inverse power iteration: find eigenvalue closest to shift mu.
 * Uses (A - mu*I)^{-1}. Converges to eigenvalue nearest mu.
 * (L5: Refine eigenvalue estimates, compute specific eigenpairs)
 */
int eigen_inverse_iteration(const Matrix *A, double mu, Vector *v,
                             double *lambda, double tol, int max_iter);

/**
 * @brief Rayleigh quotient iteration: cubic convergence for symmetric matrices.
 * Combines inverse iteration with Rayleigh quotient shift updates.
 * (L5: Fastest local convergence for symmetric eigenproblems)
 */
int eigen_rayleigh_quotient_iteration(const Matrix *A, Vector *v,
                                       double *lambda, double tol, int max_iter);

/**
 * @brief Subspace iteration (simultaneous iteration): find k dominant
 * eigenpairs simultaneously. Block version of power iteration.
 * (L5: Compute multiple eigenvalues)
 */
int eigen_subspace_iteration(const Matrix *A, Matrix *V, Vector *lambdas,
                              int k, double tol, int max_iter);

/* =========================================================================
 * L5: QR Algorithm for Full Eigenvalue Computation
 * ========================================================================= */

/**
 * @brief QR algorithm for all eigenvalues of a symmetric matrix.
 * Reduces to tridiagonal form, then applies implicit shifted QR.
 * Returns eigenvalues in descending order.
 * (L5: The standard algorithm for dense symmetric eigenvalue problems)
 */
int eigen_qr_symmetric(const Matrix *A, Vector *eigenvalues, double tol, int max_iter);

/**
 * @brief QR algorithm for eigenvalues of a general (non-symmetric) matrix.
 * Uses Hessenberg reduction + double-shift implicit QR.
 * Returns real and imaginary parts of eigenvalues.
 * (L5: The workhorse algorithm for general eigenvalue problems)
 */
int eigen_qr_general(const Matrix *A, Vector *real_parts,
                      Vector *imag_parts, double tol, int max_iter);

/* =========================================================================
 * L5: Divide and Conquer (Symmetric Tridiagonal)
 * ========================================================================= */

/**
 * @brief Divide-and-conquer algorithm for symmetric tridiagonal eigenproblems.
 * Often faster than QR for large matrices on modern architectures.
 * (L5: Advanced eigenvalue algorithm with O(n^3) worst-case but faster average)
 */
int eigen_divide_conquer(const Matrix *A, Vector *eigenvalues,
                          Matrix *eigenvectors, double tol);

/* =========================================================================
 * L2/L6: Spectral Analysis
 * ========================================================================= */

/**
 * @brief Spectral radius: rho(A) = max |lambda_i(A)|.
 * (L3: Engineering quantity - determines asymptotic behavior of discrete systems)
 */
double eigen_spectral_radius(const Matrix *A);

/**
 * @brief Spectral abscissa: alpha(A) = max Re(lambda_i(A)).
 * (L3: Determines stability of continuous-time systems: stable if alpha < 0)
 */
double eigen_spectral_abscissa(const Matrix *A);

/**
 * @brief Check if matrix is Hurwitz (all eigenvalues have negative real parts).
 * (L6: Stability criterion for continuous-time linear systems)
 */
int eigen_is_hurwitz(const Matrix *A, double tol);

/**
 * @brief Check if matrix is Schur (all eigenvalues inside unit circle).
 * (L6: Stability criterion for discrete-time linear systems)
 */
int eigen_is_schur(const Matrix *A, double tol);

/**
 * @brief Compute eigenvector from known eigenvalue.
 * Finds null space of (A - lambda*I) via SVD.
 */
int eigen_eigenvector(const Matrix *A, double lambda_real,
                       double lambda_imag, Vector *vr, Vector *vi);

/**
 * @brief Modal matrix: columns are eigenvectors of A.
 * V such that V^{-1} A V = diag(lambda_i) (when diagonalizable).
 */
int eigen_modal_matrix(const Matrix *A, Matrix *V);

/**
 * @brief Generalized eigenvalue problem: A x = lambda B x.
 * (L6: Arises in mechanical vibrations, structural dynamics)
 */
int eigen_generalized(const Matrix *A, const Matrix *B,
                       Vector *real_parts, Vector *imag_parts,
                       double tol, int max_iter);

/* =========================================================================
 * L4: Spectral Theorem for Symmetric Matrices
 * ========================================================================= */

/**
 * @brief Full eigendecomposition of symmetric matrix: A = Q Lambda Q^T.
 *
 * Theorem (L4 - Spectral): Every real symmetric matrix A can be diagonalized
 * by an orthogonal matrix Q: A = Q Lambda Q^T, where Lambda = diag(lambda_i)
 * and lambda_i are real.
 *
 * @param A         Input nxn symmetric matrix
 * @param Q         Output orthogonal matrix of eigenvectors
 * @param lambdas   Output vector of eigenvalues (sorted descending)
 * @return          0 on success
 */
int eigen_symmetric_decompose(const Matrix *A, Matrix *Q, Vector *lambdas);

/**
 * @brief Compute matrix sign function: sign(A) from eigendecomposition.
 * sign(A) = Q * diag(sign(lambda_i)) * Q^T.
 * Used in solving Lyapunov and Riccati equations (L8).
 */
int matrix_sign_function(const Matrix *A, Matrix *S, double tol);

/**
 * @brief Compute matrix square root: A^{1/2} for SPD A.
 * A^{1/2} = Q * diag(sqrt(lambda_i)) * Q^T.
 */
int matrix_sqrt(const Matrix *A, Matrix *S);

#endif /* MINI_LINALG_EIGEN_H */