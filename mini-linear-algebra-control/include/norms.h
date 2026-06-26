/**
 * @file norms.h
 * @brief Matrix and vector norms, condition numbers, error analysis
 * @module mini-linear-algebra-control
 *
 * Knowledge coverage:
 *   L1: Definitions of vector norms (L1, L2, Linf, Lp) and matrix norms
 *       (Frobenius, induced L1, induced Linf, spectral)
 *   L2: Norm equivalence in finite dimensions, submultiplicativity
 *   L3: Condition number as sensitivity measure for linear systems
 *   L4: Perturbation bounds for linear systems
 *
 * Reference: Golub & Van Loan, Matrix Computations (4th Ed, 2013)
 *            Higham, Accuracy and Stability of Numerical Algorithms (2nd Ed, 2002)
 */

#ifndef MINI_LINALG_NORMS_H
#define MINI_LINALG_NORMS_H

#include "matrix.h"
#include "vector.h"

/* =========================================================================
 * L1: Vector Norms
 * ========================================================================= */

double vector_norm_l1(const Vector *v);
double vector_norm_l2(const Vector *v);
double vector_norm_linf(const Vector *v);
double vector_norm_lp(const Vector *v, double p);

/* =========================================================================
 * L1: Matrix Norms
 * ========================================================================= */

/**
 * @brief Frobenius norm: ||A||_F = sqrt(sum_{i,j} |a_ij|^2).
 * (L1: Entry-wise norm, not induced by any vector norm)
 */
double matrix_norm_frobenius(const Matrix *A);

/**
 * @brief Induced L1 norm (maximum absolute column sum):
 * ||A||_1 = max_{1<=j<=n} sum_{i=1}^m |a_ij|
 */
double matrix_norm_induced_l1(const Matrix *A);

/**
 * @brief Induced Linf norm (maximum absolute row sum):
 * ||A||_inf = max_{1<=i<=m} sum_{j=1}^n |a_ij|
 */
double matrix_norm_induced_linf(const Matrix *A);

/**
 * @brief Spectral norm: ||A||_2 = sigma_max(A).
 * Largest singular value. Computed via power iteration or from SVD.
 */
double matrix_norm_spectral(const Matrix *A);

/**
 * @brief Nuclear norm (trace norm): ||A||_* = sum sigma_i(A).
 * Sum of singular values. Used in low-rank matrix recovery (L8).
 */
double matrix_norm_nuclear(const Matrix *A);

/**
 * @brief L2,1 norm: ||A||_{2,1} = sum_j ||a_j||_2.
 * Sum of L2 norms of columns. Promotes column sparsity.
 * (L8: Structured sparsity in system identification)
 */
double matrix_norm_l21(const Matrix *A);

/* =========================================================================
 * L3: Condition Numbers
 * ========================================================================= */

/**
 * @brief Condition number in L2 norm: kappa_2(A) = ||A||_2 * ||A^{-1}||_2.
 * (L3: Engineering quantity measuring sensitivity of linear system solutions)
 */
double matrix_condition_l2(const Matrix *A);

/**
 * @brief Condition number in L1 norm: kappa_1(A) = ||A||_1 * ||A^{-1}||_1.
 */
double matrix_condition_l1(const Matrix *A);

/**
 * @brief Condition number in Linf norm.
 */
double matrix_condition_linf(const Matrix *A);

/**
 * @brief Condition number in Frobenius norm.
 */
double matrix_condition_frobenius(const Matrix *A);

/**
 * @brief Estimate L2 condition number via power iteration (no inverse needed).
 * Uses Hager-Higham estimator. (L5: Engineering method)
 */
double matrix_condition_l2_estimate(const Matrix *A);

/**
 * @brief Compute condition number for controllability matrix.
 * (L6: Measures numerical robustness of controllability analysis)
 */
double controllability_condition_number(const Matrix *A, const Matrix *B);

/**
 * @brief Condition number for Lyapunov equation solution.
 * (L6: Sensitivity of Lyapunov solution to perturbations)
 */
double lyapunov_condition_number(const Matrix *A, double tol);

/* =========================================================================
 * L3: Numerical Stability Measures
 * ========================================================================= */

/**
 * @brief Compute growth factor for LU decomposition.
 * rho = max_{i,j,k} |a_{ij}^{(k)}| / max_{i,j} |a_{ij}|.
 * (L3: Measures numerical stability of Gaussian elimination)
 */
double lu_growth_factor(const Matrix *A);

/**
 * @brief Orthogonality measure: ||Q^T Q - I||_F for a matrix Q.
 * (L3: Assesses quality of orthogonal transformation)
 */
double orthogonality_error(const Matrix *Q);

/**
 * @brief Residual norm for linear system: r = b - Ax, returns ||r||.
 * (L3: Backward error measure)
 */
double linear_system_residual(const Matrix *A, const Vector *x, const Vector *b);

/**
 * @brief Residual for Lyapunov equation: ||A X + X A^T + Q||.
 */
double lyapunov_residual(const Matrix *A, const Matrix *X, const Matrix *Q);

#endif /* MINI_LINALG_NORMS_H */