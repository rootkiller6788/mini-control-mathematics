/**
 * @file solvers.h
 * @brief Linear system solvers and least-squares methods
 * @module mini-linear-algebra-control
 *
 * Knowledge coverage:
 *   L1: Linear system definition, solution existence and uniqueness
 *   L2: Forward/backward substitution for triangular systems
 *   L5: Gaussian elimination, least-squares methods (normal equations, QR, SVD)
 *   L6: Underdetermined/overdetermined systems, minimum-norm solutions
 *
 * Reference: Golub & Van Loan, Matrix Computations (4th Ed, 2013)
 *            Bjorck, Numerical Methods for Least Squares Problems (1996)
 */

#ifndef MINI_LINALG_SOLVERS_H
#define MINI_LINALG_SOLVERS_H

#include "matrix.h"
#include "vector.h"

/* =========================================================================
 * L2: Triangular System Solvers
 * ========================================================================= */

/**
 * @brief Forward substitution: solve Lx = b where L is lower triangular.
 * Complexity: O(n^2).
 */
int solve_forward_substitution(const Matrix *L, const Vector *b, Vector *x);

/**
 * @brief Backward substitution: solve Ux = b where U is upper triangular.
 * Complexity: O(n^2).
 */
int solve_backward_substitution(const Matrix *U, const Vector *b, Vector *x);

/* =========================================================================
 * L5: Gaussian Elimination
 * ========================================================================= */

/**
 * @brief Solve Ax = b via Gaussian elimination with partial pivoting.
 * A is modified in-place. Returns solution in x.
 * Complexity: O(n^3) for elimination + O(n^2) for back-substitution.
 */
int solve_gaussian_elimination(Matrix *A, Vector *b, Vector *x);

/**
 * @brief Solve Ax = b via LU decomposition.
 * Factorizes A once, solves for each RHS in O(n^2).
 * For multiple RHS: use lu_solve() from decompositions.h.
 */
int solve_via_lu(const Matrix *A, const Vector *b, Vector *x);

/**
 * @brief Solve Ax = b via Cholesky (A must be SPD).
 */
int solve_via_cholesky(const Matrix *A, const Vector *b, Vector *x);

/**
 * @brief Solve multiple RHS: AX = B where B has m columns.
 * Uses LU decomposition. Complexity: O(n^3 + mn^2).
 */
int solve_multiple_rhs(const Matrix *A, const Matrix *B, Matrix *X);

/* =========================================================================
 * L5/L6: Least-Squares Problems
 * ========================================================================= */

/**
 * @brief Solve least-squares min ||Ax - b||_2 via normal equations.
 * Solves A^T A x = A^T b. Fast but potentially ill-conditioned
 * (squares the condition number).
 *
 * (L5: Normal equations method - simplest least-squares algorithm)
 */
int least_squares_normal_equations(const Matrix *A, const Vector *b, Vector *x);

/**
 * @brief Solve least-squares via QR decomposition.
 * Factor A = QR, solve Rx = Q^T b. Numerically more stable than normal eqs.
 * Complexity: O(mn^2) for QR + O(n^2) for solve.
 *
 * (L5: QR method - recommended for general least squares)
 */
int least_squares_qr(const Matrix *A, const Vector *b, Vector *x);

/**
 * @brief Solve least-squares via SVD.
 * Most robust method, handles rank-deficient A correctly.
 * x = V * S^+ * U^T * b where S^+ is pseudoinverse of singular values.
 *
 * (L5: SVD method - gold standard for ill-conditioned least squares)
 */
int least_squares_svd(const Matrix *A, const Vector *b, Vector *x, double tol);

/**
 * @brief Ridge regression (Tikhonov regularization):
 * min ||Ax - b||_2^2 + lambda^2 ||x||_2^2.
 * Equivalent to augmented system or normal equations with regularization.
 * (L8: Regularization for ill-posed problems)
 */
int least_squares_ridge(const Matrix *A, const Vector *b, double lambda,
                         Vector *x);

/**
 * @brief Total least squares (TLS): min ||[Delta_A  Delta_b]||_F
 * subject to (A + Delta_A)x = b + Delta_b.
 * Uses SVD of augmented matrix [A b].
 * (L8: Error-in-variables model, more general than ordinary LS)
 */
int least_squares_total(const Matrix *A, const Vector *b, Vector *x);

/**
 * @brief Minimum-norm solution for underdetermined systems (m < n).
 * min ||x||_2 subject to Ax = b. Uses Lagrange multipliers or SVD.
 * (L6: Canonical underdetermined problem)
 */
int solve_minimum_norm(const Matrix *A, const Vector *b, Vector *x);

/**
 * @brief Non-negative least squares (NNLS): min ||Ax - b||_2 subject to x >= 0.
 * Lawson-Hanson active set algorithm.
 * (L7: Application - constrained estimation in control)
 */
int least_squares_nonnegative(const Matrix *A, const Vector *b, Vector *x,
                               double tol, int max_iter);

/* =========================================================================
 * L7: Iterative Solvers
 * ========================================================================= */

/**
 * @brief Conjugate Gradient method for SPD systems: Ax = b.
 * Iterative method, complexity O(k n^2) for k iterations.
 * (L7: Large sparse systems, model predictive control)
 */
int solve_conjugate_gradient(const Matrix *A, const Vector *b, Vector *x,
                              double tol, int max_iter);

/**
 * @brief Preconditioned Conjugate Gradient.
 * Uses Jacobi (diagonal) preconditioner M = diag(A).
 * (L7: Accelerated CG for ill-conditioned SPD systems)
 */
int solve_cg_jacobi_preconditioned(const Matrix *A, const Vector *b,
                                    Vector *x, double tol, int max_iter);

/**
 * @brief GMRES for general (non-symmetric) systems.
 * Generalized Minimal Residual method.
 * (L7: Large sparse non-symmetric systems, discretized PDEs in control)
 */
int solve_gmres(const Matrix *A, const Vector *b, Vector *x,
                int restart, double tol, int max_iter);

#endif /* MINI_LINALG_SOLVERS_H */