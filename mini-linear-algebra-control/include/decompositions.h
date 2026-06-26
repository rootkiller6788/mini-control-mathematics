/**
 * @file decompositions.h
 * @brief Matrix factorizations and decompositions
 * @module mini-linear-algebra-control
 *
 * Knowledge coverage:
 *   L1: Definition of LU, QR, Cholesky, SVD, Schur decompositions
 *   L2: Concept of matrix factorization as product of structured matrices
 *   L3: Computational complexity of each decomposition
 *   L4: Existence theorems (LU for invertible, Cholesky for SPD, etc.)
 *   L5: Algorithms - partial pivoting, Householder reflections, Givens rotations
 *
 * Reference: Golub & Van Loan, Matrix Computations (4th Ed, 2013)
 *            Trefethen & Bau, Numerical Linear Algebra (1997)
 */

#ifndef MINI_LINALG_DECOMPOSITIONS_H
#define MINI_LINALG_DECOMPOSITIONS_H

#include "matrix.h"
#include "vector.h"

/* LU Decomposition */
int lu_decompose(Matrix *A, int *piv);
int lu_solve(const Matrix *LU, const int *piv, const Vector *b, Vector *x);
int lu_extract_l(const Matrix *LU, Matrix *L);
int lu_extract_u(const Matrix *LU, Matrix *U);
double lu_determinant(const Matrix *LU, const int *piv);
int lu_inverse(const Matrix *A, Matrix *Ainv);

/* QR Decomposition */
int qr_decompose_householder(Matrix *A, Matrix *Q);
int qr_extract_r(const Matrix *A, Matrix *R);
int qr_solve(const Matrix *Q, const Matrix *R, const Vector *b, Vector *x);
int qr_decompose_mgs(const Matrix *A, Matrix *Q, Matrix *R);
int qr_givens_rotate(Matrix *A, int i, int j, int pivot_row, double *c, double *s);

/* Cholesky Decomposition */
int cholesky_decompose(Matrix *A);
int cholesky_solve(const Matrix *L, const Vector *b, Vector *x);
int ldl_decompose(Matrix *A, Vector *D);

/* SVD */
int svd_jacobi(Matrix *A, Matrix *U, Vector *S, Matrix *Vt, double tol, int max_iter);
double svd_condition_number(const Vector *S);
int svd_effective_rank(const Vector *S, double tol);
int svd_truncate(const Matrix *U, const Vector *S, const Matrix *Vt, int k, Matrix *U_k, Vector *S_k, Matrix *Vt_k);

/* Schur Decomposition */
int schur_decompose(const Matrix *A, Matrix *Q, Matrix *T, double tol, int max_iter);
int schur_eigenvalues(const Matrix *T, Vector *real_parts, Vector *imag_parts);
int hessenberg_reduce(const Matrix *A, Matrix *H, Matrix *Q);

#endif /* MINI_LINALG_DECOMPOSITIONS_H */