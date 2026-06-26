/**
 * @file matrix.h
 * @brief Core matrix data structures and fundamental operations
 * @module mini-linear-algebra-control
 *
 * Knowledge coverage:
 *   L1: Matrix data type, row/column/diagonal access, determinant, trace
 *   L2: Linear transformations represented as matrices
 *   L3: Dense matrix storage, stride-based access, submatrix extraction
 *   L4: Cayley-Hamilton theorem verification via matrix polynomial
 *
 * Reference: Golub & Van Loan, Matrix Computations (4th Ed, 2013)
 *            Strang, Introduction to Linear Algebra (5th Ed, 2016)
 */

#ifndef MINI_LINALG_MATRIX_H
#define MINI_LINALG_MATRIX_H

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* =========================================================================
 * L1: Core Data Structure — Dense Matrix
 * ========================================================================= */

/**
 * @brief Dense matrix in row-major order.
 * data[i * cols + j] = element at row i, column j (0-indexed).
 */
typedef struct {
    double *data;    /**< Row-major flat array of elements */
    int     rows;    /**< Number of rows */
    int     cols;    /**< Number of columns */
    int     owner;   /**< 1 if this struct owns the data allocation, 0 if borrowed */
} Matrix;

/* L1: Construction and destruction */
Matrix* matrix_create(int rows, int cols);
Matrix* matrix_create_zero(int rows, int cols);
Matrix* matrix_create_identity(int n);
Matrix* matrix_create_ones(int rows, int cols);
Matrix* matrix_create_diag(const double *diag_vals, int n);
Matrix* matrix_create_copy(const Matrix *src);
Matrix* matrix_create_from_array(const double *vals, int rows, int cols);
void    matrix_free(Matrix *m);

/* L1: Element access */
double  matrix_get(const Matrix *m, int i, int j);
void    matrix_set(Matrix *m, int i, int j, double val);
double* matrix_row_ptr(Matrix *m, int i);
double* matrix_col_ptr(Matrix *m, int j);

/* L1: Dimension queries */
int matrix_is_square(const Matrix *m);
int matrix_is_symmetric(const Matrix *m, double tol);
int matrix_is_diagonal(const Matrix *m, double tol);
int matrix_is_upper_triangular(const Matrix *m, double tol);
int matrix_is_lower_triangular(const Matrix *m, double tol);
int matrix_is_orthogonal(const Matrix *m, double tol);

/* =========================================================================
 * L1/L2: Fundamental matrix operations
 * ========================================================================= */

/* L2: Scalar operations */
void matrix_scale(Matrix *m, double alpha);
void matrix_add_scalar(Matrix *m, double alpha);
void matrix_negate(Matrix *m);

/* L2: Element-wise operations */
int matrix_add(const Matrix *A, const Matrix *B, Matrix *result);
int matrix_subtract(const Matrix *A, const Matrix *B, Matrix *result);
int matrix_multiply(const Matrix *A, const Matrix *B, Matrix *result);
int matrix_hadamard(const Matrix *A, const Matrix *B, Matrix *result);

/* L2: Transpose */
int matrix_transpose(const Matrix *A, Matrix *result);
void matrix_transpose_inplace_square(Matrix *A);

/* L1: Trace — sum of diagonal elements */
double matrix_trace(const Matrix *A);

/* L1: Determinant — recursive Laplace expansion (small matrices)
 * and LU-based for general case */
double matrix_determinant(const Matrix *A);

/* L1: Rank via QR decomposition */
int matrix_rank(const Matrix *A, double tol);

/* L1: Matrix power — A^k for square matrices */
int matrix_power(const Matrix *A, int k, Matrix *result);

/* L4: Matrix polynomial evaluation — p(A) = c0*I + c1*A + c2*A^2 + ...
 * Used for Cayley-Hamilton verification */
int matrix_polynomial(const Matrix *A, const double *coeffs, int degree, Matrix *result);

/* L2: Inverse via LU decomposition */
int matrix_inverse(const Matrix *A, Matrix *result);

/* L2: Pseudoinverse via SVD (Moore-Penrose) */
int matrix_pseudoinverse(const Matrix *A, Matrix *result, double tol);

/* L2: Kronecker product */
int matrix_kronecker(const Matrix *A, const Matrix *B, Matrix *result);

/* L2: Vectorization (stack columns) */
int matrix_vec(const Matrix *A, Matrix *result);

/* L2: Submatrix extraction */
Matrix* matrix_submatrix(const Matrix *A, int row_start, int col_start,
                          int num_rows, int num_cols);

/* L2: Block matrix operations */
int matrix_block_multiply(const Matrix *A11, const Matrix *A12,
                          const Matrix *A21, const Matrix *A22,
                          const Matrix *B11, const Matrix *B12,
                          const Matrix *B21, const Matrix *B22,
                          Matrix *C11, Matrix *C12,
                          Matrix *C21, Matrix *C22);

/* L2: Concatenation */
int matrix_horzcat(const Matrix *A, const Matrix *B, Matrix *result);
int matrix_vertcat(const Matrix *A, const Matrix *B, Matrix *result);

/* L2: Matrix comparison */
int matrix_equal(const Matrix *A, const Matrix *B, double tol);

/* L2: Print matrix for debugging */
void matrix_print(const Matrix *A, const char *name);

#endif /* MINI_LINALG_MATRIX_H */
