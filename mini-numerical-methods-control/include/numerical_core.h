/**
 * numerical_core.h — Core Numerical Definitions for Control Engineering
 *
 * Covers: L1 Definitions, L2 Core Concepts, L3 Engineering Quantities
 *
 * Defines fundamental vector/matrix data structures, norms, and conditioning
 * metrics essential for all numerical methods in control theory.
 *
 * Reference: Golub & Van Loan, Matrix Computations, 4th ed. (2013)
 *            Trefethen & Bau, Numerical Linear Algebra (1997)
 *            Higham, Accuracy and Stability of Numerical Algorithms (2002)
 */

#ifndef NUMERICAL_CORE_H
#define NUMERICAL_CORE_H

#include <stddef.h>
#include <math.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Core Data Structures — Vector and Matrix types
 * ========================================================================== */

/** Vector — dynamic-length real vector (column vector by convention) */
typedef struct {
    double *data;
    size_t  n;
} Vector;

/** Matrix — dynamic-size real matrix, row-major storage */
typedef struct {
    double *data;
    size_t  rows;
    size_t  cols;
} Matrix;

/** SVDResult — Singular Value Decomposition output A = U Σ V^T */
typedef struct {
    Matrix *U;
    Matrix *Vt;
    double *S;
    size_t  rank;
} SVDResult;

/* ==========================================================================
 * L2: Core Concepts — Condition Number and Numerical Stability
 * ========================================================================== */

/** ConditionNumber — stores various condition number estimates */
typedef struct {
    double cond_1;
    double cond_2;
    double cond_inf;
    double rcond;
} ConditionNumber;

/* ==========================================================================
 * L3: Engineering Quantities — Machine Constants
 * ========================================================================== */

#define MNC_EPSILON      2.220446049250313e-16
#define MNC_SQRT_EPSILON 1.490116119384766e-08
#define MNC_TINY         2.225073858507201e-308

/* ==========================================================================
 * Vector Operations
 * ========================================================================== */

Vector* vector_create(size_t n);

Vector* vector_create_from(size_t n, const double *values);

void vector_free(Vector *v);

Vector* vector_copy(const Vector *v);

Vector* vector_add(const Vector *u, const Vector *v);

Vector* vector_sub(const Vector *u, const Vector *v);

Vector* vector_scale(double alpha, const Vector *v);

void vector_axpy(double alpha, const Vector *x, Vector *y);

double vector_dot(const Vector *u, const Vector *v);

double vector_norm_l1(const Vector *v);

double vector_norm_l2(const Vector *v);

double vector_norm_inf(const Vector *v);

Vector* vector_normalize(const Vector *v);

Vector* vector_cross(const Vector *u, const Vector *v);

size_t vector_max_element(const Vector *v);

int vector_equals(const Vector *a, const Vector *b, double tol);

/* ==========================================================================
 * Matrix Operations
 * ========================================================================== */

Matrix* matrix_create(size_t rows, size_t cols);

Matrix* matrix_create_zeros(size_t rows, size_t cols);

Matrix* matrix_create_identity(size_t n);

Matrix* matrix_create_diag(const Vector *diagonal);

Matrix* matrix_create_from(size_t rows, size_t cols, const double *values);

void matrix_free(Matrix *A);

Matrix* matrix_copy(const Matrix *A);

double matrix_get(const Matrix *A, size_t i, size_t j);

void matrix_set(Matrix *A, size_t i, size_t j, double value);

Matrix* matrix_add(const Matrix *A, const Matrix *B);

Matrix* matrix_sub(const Matrix *A, const Matrix *B);

Matrix* matrix_scale(double alpha, const Matrix *A);

Matrix* matrix_multiply(const Matrix *A, const Matrix *B);

Matrix* matrix_transpose(const Matrix *A);

double matrix_trace(const Matrix *A);

double matrix_frobenius_norm(const Matrix *A);

double matrix_norm_1(const Matrix *A);

double matrix_norm_inf(const Matrix *A);

double matrix_norm_2_estimate(const Matrix *A, size_t max_iter, double tol);

double matrix_determinant(const Matrix *A);

int matrix_is_symmetric(const Matrix *A, double tol);

int matrix_is_positive_definite(const Matrix *A);

int matrix_equals(const Matrix *A, const Matrix *B, double tol);

/* ==========================================================================
 * Condition Number and Error Estimation
 * ========================================================================== */

ConditionNumber* condition_number_estimate(const Matrix *A);

void condition_number_free(ConditionNumber *cn);

double residual_norm(const Matrix *A, const Vector *x, const Vector *b);

double forward_error_bound(const ConditionNumber *cn, double backward_error);

#ifdef __cplusplus
}
#endif

#endif /* NUMERICAL_CORE_H */
