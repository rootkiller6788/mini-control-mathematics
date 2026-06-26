/**
 * @file vector.h
 * @brief Vector operations and linear algebra primitives
 * @module mini-linear-algebra-control
 *
 * Knowledge coverage:
 *   L1: Vector data type, dot product, cross product, norm, unit vector
 *   L2: Linear combinations, span, linear independence testing
 *   L3: Vector spaces, basis construction, orthogonalization
 *   L4: Gram-Schmidt orthogonalization (constructive proof of basis existence)
 *
 * Reference: Strang, Introduction to Linear Algebra (5th Ed, 2016)
 *            Golub & Van Loan, Matrix Computations (4th Ed, 2013)
 */

#ifndef MINI_LINALG_VECTOR_H
#define MINI_LINALG_VECTOR_H

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* =========================================================================
 * L1: Core Data Structure — Vector
 * ========================================================================= */

/**
 * @brief Dense vector of double-precision values.
 */
typedef struct {
    double *data;  /**< Array of elements */
    int     size;  /**< Number of elements */
    int     owner; /**< 1 if this struct owns data allocation */
} Vector;

/* L1: Construction and destruction */
Vector* vector_create(int size);
Vector* vector_create_zero(int size);
Vector* vector_create_ones(int size);
Vector* vector_create_copy(const Vector *src);
Vector* vector_create_from_array(const double *vals, int size);
Vector* vector_create_linspace(double start, double end, int n);
Vector* vector_create_range(int start, int end);
void    vector_free(Vector *v);

/* L1: Element access */
double  vector_get(const Vector *v, int i);
void    vector_set(Vector *v, int i, double val);

/* L1: Scalar operations */
void vector_scale(Vector *v, double alpha);
void vector_negate(Vector *v);
void vector_add_scalar(Vector *v, double alpha);

/* L1: Element-wise vector operations */
int  vector_add(const Vector *a, const Vector *b, Vector *result);
int  vector_subtract(const Vector *a, const Vector *b, Vector *result);
int  vector_hadamard(const Vector *a, const Vector *b, Vector *result);
int  vector_axpy(double alpha, const Vector *x, Vector *y);

/* =========================================================================
 * L1: Inner product space operations
 * ========================================================================= */

/* L1: Dot product (Euclidean inner product) */
double vector_dot(const Vector *a, const Vector *b);

/* L1: Cross product (3D only) */
int vector_cross(const Vector *a, const Vector *b, Vector *result);

/* L1: Vector norms */
double vector_norm_l1(const Vector *v);
double vector_norm_l2(const Vector *v);
double vector_norm_linf(const Vector *v);
double vector_norm_lp(const Vector *v, double p);

/* L1: Normalize to unit vector */
int vector_normalize(Vector *v);

/* L1: Distance between vectors */
double vector_distance(const Vector *a, const Vector *b);

/* L1: Angle between two vectors (radians) */
double vector_angle(const Vector *a, const Vector *b);

/* L1: Projection — proj of b onto a */
int vector_projection(const Vector *a, const Vector *b, Vector *result);

/* =========================================================================
 * L2: Linear algebra vector operations
 * ========================================================================= */

/* L2: Linear combination — result = sum(c[i] * vectors[i])
 * Implements the concept of linear combination in vector spaces */
int vector_linear_combination(const Vector **vectors, const double *coeffs,
                               int n, Vector *result);

/* L2: Check if a set of vectors is linearly independent
 * Uses Gaussian elimination on the matrix formed by vectors as columns */
int vector_is_linearly_independent(const Vector **vectors, int n, double tol);

/* L2: Find the span dimension of a set of vectors
 * Returns the rank of the matrix formed by these vectors */
int vector_span_dimension(const Vector **vectors, int n, double tol);

/* L2: Express a vector as a linear combination of basis vectors
 * Solves the linear system for coefficients */
int vector_coordinates(const Vector *v, const Vector **basis, int n,
                        double *coeffs, double tol);

/* =========================================================================
 * L4: Gram-Schmidt Orthogonalization
 * ========================================================================= */

/**
 * @brief Classical Gram-Schmidt orthogonalization.
 * Given a set of vectors {v1, ..., vn}, produces orthonormal set {q1, ..., qn}
 * such that span{q1,...,qk} = span{v1,...,vk} for all k.
 *
 * This is a constructive proof that every finite-dimensional inner product
 * space has an orthonormal basis (L4).
 *
 * @param vectors  Input vectors (will be modified in-place to orthonormal)
 * @param n        Number of vectors
 * @param tol      Tolerance for zero-norm detection
 * @return 0 on success, -1 if vectors are linearly dependent
 */
int vector_gram_schmidt(Vector **vectors, int n, double tol);

/**
 * @brief Modified Gram-Schmidt (numerically more stable).
 * Same mathematical result as classical GS but better numerical properties.
 * (L5: Engineering method for numerical stability)
 */
int vector_gram_schmidt_modified(Vector **vectors, int n, double tol);

/* L2: Vector comparison */
int vector_equal(const Vector *a, const Vector *b, double tol);

/* L2: Print for debugging */
void vector_print(const Vector *v, const char *name);

/* L2: Convert between Vector and Matrix (column vector) */
int  vector_to_matrix_col(const Vector *v, void *matrix);
int  vector_from_matrix_col(void *matrix, int col, Vector *v);
void vector_print_sparse(const Vector *v, const char *name);

#endif /* MINI_LINALG_VECTOR_H */
