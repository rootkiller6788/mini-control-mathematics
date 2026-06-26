/**
 * @file vector.c
 * @brief Vector operations implementation
 * @module mini-linear-algebra-control
 *
 * Knowledge coverage:
 *   L1: Vector construction, dot/cross products, norms, angle, projection
 *   L2: Linear combination, linear independence testing, span dimension
 *   L4: Gram-Schmidt orthogonalization (constructive proof)
 *   L5: Modified Gram-Schmidt (numerically stable variant)
 *
 * Reference: Strang, Introduction to Linear Algebra (5th Ed, 2016)
 *            Golub & Van Loan, Matrix Computations (4th Ed, 2013)
 */

#include "vector.h"
#include "matrix.h"
#include "solvers.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

/* =========================================================================
 * L1: Construction and Destruction
 * ========================================================================= */

Vector* vector_create(int size) {
    if (size <= 0) return NULL;
    Vector *v = (Vector*)malloc(sizeof(Vector));
    if (!v) return NULL;
    v->data = (double*)calloc((size_t)size, sizeof(double));
    if (!v->data) { free(v); return NULL; }
    v->size = size;
    v->owner = 1;
    return v;
}

Vector* vector_create_zero(int size) {
    return vector_create(size);
}

Vector* vector_create_ones(int size) {
    Vector *v = vector_create(size);
    if (!v) return NULL;
    for (int i = 0; i < size; i++) {
        v->data[i] = 1.0;
    }
    return v;
}

Vector* vector_create_copy(const Vector *src) {
    if (!src) return NULL;
    Vector *v = vector_create(src->size);
    if (!v) return NULL;
    memcpy(v->data, src->data, (size_t)src->size * sizeof(double));
    return v;
}

Vector* vector_create_from_array(const double *vals, int size) {
    if (!vals || size <= 0) return NULL;
    Vector *v = vector_create(size);
    if (!v) return NULL;
    memcpy(v->data, vals, (size_t)size * sizeof(double));
    return v;
}

Vector* vector_create_linspace(double start, double end, int n) {
    if (n <= 1) return NULL;
    Vector *v = vector_create(n);
    if (!v) return NULL;
    double step = (end - start) / (n - 1);
    for (int i = 0; i < n; i++) {
        v->data[i] = start + i * step;
    }
    return v;
}

Vector* vector_create_range(int start, int end) {
    int n = end - start + 1;
    if (n <= 0) return NULL;
    Vector *v = vector_create(n);
    if (!v) return NULL;
    for (int i = 0; i < n; i++) {
        v->data[i] = (double)(start + i);
    }
    return v;
}

void vector_free(Vector *v) {
    if (!v) return;
    if (v->owner && v->data) {
        free(v->data);
    }
    free(v);
}

/* =========================================================================
 * L1: Element Access
 * ========================================================================= */

double vector_get(const Vector *v, int i) {
    assert(v && v->data);
    assert(i >= 0 && i < v->size);
    return v->data[i];
}

void vector_set(Vector *v, int i, double val) {
    assert(v && v->data);
    assert(i >= 0 && i < v->size);
    v->data[i] = val;
}

/* =========================================================================
 * L1: Scalar Operations
 * ========================================================================= */

void vector_scale(Vector *v, double alpha) {
    if (!v || !v->data) return;
    for (int i = 0; i < v->size; i++) {
        v->data[i] *= alpha;
    }
}

void vector_negate(Vector *v) {
    vector_scale(v, -1.0);
}

void vector_add_scalar(Vector *v, double alpha) {
    if (!v || !v->data) return;
    for (int i = 0; i < v->size; i++) {
        v->data[i] += alpha;
    }
}

/* =========================================================================
 * L1: Element-wise Vector Operations
 * ========================================================================= */

int vector_add(const Vector *a, const Vector *b, Vector *result) {
    if (!a || !b || !result) return -1;
    if (a->size != b->size || a->size != result->size) return -1;
    for (int i = 0; i < a->size; i++) {
        result->data[i] = a->data[i] + b->data[i];
    }
    return 0;
}

int vector_subtract(const Vector *a, const Vector *b, Vector *result) {
    if (!a || !b || !result) return -1;
    if (a->size != b->size || a->size != result->size) return -1;
    for (int i = 0; i < a->size; i++) {
        result->data[i] = a->data[i] - b->data[i];
    }
    return 0;
}

int vector_hadamard(const Vector *a, const Vector *b, Vector *result) {
    if (!a || !b || !result) return -1;
    if (a->size != b->size || a->size != result->size) return -1;
    for (int i = 0; i < a->size; i++) {
        result->data[i] = a->data[i] * b->data[i];
    }
    return 0;
}

int vector_axpy(double alpha, const Vector *x, Vector *y) {
    if (!x || !y || x->size != y->size) return -1;
    for (int i = 0; i < x->size; i++) {
        y->data[i] += alpha * x->data[i];
    }
    return 0;
}

/* =========================================================================
 * L1: Inner Product Space Operations
 * ========================================================================= */

double vector_dot(const Vector *a, const Vector *b) {
    if (!a || !b || a->size != b->size) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < a->size; i++) {
        sum += a->data[i] * b->data[i];
    }
    return sum;
}

int vector_cross(const Vector *a, const Vector *b, Vector *result) {
    if (!a || !b || !result) return -1;
    if (a->size != 3 || b->size != 3 || result->size != 3) return -1;
    result->data[0] = a->data[1] * b->data[2] - a->data[2] * b->data[1];
    result->data[1] = a->data[2] * b->data[0] - a->data[0] * b->data[2];
    result->data[2] = a->data[0] * b->data[1] - a->data[1] * b->data[0];
    return 0;
}

double vector_norm_l1(const Vector *v) {
    if (!v) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < v->size; i++) {
        sum += fabs(v->data[i]);
    }
    return sum;
}

double vector_norm_l2(const Vector *v) {
    if (!v) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < v->size; i++) {
        sum += v->data[i] * v->data[i];
    }
    return sqrt(sum);
}

double vector_norm_linf(const Vector *v) {
    if (!v) return 0.0;
    double max_val = 0.0;
    for (int i = 0; i < v->size; i++) {
        double abs_val = fabs(v->data[i]);
        if (abs_val > max_val) max_val = abs_val;
    }
    return max_val;
}

double vector_norm_lp(const Vector *v, double p) {
    if (!v || p < 1.0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < v->size; i++) {
        sum += pow(fabs(v->data[i]), p);
    }
    return pow(sum, 1.0 / p);
}

int vector_normalize(Vector *v) {
    if (!v) return -1;
    double norm = vector_norm_l2(v);
    if (norm < 1e-300) return -1;  /* Cannot normalize zero vector */
    vector_scale(v, 1.0 / norm);
    return 0;
}

double vector_distance(const Vector *a, const Vector *b) {
    if (!a || !b || a->size != b->size) return -1.0;
    double sum = 0.0;
    for (int i = 0; i < a->size; i++) {
        double diff = a->data[i] - b->data[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

double vector_angle(const Vector *a, const Vector *b) {
    if (!a || !b) return 0.0;
    double dot = vector_dot(a, b);
    double na = vector_norm_l2(a);
    double nb = vector_norm_l2(b);
    if (na < 1e-300 || nb < 1e-300) return 0.0;
    double cos_angle = dot / (na * nb);
    /* Clamp to [-1, 1] for numerical stability */
    if (cos_angle > 1.0) cos_angle = 1.0;
    if (cos_angle < -1.0) cos_angle = -1.0;
    return acos(cos_angle);
}

int vector_projection(const Vector *a, const Vector *b, Vector *result) {
    /* Project b onto a: proj_a(b) = (a.b / a.a) * a */
    if (!a || !b || !result || a->size != b->size || a->size != result->size) {
        return -1;
    }
    double dot_ab = vector_dot(a, b);
    double dot_aa = vector_dot(a, a);
    if (dot_aa < 1e-300) return -1;  /* Zero vector projection undefined */
    double scalar = dot_ab / dot_aa;
    for (int i = 0; i < a->size; i++) {
        result->data[i] = scalar * a->data[i];
    }
    return 0;
}

/* =========================================================================
 * L2: Linear Combination and Linear Independence
 * ========================================================================= */

int vector_linear_combination(const Vector **vectors, const double *coeffs,
                               int n, Vector *result) {
    if (!vectors || !coeffs || !result || n <= 0) return -1;
    /* Initialize result to zero */
    memset(result->data, 0, (size_t)result->size * sizeof(double));
    for (int k = 0; k < n; k++) {
        if (!vectors[k]) return -1;
        if (vectors[k]->size != result->size) return -1;
        if (fabs(coeffs[k]) < 1e-300) continue;
        for (int i = 0; i < result->size; i++) {
            result->data[i] += coeffs[k] * vectors[k]->data[i];
        }
    }
    return 0;
}

int vector_is_linearly_independent(const Vector **vectors, int n, double tol) {
    if (!vectors || n <= 0) return 0;
    if (n == 1) return (vector_norm_l2(vectors[0]) > tol) ? 1 : 0;
    if (tol <= 0) tol = 1e-10;

    /* Form matrix with vectors as columns, compute rank */
    int dim = vectors[0]->size;
    Matrix *M = matrix_create(dim, n);
    if (!M) return 0;

    for (int j = 0; j < n; j++) {
        if (vectors[j]->size != dim) { matrix_free(M); return 0; }
        for (int i = 0; i < dim; i++) {
            matrix_set(M, i, j, vector_get(vectors[j], i));
        }
    }

    int rank = matrix_rank(M, tol);
    matrix_free(M);
    return (rank == n) ? 1 : 0;
}

int vector_span_dimension(const Vector **vectors, int n, double tol) {
    if (!vectors || n <= 0) return 0;
    if (tol <= 0) tol = 1e-10;

    int dim = vectors[0]->size;
    Matrix *M = matrix_create(dim, n);
    if (!M) return 0;

    for (int j = 0; j < n; j++) {
        if (!vectors[j]) { matrix_free(M); return 0; }
        if (vectors[j]->size != dim) { matrix_free(M); return 0; }
        for (int i = 0; i < dim; i++) {
            matrix_set(M, i, j, vector_get(vectors[j], i));
        }
    }

    int rank = matrix_rank(M, tol);
    matrix_free(M);
    return rank;
}

int vector_coordinates(const Vector *v, const Vector **basis, int n,
                        double *coeffs, double tol) {
    if (!v || !basis || !coeffs || n <= 0) return -1;
    if (tol <= 0) tol = 1e-10;

    /* Set up linear system: V * c = v where V columns are basis vectors */
    int dim = v->size;
    Matrix *Vmat = matrix_create(dim, n);
    Matrix *vcol = matrix_create(dim, 1);
    if (!Vmat || !vcol) { matrix_free(Vmat); matrix_free(vcol); return -1; }

    for (int j = 0; j < n; j++) {
        if (!basis[j] || basis[j]->size != dim) {
            matrix_free(Vmat); matrix_free(vcol); return -1;
        }
        for (int i = 0; i < dim; i++) {
            matrix_set(Vmat, i, j, vector_get(basis[j], i));
        }
        matrix_set(vcol, j, 0, vector_get(v, j));
    }

    /* Solve via least squares if overdetermined */
    Vector *bvec = vector_create_from_array(v->data, dim);
    Vector *xvec = vector_create(n);
    if (!bvec || !xvec) {
        vector_free(bvec); vector_free(xvec);
        matrix_free(Vmat); matrix_free(vcol);
        return -1;
    }

    int ret = least_squares_qr(Vmat, bvec, xvec);
    if (ret == 0) {
        for (int i = 0; i < n; i++) {
            coeffs[i] = vector_get(xvec, i);
        }
    }

    vector_free(bvec); vector_free(xvec);
    matrix_free(Vmat); matrix_free(vcol);
    return ret;
}

/* =========================================================================
 * L4: Classical Gram-Schmidt Orthogonalization
 *
 * Theorem: Every finite-dimensional inner product space has an orthonormal
 * basis. Constructive proof via Gram-Schmidt process.
 *
 * Algorithm: For k = 1..n:
 *   q_k = v_k - sum_{j=1}^{k-1} proj_{q_j}(v_k)
 *   q_k = q_k / ||q_k||
 * ========================================================================= */

int vector_gram_schmidt(Vector **vectors, int n, double tol) {
    if (!vectors || n <= 0) return -1;
    if (tol <= 0) tol = 1e-10;

    for (int k = 0; k < n; k++) {
        if (!vectors[k]) return -1;
        int dim = vectors[k]->size;

        /* Subtract projections onto previous orthonormal vectors */
        for (int j = 0; j < k; j++) {
            double dot = vector_dot(vectors[k], vectors[j]);
            for (int i = 0; i < dim; i++) {
                vectors[k]->data[i] -= dot * vectors[j]->data[i];
            }
        }

        /* Normalize */
        double norm = vector_norm_l2(vectors[k]);
        if (norm < tol) return -1;  /* Linear dependence detected */
        vector_scale(vectors[k], 1.0 / norm);
    }
    return 0;
}

/* =========================================================================
 * L5: Modified Gram-Schmidt (numerically more stable)
 *
 * Instead of subtracting all previous projections at once,
 * we apply projections one at a time in an inner loop.
 * This reduces accumulation of rounding errors.
 * ========================================================================= */

int vector_gram_schmidt_modified(Vector **vectors, int n, double tol) {
    if (!vectors || n <= 0) return -1;
    if (tol <= 0) tol = 1e-10;

    for (int k = 0; k < n; k++) {
        if (!vectors[k]) return -1;
        int dim = vectors[k]->size;

        /* Modified GS: process sequentially */
        for (int j = 0; j < k; j++) {
            double dot = vector_dot(vectors[k], vectors[j]);
            for (int i = 0; i < dim; i++) {
                vectors[k]->data[i] -= dot * vectors[j]->data[i];
            }
        }

        double norm = vector_norm_l2(vectors[k]);
        if (norm < tol) return -1;
        vector_scale(vectors[k], 1.0 / norm);
    }
    return 0;
}

/* =========================================================================
 * L2: Vector Equality and Printing
 * ========================================================================= */

int vector_equal(const Vector *a, const Vector *b, double tol) {
    if (!a || !b) return 0;
    if (a->size != b->size) return 0;
    for (int i = 0; i < a->size; i++) {
        if (fabs(a->data[i] - b->data[i]) > tol) return 0;
    }
    return 1;
}

void vector_print(const Vector *v, const char *name) {
    if (!v) { printf("%s = NULL\n", name ? name : "vector"); return; }
    printf("%s (size=%d): [", name ? name : "vector", v->size);
    for (int i = 0; i < v->size; i++) {
        printf("%s%12.6f", (i > 0) ? ", " : "", v->data[i]);
    }
    printf(" ]\n");
}

void vector_print_sparse(const Vector *v, const char *name) {
    if (!v) { printf("%s = NULL\n", name ? name : "vector"); return; }
    printf("%s (size=%d, nonzeros): [", name ? name : "vector", v->size);
    int first = 1;
    for (int i = 0; i < v->size; i++) {
        if (fabs(v->data[i]) > 1e-300) {
            printf("%s[%d]=%.6f", first ? "" : ", ", i, v->data[i]);
            first = 0;
        }
    }
    printf(" ]\n");
}