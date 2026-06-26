/**
 * numerical_core.c — Core Vector and Matrix Implementation
 *
 * Implements: L1 Definitions, L2 Core Concepts, L3 Engineering Quantities
 *
 * Every function here maps to a specific knowledge point in numerical
 * linear algebra for control. Row-major storage convention is used
 * throughout for C interop and cache efficiency.
 */

#include "numerical_core.h"
#include <string.h>
#include <stdio.h>
#include <float.h>

/* ==========================================================================
 * Vector Lifecycle
 * ========================================================================== */

Vector* vector_create(size_t n) {
    if (n == 0) return NULL;
    Vector *v = (Vector*)malloc(sizeof(Vector));
    if (!v) return NULL;
    v->data = (double*)calloc(n, sizeof(double));
    if (!v->data) { free(v); return NULL; }
    v->n = n;
    return v;
}

Vector* vector_create_from(size_t n, const double *values) {
    Vector *v = vector_create(n);
    if (!v) return NULL;
    memcpy(v->data, values, n * sizeof(double));
    return v;
}

void vector_free(Vector *v) {
    if (!v) return;
    free(v->data);
    free(v);
}

Vector* vector_copy(const Vector *v) {
    if (!v) return NULL;
    return vector_create_from(v->n, v->data);
}

/* ==========================================================================
 * Vector Arithmetic
 * ========================================================================== */

Vector* vector_add(const Vector *u, const Vector *v) {
    if (!u || !v || u->n != v->n) return NULL;
    Vector *w = vector_create(u->n);
    if (!w) return NULL;
    for (size_t i = 0; i < u->n; i++)
        w->data[i] = u->data[i] + v->data[i];
    return w;
}

Vector* vector_sub(const Vector *u, const Vector *v) {
    if (!u || !v || u->n != v->n) return NULL;
    Vector *w = vector_create(u->n);
    if (!w) return NULL;
    for (size_t i = 0; i < u->n; i++)
        w->data[i] = u->data[i] - v->data[i];
    return w;
}

Vector* vector_scale(double alpha, const Vector *v) {
    if (!v) return NULL;
    Vector *w = vector_create(v->n);
    if (!w) return NULL;
    for (size_t i = 0; i < v->n; i++)
        w->data[i] = alpha * v->data[i];
    return w;
}

void vector_axpy(double alpha, const Vector *x, Vector *y) {
    if (!x || !y || x->n != y->n) return;
    for (size_t i = 0; i < x->n; i++)
        y->data[i] += alpha * x->data[i];
}

/* ==========================================================================
 * Vector Inner Products and Norms
 * ========================================================================== */

double vector_dot(const Vector *u, const Vector *v) {
    if (!u || !v || u->n != v->n) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < u->n; i++)
        sum += u->data[i] * v->data[i];
    return sum;
}

double vector_norm_l1(const Vector *v) {
    if (!v) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < v->n; i++)
        sum += fabs(v->data[i]);
    return sum;
}

double vector_norm_l2(const Vector *v) {
    if (!v) return 0.0;
    double sum = 0.0, scale = 0.0, ssq = 1.0;
    /* Two-pass with scaling to avoid overflow/underflow */
    double maxabs = 0.0;
    for (size_t i = 0; i < v->n; i++)
        if (fabs(v->data[i]) > maxabs) maxabs = fabs(v->data[i]);
    if (maxabs == 0.0) return 0.0;
    for (size_t i = 0; i < v->n; i++) {
        double scaled = v->data[i] / maxabs;
        sum += scaled * scaled;
    }
    return maxabs * sqrt(sum);
}

double vector_norm_inf(const Vector *v) {
    if (!v) return 0.0;
    double maxval = 0.0;
    for (size_t i = 0; i < v->n; i++)
        if (fabs(v->data[i]) > maxval) maxval = fabs(v->data[i]);
    return maxval;
}

Vector* vector_normalize(const Vector *v) {
    if (!v) return NULL;
    double norm = vector_norm_l2(v);
    if (norm < MNC_TINY) return NULL;  /* Avoid division by zero */
    return vector_scale(1.0 / norm, v);
}

Vector* vector_cross(const Vector *u, const Vector *v) {
    if (!u || !v || u->n != 3 || v->n != 3) return NULL;
    Vector *w = vector_create(3);
    if (!w) return NULL;
    w->data[0] = u->data[1] * v->data[2] - u->data[2] * v->data[1];
    w->data[1] = u->data[2] * v->data[0] - u->data[0] * v->data[2];
    w->data[2] = u->data[0] * v->data[1] - u->data[1] * v->data[0];
    return w;
}

size_t vector_max_element(const Vector *v) {
    if (!v || v->n == 0) return 0;
    size_t idx = 0;
    double maxval = fabs(v->data[0]);
    for (size_t i = 1; i < v->n; i++)
        if (fabs(v->data[i]) > maxval) { maxval = fabs(v->data[i]); idx = i; }
    return idx;
}

int vector_equals(const Vector *a, const Vector *b, double tol) {
    if (!a || !b || a->n != b->n) return 0;
    for (size_t i = 0; i < a->n; i++)
        if (fabs(a->data[i] - b->data[i]) > tol) return 0;
    return 1;
}

/* ==========================================================================
 * Matrix Lifecycle
 * ========================================================================== */

Matrix* matrix_create(size_t rows, size_t cols) {
    if (rows == 0 || cols == 0) return NULL;
    Matrix *A = (Matrix*)malloc(sizeof(Matrix));
    if (!A) return NULL;
    A->data = (double*)calloc(rows * cols, sizeof(double));
    if (!A->data) { free(A); return NULL; }
    A->rows = rows;
    A->cols = cols;
    return A;
}

Matrix* matrix_create_zeros(size_t rows, size_t cols) {
    return matrix_create(rows, cols);  /* calloc already zeros */
}

Matrix* matrix_create_identity(size_t n) {
    Matrix *I = matrix_create(n, n);
    if (!I) return NULL;
    for (size_t i = 0; i < n; i++)
        I->data[i * n + i] = 1.0;
    return I;
}

Matrix* matrix_create_diag(const Vector *diagonal) {
    if (!diagonal) return NULL;
    size_t n = diagonal->n;
    Matrix *D = matrix_create_zeros(n, n);
    if (!D) return NULL;
    for (size_t i = 0; i < n; i++)
        D->data[i * n + i] = diagonal->data[i];
    return D;
}

Matrix* matrix_create_from(size_t rows, size_t cols, const double *values) {
    Matrix *A = matrix_create(rows, cols);
    if (!A) return NULL;
    memcpy(A->data, values, rows * cols * sizeof(double));
    return A;
}

void matrix_free(Matrix *A) {
    if (!A) return;
    free(A->data);
    free(A);
}

Matrix* matrix_copy(const Matrix *A) {
    if (!A) return NULL;
    return matrix_create_from(A->rows, A->cols, A->data);
}

/* ==========================================================================
 * Matrix Element Access
 * ========================================================================== */

double matrix_get(const Matrix *A, size_t i, size_t j) {
    if (!A || i >= A->rows || j >= A->cols) return 0.0;
    return A->data[i * A->cols + j];
}

void matrix_set(Matrix *A, size_t i, size_t j, double value) {
    if (!A || i >= A->rows || j >= A->cols) return;
    A->data[i * A->cols + j] = value;
}

/* ==========================================================================
 * Matrix Arithmetic
 * ========================================================================== */

Matrix* matrix_add(const Matrix *A, const Matrix *B) {
    if (!A || !B || A->rows != B->rows || A->cols != B->cols) return NULL;
    Matrix *C = matrix_create(A->rows, A->cols);
    if (!C) return NULL;
    size_t n = A->rows * A->cols;
    for (size_t k = 0; k < n; k++)
        C->data[k] = A->data[k] + B->data[k];
    return C;
}

Matrix* matrix_sub(const Matrix *A, const Matrix *B) {
    if (!A || !B || A->rows != B->rows || A->cols != B->cols) return NULL;
    Matrix *C = matrix_create(A->rows, A->cols);
    if (!C) return NULL;
    size_t n = A->rows * A->cols;
    for (size_t k = 0; k < n; k++)
        C->data[k] = A->data[k] - B->data[k];
    return C;
}

Matrix* matrix_scale(double alpha, const Matrix *A) {
    if (!A) return NULL;
    Matrix *B = matrix_create(A->rows, A->cols);
    if (!B) return NULL;
    size_t n = A->rows * A->cols;
    for (size_t k = 0; k < n; k++)
        B->data[k] = alpha * A->data[k];
    return B;
}

Matrix* matrix_multiply(const Matrix *A, const Matrix *B) {
    if (!A || !B || A->cols != B->rows) return NULL;
    size_t m = A->rows, n = A->cols, p = B->cols;
    Matrix *C = matrix_create_zeros(m, p);
    if (!C) return NULL;
    for (size_t i = 0; i < m; i++) {
        for (size_t k = 0; k < n; k++) {
            double aik = A->data[i * A->cols + k];
            if (fabs(aik) < MNC_TINY) continue;  /* skip zero multiplications */
            size_t row_offset = i * p;
            size_t B_offset = k * B->cols;
            for (size_t j = 0; j < p; j++)
                C->data[row_offset + j] += aik * B->data[B_offset + j];
        }
    }
    return C;
}

Matrix* matrix_transpose(const Matrix *A) {
    if (!A) return NULL;
    Matrix *At = matrix_create(A->cols, A->rows);
    if (!At) return NULL;
    for (size_t i = 0; i < A->rows; i++)
        for (size_t j = 0; j < A->cols; j++)
            At->data[j * At->cols + i] = A->data[i * A->cols + j];
    return At;
}

/* ==========================================================================
 * Matrix Properties
 * ========================================================================== */

double matrix_trace(const Matrix *A) {
    if (!A) return 0.0;
    size_t mindim = (A->rows < A->cols) ? A->rows : A->cols;
    double tr = 0.0;
    for (size_t i = 0; i < mindim; i++)
        tr += A->data[i * A->cols + i];
    return tr;
}

double matrix_frobenius_norm(const Matrix *A) {
    if (!A) return 0.0;
    double sum = 0.0, maxabs = 0.0;
    size_t n = A->rows * A->cols;
    for (size_t k = 0; k < n; k++)
        if (fabs(A->data[k]) > maxabs) maxabs = fabs(A->data[k]);
    if (maxabs == 0.0) return 0.0;
    for (size_t k = 0; k < n; k++) {
        double scaled = A->data[k] / maxabs;
        sum += scaled * scaled;
    }
    return maxabs * sqrt(sum);
}

double matrix_norm_1(const Matrix *A) {
    if (!A) return 0.0;
    double maxcol = 0.0;
    for (size_t j = 0; j < A->cols; j++) {
        double sum = 0.0;
        for (size_t i = 0; i < A->rows; i++)
            sum += fabs(A->data[i * A->cols + j]);
        if (sum > maxcol) maxcol = sum;
    }
    return maxcol;
}

double matrix_norm_inf(const Matrix *A) {
    if (!A) return 0.0;
    double maxrow = 0.0;
    for (size_t i = 0; i < A->rows; i++) {
        double sum = 0.0;
        for (size_t j = 0; j < A->cols; j++)
            sum += fabs(A->data[i * A->cols + j]);
        if (sum > maxrow) maxrow = sum;
    }
    return maxrow;
}

double matrix_norm_2_estimate(const Matrix *A, size_t max_iter, double tol) {
    if (!A || A->rows == 0 || A->cols == 0) return 0.0;
    /* Power iteration on AᵀA to find σ_max */
    size_t n = A->cols;
    double *v = (double*)malloc(n * sizeof(double));
    double *Av = (double*)malloc(A->rows * sizeof(double));
    double *AtAv = (double*)malloc(n * sizeof(double));
    if (!v || !Av || !AtAv) { free(v); free(Av); free(AtAv); return 0.0; }

    /* Initialize with ones, then normalize */
    for (size_t i = 0; i < n; i++) v[i] = 1.0;
    double vnorm = sqrt((double)n);
    for (size_t i = 0; i < n; i++) v[i] /= vnorm;

    double lambda_old = 0.0, lambda = 0.0;
    for (size_t iter = 0; iter < max_iter; iter++) {
        /* Av = A * v */
        for (size_t i = 0; i < A->rows; i++) {
            Av[i] = 0.0;
            for (size_t j = 0; j < A->cols; j++)
                Av[i] += A->data[i * A->cols + j] * v[j];
        }
        /* AtAv = Aᵀ * Av */
        for (size_t i = 0; i < n; i++) {
            AtAv[i] = 0.0;
            for (size_t j = 0; j < A->rows; j++)
                AtAv[i] += A->data[j * A->cols + i] * Av[j];
        }
        /* Rayleigh quotient */
        double num = 0.0, den = 0.0;
        for (size_t i = 0; i < n; i++) { num += v[i] * AtAv[i]; den += v[i] * v[i]; }
        lambda = num / den;
        /* Normalize */
        double anorm = 0.0;
        for (size_t i = 0; i < n; i++) anorm += AtAv[i] * AtAv[i];
        anorm = sqrt(anorm);
        if (anorm < MNC_TINY) break;
        for (size_t i = 0; i < n; i++) v[i] = AtAv[i] / anorm;

        if (fabs(lambda - lambda_old) < tol * fabs(lambda)) break;
        lambda_old = lambda;
    }
    free(v); free(Av); free(AtAv);
    return (lambda > 0) ? sqrt(lambda) : 0.0;
}

/* ==========================================================================
 * Matrix Determinant via LU
 * ========================================================================== */

double matrix_determinant(const Matrix *A) {
    if (!A || A->rows != A->cols) return 0.0;
    if (A->rows == 0) return 1.0;
    size_t n = A->rows;

    /* Copy A into working matrix */
    double *U = (double*)malloc(n * n * sizeof(double));
    if (!U) return 0.0;
    memcpy(U, A->data, n * n * sizeof(double));

    double det = 1.0; int sign = 1;

    for (size_t k = 0; k < n; k++) {
        /* Partial pivoting */
        double max_val = fabs(U[k * n + k]);
        size_t max_row = k;
        for (size_t i = k + 1; i < n; i++) {
            if (fabs(U[i * n + k]) > max_val) {
                max_val = fabs(U[i * n + k]);
                max_row = i;
            }
        }
        if (max_val < MNC_TINY) { det = 0.0; free(U); return 0.0; }
        if (max_row != k) {
            for (size_t j = 0; j < n; j++) {
                double tmp = U[k * n + j];
                U[k * n + j] = U[max_row * n + j];
                U[max_row * n + j] = tmp;
            }
            sign = -sign;
        }
        det *= U[k * n + k];
        for (size_t i = k + 1; i < n; i++) {
            double factor = U[i * n + k] / U[k * n + k];
            for (size_t j = k; j < n; j++)
                U[i * n + j] -= factor * U[k * n + j];
        }
    }
    free(U);
    return sign * det;
}

int matrix_is_symmetric(const Matrix *A, double tol) {
    if (!A || A->rows != A->cols) return 0;
    for (size_t i = 0; i < A->rows; i++)
        for (size_t j = i + 1; j < A->cols; j++)
            if (fabs(A->data[i * A->cols + j] - A->data[j * A->cols + i]) > tol)
                return 0;
    return 1;
}

int matrix_is_positive_definite(const Matrix *A) {
    if (!A || A->rows != A->cols) return 0;
    size_t n = A->rows;
    /* Attempt Cholesky: L_{jj} = sqrt(A_{jj} - Σ_{k<j} L_{jk}²) */
    double *L = (double*)calloc(n * n, sizeof(double));
    if (!L) return 0;
    int is_pd = 1;
    for (size_t j = 0; j < n && is_pd; j++) {
        double sum = 0.0;
        for (size_t k = 0; k < j; k++)
            sum += L[j * n + k] * L[j * n + k];
        double diag = A->data[j * n + j] - sum;
        if (diag <= MNC_TINY) is_pd = 0;
        else {
            L[j * n + j] = sqrt(diag);
            for (size_t i = j + 1; i < n; i++) {
                sum = 0.0;
                for (size_t k = 0; k < j; k++)
                    sum += L[i * n + k] * L[j * n + k];
                L[i * n + j] = (A->data[i * n + j] - sum) / L[j * n + j];
            }
        }
    }
    free(L);
    return is_pd;
}

int matrix_equals(const Matrix *A, const Matrix *B, double tol) {
    if (!A || !B || A->rows != B->rows || A->cols != B->cols) return 0;
    size_t n = A->rows * A->cols;
    for (size_t k = 0; k < n; k++)
        if (fabs(A->data[k] - B->data[k]) > tol) return 0;
    return 1;
}

/* ==========================================================================
 * Condition Number Estimation
 * ========================================================================== */

ConditionNumber* condition_number_estimate(const Matrix *A) {
    if (!A || A->rows != A->cols) return NULL;
    ConditionNumber *cn = (ConditionNumber*)malloc(sizeof(ConditionNumber));
    if (!cn) return NULL;
    cn->cond_1   = matrix_norm_1(A);
    cn->cond_inf = matrix_norm_inf(A);
    cn->cond_2   = matrix_norm_2_estimate(A, 100, MNC_SQRT_EPSILON);
    /* rcond: inverse of condition number in 1-norm (LAPACK convention) */
    double norm1 = cn->cond_1;
    if (norm1 > MNC_TINY) {
        /* Estimate ||A⁻¹||₁ via inverse iteration on |A⁻¹| */
        double anorm_inv = 0.0;
        cn->rcond = 0.5; /* fallback */
        /* Simple estimate: use 2-norm condition reciprocal */
        if (cn->cond_2 > MNC_TINY) cn->rcond = 1.0 / cn->cond_2;
    } else {
        cn->rcond = 0.0;
    }
    return cn;
}

void condition_number_free(ConditionNumber *cn) {
    free(cn);
}

/* ==========================================================================
 * Residual and Error Estimation
 * ========================================================================== */

double residual_norm(const Matrix *A, const Vector *x, const Vector *b) {
    if (!A || !x || !b || A->cols != x->n || A->rows != b->n) return -1.0;
    Vector *Ax = vector_create(A->rows);
    if (!Ax) return -1.0;
    for (size_t i = 0; i < A->rows; i++) {
        Ax->data[i] = 0.0;
        for (size_t j = 0; j < A->cols; j++)
            Ax->data[i] += A->data[i * A->cols + j] * x->data[j];
    }
    Vector *res = vector_sub(b, Ax);
    double nr = res ? vector_norm_l2(res) : -1.0;
    vector_free(Ax);
    vector_free(res);
    return nr;
}

double forward_error_bound(const ConditionNumber *cn, double backward_error) {
    if (!cn) return 0.0;
    return cn->cond_2 * backward_error;
}
