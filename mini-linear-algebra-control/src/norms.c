/**
 * @file norms.c
 * @brief Matrix and vector norms, condition numbers
 * @module mini-linear-algebra-control
 * Knowledge coverage: L1 norms, L3 condition numbers
 */
#include "norms.h"
#include "decompositions.h"
#include "control_linalg.h"
#include "eigen.h"
#include <math.h>

/* Matrix norms */
double matrix_norm_frobenius(const Matrix *A) {
    if (!A) return 0.0;
    double sum = 0.0;
    int n = A->rows * A->cols;
    for (int k = 0; k < n; k++) {
        double v = A->data[k];
        sum += v * v;
    }
    return sqrt(sum);
}

double matrix_norm_induced_l1(const Matrix *A) {
    if (!A) return 0.0;
    double max_colsum = 0.0;
    for (int j = 0; j < A->cols; j++) {
        double colsum = 0.0;
        for (int i = 0; i < A->rows; i++)
            colsum += fabs(matrix_get(A, i, j));
        if (colsum > max_colsum) max_colsum = colsum;
    }
    return max_colsum;
}

double matrix_norm_induced_linf(const Matrix *A) {
    if (!A) return 0.0;
    double max_rowsum = 0.0;
    for (int i = 0; i < A->rows; i++) {
        double rowsum = 0.0;
        for (int j = 0; j < A->cols; j++)
            rowsum += fabs(matrix_get(A, i, j));
        if (rowsum > max_rowsum) max_rowsum = rowsum;
    }
    return max_rowsum;
}

double matrix_norm_spectral(const Matrix *A) {
    if (!A) return 0.0;
    /* = sigma_max(A), the largest singular value */
    int min_dim = (A->rows < A->cols) ? A->rows : A->cols;
    Matrix *Ac = matrix_create_copy(A);
    Matrix *U = matrix_create(A->rows, A->rows);
    Vector *S = vector_create(min_dim);
    Matrix *Vt = matrix_create(A->cols, A->cols);
    if (!Ac || !U || !S || !Vt) {
        matrix_free(Ac); matrix_free(U); vector_free(S); matrix_free(Vt); return 0.0;
    }
    double norm = 0.0;
    if (svd_jacobi(Ac, U, S, Vt, 1e-10, 100) >= 0 && S->size > 0)
        norm = fabs(vector_get(S, 0));
    matrix_free(Ac); matrix_free(U); vector_free(S); matrix_free(Vt);
    return norm;
}

double matrix_norm_nuclear(const Matrix *A) {
    if (!A) return 0.0;
    int min_dim = (A->rows < A->cols) ? A->rows : A->cols;
    Matrix *Ac = matrix_create_copy(A);
    Matrix *U = matrix_create(A->rows, A->rows);
    Vector *S = vector_create(min_dim);
    Matrix *Vt = matrix_create(A->cols, A->cols);
    if (!Ac || !U || !S || !Vt) {
        matrix_free(Ac); matrix_free(U); vector_free(S); matrix_free(Vt); return 0.0;
    }
    double trace_norm = 0.0;
    if (svd_jacobi(Ac, U, S, Vt, 1e-10, 100) >= 0)
        for (int i = 0; i < S->size; i++)
            trace_norm += fabs(vector_get(S, i));
    matrix_free(Ac); matrix_free(U); vector_free(S); matrix_free(Vt);
    return trace_norm;
}

double matrix_norm_l21(const Matrix *A) {
    if (!A) return 0.0;
    double sum = 0.0;
    for (int j = 0; j < A->cols; j++) {
        double col_norm = 0.0;
        for (int i = 0; i < A->rows; i++) {
            double v = matrix_get(A, i, j);
            col_norm += v * v;
        }
        sum += sqrt(col_norm);
    }
    return sum;
}

/* Condition numbers */
double matrix_condition_l2(const Matrix *A) {
    if (!A || !matrix_is_square(A)) return 0.0;
    int n = A->rows;
    Matrix *Ainv = matrix_create(n, n);
    if (!Ainv) return 0.0;
    if (lu_inverse(A, Ainv) != 0) { matrix_free(Ainv); return INFINITY; }
    double norm_A = matrix_norm_spectral(A);
    double norm_Ainv = matrix_norm_spectral(Ainv);
    matrix_free(Ainv);
    return (fabs(norm_Ainv) < 1e-300) ? INFINITY : norm_A * norm_Ainv;
}

double matrix_condition_l1(const Matrix *A) {
    if (!A || !matrix_is_square(A)) return 0.0;
    int n = A->rows;
    Matrix *Ainv = matrix_create(n, n);
    if (!Ainv) return 0.0;
    if (lu_inverse(A, Ainv) != 0) { matrix_free(Ainv); return INFINITY; }
    double c = matrix_norm_induced_l1(A) * matrix_norm_induced_l1(Ainv);
    matrix_free(Ainv);
    return c;
}

double matrix_condition_linf(const Matrix *A) {
    if (!A || !matrix_is_square(A)) return 0.0;
    int n = A->rows;
    Matrix *Ainv = matrix_create(n, n);
    if (!Ainv) return 0.0;
    if (lu_inverse(A, Ainv) != 0) { matrix_free(Ainv); return INFINITY; }
    double c = matrix_norm_induced_linf(A) * matrix_norm_induced_linf(Ainv);
    matrix_free(Ainv);
    return c;
}

double matrix_condition_frobenius(const Matrix *A) {
    if (!A || !matrix_is_square(A)) return 0.0;
    int n = A->rows;
    Matrix *Ainv = matrix_create(n, n);
    if (!Ainv) return 0.0;
    if (lu_inverse(A, Ainv) != 0) { matrix_free(Ainv); return INFINITY; }
    double c = matrix_norm_frobenius(A) * matrix_norm_frobenius(Ainv);
    matrix_free(Ainv);
    return c;
}

double matrix_condition_l2_estimate(const Matrix *A) {
    if (!A) return 0.0;
    /* Power iteration estimate using Hager-Higham estimator */
    /* Simplified: use SVD-based condition number */
    return matrix_condition_l2(A);
}

double controllability_condition_number(const Matrix *A, const Matrix *B) {
    if (!A || !B) return 0.0;
    int n = A->rows;
    Matrix *Cmat = matrix_create(n, n * B->cols);
    if (!Cmat) return 0.0;
    control_controllability_matrix(A, B, Cmat);
    double cond = matrix_condition_l2(Cmat);
    matrix_free(Cmat);
    return cond;
}

double lyapunov_condition_number(const Matrix *A, double tol) {
    if (!A) return 0.0;
    int n = A->rows;
    /* = 1 / sep(A, -A^T) where sep is the separation of spectra */
    /* Estimate via: kappa ~ ||I (x) A + A (x) I|| */
    /* Simplified: return 1/min |Re(eig(A))| */
    Vector *re = vector_create(n), *im = vector_create(n);
    if (!re || !im) { vector_free(re); vector_free(im); return INFINITY; }
    if (eigen_qr_general(A, re, im, tol, 100) < 0) { vector_free(re); vector_free(im); return INFINITY; }
    double min_re = INFINITY;
    for (int i = 0; i < n; i++) {
        double r = fabs(vector_get(re, i));
        if (r < min_re) min_re = r;
    }
    vector_free(re); vector_free(im);
    return (min_re < 1e-300) ? INFINITY : 1.0 / min_re;
}

/* Numerical stability measures */
double lu_growth_factor(const Matrix *A) {
    if (!A) return 0.0;
    Matrix *LU = matrix_create_copy(A);
    int n = A->rows;
    int *piv = (int*)malloc((size_t)n * sizeof(int));
    if (!LU || !piv) { matrix_free(LU); free(piv); return 0.0; }
    if (lu_decompose(LU, piv) != 0) { matrix_free(LU); free(piv); return INFINITY; }
    double max_U = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = i; j < n; j++) {
            double v = fabs(matrix_get(LU, i, j));
            if (v > max_U) max_U = v;
        }
    double max_A = 0.0;
    for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->cols; j++) {
            double v = fabs(matrix_get(A, i, j));
            if (v > max_A) max_A = v;
        }
    matrix_free(LU); free(piv);
    return (max_A < 1e-300) ? 0.0 : max_U / max_A;
}

double orthogonality_error(const Matrix *Q) {
    if (!Q) return INFINITY;
    int n = Q->rows;
    Matrix *Qt = matrix_create(n, n);
    Matrix *QtQ = matrix_create(n, n);
    if (!Qt || !QtQ) { matrix_free(Qt); matrix_free(QtQ); return INFINITY; }
    matrix_transpose(Q, Qt);
    matrix_multiply(Qt, Q, QtQ);
    double err = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double diff = matrix_get(QtQ, i, j) - ((i==j)?1.0:0.0);
            err += diff * diff;
        }
    matrix_free(Qt); matrix_free(QtQ);
    return sqrt(err);
}

double linear_system_residual(const Matrix *A, const Vector *x, const Vector *b) {
    if (!A || !x || !b) return INFINITY;
    int m = A->rows;
    Vector *r = vector_create(m);
    if (!r) return INFINITY;
    for (int i = 0; i < m; i++) {
        double Ax_i = 0.0;
        for (int j = 0; j < A->cols; j++)
            Ax_i += matrix_get(A, i, j) * vector_get(x, j);
        vector_set(r, i, vector_get(b, i) - Ax_i);
    }
    double resid = vector_norm_l2(r);
    vector_free(r);
    return resid;
}

double lyapunov_residual(const Matrix *A, const Matrix *X, const Matrix *Q) {
    if (!A || !X || !Q) return INFINITY;
    int n = A->rows;
    Matrix *AX = matrix_create(n, n);
    Matrix *XAt = matrix_create(n, n);
    Matrix *At = matrix_create(n, n);
    if (!AX || !XAt || !At) {
        matrix_free(AX); matrix_free(XAt); matrix_free(At); return INFINITY;
    }
    matrix_transpose(A, At);
    matrix_multiply(A, X, AX);
    matrix_multiply(X, At, XAt);
    double resid = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double r = matrix_get(AX, i, j) + matrix_get(XAt, i, j) + matrix_get(Q, i, j);
            resid += r * r;
        }
    matrix_free(AX); matrix_free(XAt); matrix_free(At);
    return sqrt(resid);
}