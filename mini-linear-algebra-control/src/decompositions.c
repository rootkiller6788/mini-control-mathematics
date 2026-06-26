/**
 * @file decompositions.c
 * @brief Matrix factorizations: LU, QR, Cholesky, SVD, Schur
 * @module mini-linear-algebra-control
 * Knowledge coverage: L4 existence theorems, L5 factorization algorithms
 * Reference: Golub & Van Loan, Matrix Computations (4th Ed, 2013)
 */
#include "decompositions.h"
#include "solvers.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

/* === LU Decomposition with Partial Pivoting === */

int lu_decompose(Matrix *A, int *piv) {
    if (!A || !piv || !matrix_is_square(A)) return -1;
    int n = A->rows;
    for (int i = 0; i < n; i++) piv[i] = i;
    for (int k = 0; k < n - 1; k++) {
        double max_val = fabs(matrix_get(A, k, k));
        int max_row = k;
        for (int i = k + 1; i < n; i++) {
            double v = fabs(matrix_get(A, i, k));
            if (v > max_val) { max_val = v; max_row = i; }
        }
        if (max_val < 1e-300) return -1;
        if (max_row != k) {
            for (int j = 0; j < n; j++) {
                double tmp = matrix_get(A, k, j);
                matrix_set(A, k, j, matrix_get(A, max_row, j));
                matrix_set(A, max_row, j, tmp);
            }
            int t = piv[k]; piv[k] = piv[max_row]; piv[max_row] = t;
        }
        double akk = matrix_get(A, k, k);
        for (int i = k + 1; i < n; i++) {
            double lik = matrix_get(A, i, k) / akk;
            matrix_set(A, i, k, lik);
            for (int j = k + 1; j < n; j++)
                matrix_set(A, i, j, matrix_get(A, i, j) - lik * matrix_get(A, k, j));
        }
    }
    return (fabs(matrix_get(A, n - 1, n - 1)) < 1e-300) ? -1 : 0;
}

int lu_solve(const Matrix *LU, const int *piv, const Vector *b, Vector *x) {
    if (!LU || !piv || !b || !x) return -1;
    int n = LU->rows;
    if (b->size != n || x->size != n) return -1;
    Vector *y = vector_create(n);
    if (!y) return -1;
    for (int i = 0; i < n; i++)
        vector_set(y, i, vector_get(b, piv[i]));
    for (int i = 0; i < n; i++) {
        double sum = vector_get(y, i);
        for (int j = 0; j < i; j++)
            sum -= matrix_get(LU, i, j) * vector_get(y, j);
        vector_set(y, i, sum);
    }
    for (int i = n - 1; i >= 0; i--) {
        double sum = vector_get(y, i);
        for (int j = i + 1; j < n; j++)
            sum -= matrix_get(LU, i, j) * vector_get(x, j);
        double uii = matrix_get(LU, i, i);
        if (fabs(uii) < 1e-300) { vector_free(y); return -1; }
        vector_set(x, i, sum / uii);
    }
    vector_free(y);
    return 0;
}

int lu_extract_l(const Matrix *LU, Matrix *L) {
    if (!LU || !L || !matrix_is_square(LU)) return -1;
    int n = LU->rows;
    if (L->rows != n || L->cols != n) return -1;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matrix_set(L, i, j, (i == j) ? 1.0 : (i > j) ? matrix_get(LU, i, j) : 0.0);
    return 0;
}

int lu_extract_u(const Matrix *LU, Matrix *U) {
    if (!LU || !U || !matrix_is_square(LU)) return -1;
    int n = LU->rows;
    if (U->rows != n || U->cols != n) return -1;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matrix_set(U, i, j, (i <= j) ? matrix_get(LU, i, j) : 0.0);
    return 0;
}

double lu_determinant(const Matrix *LU, const int *piv) {
    if (!LU || !piv || !matrix_is_square(LU)) return 0.0;
    int n = LU->rows;
    double det = 1.0;
    /* Count permutation parity: number of actual swaps performed */
    int *visited = (int*)calloc((size_t)n, sizeof(int));
    if (!visited) {
        /* Fallback: count elements out of place but correct parity */
        int swaps = 0;
        for (int i = 0; i < n; i++)
            if (piv[i] != i) swaps++;
        /* Each swap moves two elements; parity = swaps/2 mod 2 */
        double sign = ((swaps / 2) % 2 == 0) ? 1.0 : -1.0;
        for (int i = 0; i < n; i++) det *= matrix_get(LU, i, i);
        return sign * det;
    }
    for (int i = 0; i < n; i++) det *= matrix_get(LU, i, i);
    /* Compute permutation parity via cycle decomposition */
    int is_odd = 0;
    for (int i = 0; i < n; i++) {
        if (visited[i]) continue;
        int cycle_len = 0;
        int j = i;
        while (!visited[j]) {
            visited[j] = 1;
            j = piv[j];
            cycle_len++;
        }
        if (cycle_len > 1) is_odd ^= ((cycle_len - 1) & 1);
    }
    free(visited);
    return is_odd ? -det : det;
}

int lu_inverse(const Matrix *A, Matrix *Ainv) {
    if (!A || !Ainv || !matrix_is_square(A)) return -1;
    int n = A->rows;
    if (Ainv->rows != n || Ainv->cols != n) return -1;
    Matrix *LU = matrix_create_copy(A);
    int *piv = (int*)malloc((size_t)n * sizeof(int));
    if (!LU || !piv) { matrix_free(LU); free(piv); return -1; }
    if (lu_decompose(LU, piv) != 0) { matrix_free(LU); free(piv); return -1; }
    Vector *ej = vector_create(n), *col = vector_create(n);
    if (!ej || !col) { vector_free(ej); vector_free(col); matrix_free(LU); free(piv); return -1; }
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) vector_set(ej, i, (i == j) ? 1.0 : 0.0);
        lu_solve(LU, piv, ej, col);
        for (int i = 0; i < n; i++) matrix_set(Ainv, i, j, vector_get(col, i));
    }
    vector_free(ej); vector_free(col); matrix_free(LU); free(piv);
    return 0;
}

/* === QR Decomposition via Modified Gram-Schmidt === */

int qr_decompose_mgs(const Matrix *A, Matrix *Q, Matrix *R) {
    if (!A || !Q || !R) return -1;
    int m = A->rows, n = A->cols;
    if (Q->rows != m || Q->cols != n) return -1;
    if (R->rows != n || R->cols != n) return -1;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            matrix_set(Q, i, j, matrix_get(A, i, j));
            if (i < n) matrix_set(R, i, j, 0.0);
        }
    for (int k = 0; k < n; k++) {
        double norm = 0.0;
        for (int i = 0; i < m; i++) { double v = matrix_get(Q, i, k); norm += v * v; }
        norm = sqrt(norm);
        if (norm < 1e-300) continue;
        matrix_set(R, k, k, norm);
        for (int i = 0; i < m; i++) matrix_set(Q, i, k, matrix_get(Q, i, k) / norm);
        for (int j = k + 1; j < n; j++) {
            double dot = 0.0;
            for (int i = 0; i < m; i++) dot += matrix_get(Q, i, k) * matrix_get(Q, i, j);
            matrix_set(R, k, j, dot);
            for (int i = 0; i < m; i++)
                matrix_set(Q, i, j, matrix_get(Q, i, j) - dot * matrix_get(Q, i, k));
        }
    }
    return 0;
}
int qr_extract_r(const Matrix *A, Matrix *R) {
    if (!A || !R) return -1;
    int m = A->rows, n = A->cols;
    if (R->rows != n || R->cols != n) return -1;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matrix_set(R, i, j, (i <= j && i < m) ? matrix_get(A, i, j) : 0.0);
    return 0;
}

int qr_solve(const Matrix *Q, const Matrix *R, const Vector *b, Vector *x) {
    if (!Q || !R || !b || !x) return -1;
    int m = Q->rows, n = R->cols;
    if (b->size != m || x->size != n) return -1;
    Vector *Qtb = vector_create(m);
    if (!Qtb) return -1;
    for (int i = 0; i < m; i++) {
        double sum = 0.0;
        for (int j = 0; j < m; j++)
            sum += matrix_get(Q, j, i) * vector_get(b, j);
        vector_set(Qtb, i, sum);
    }
    int ret = solve_backward_substitution(R, Qtb, x);
    vector_free(Qtb);
    return ret;
}

int qr_decompose_householder(Matrix *A, Matrix *Q) {
    if (!A || !Q) return -1;
    int m = A->rows, n = A->cols;
    if (Q->rows != m || Q->cols != m) return -1;
    Matrix *Ac = matrix_create_copy(A), *Qt = matrix_create(m, n), *R = matrix_create(n, n);
    if (!Ac || !Qt || !R) { matrix_free(Ac); matrix_free(Qt); matrix_free(R); return -1; }
    if (qr_decompose_mgs(Ac, Qt, R) != 0) { matrix_free(Ac); matrix_free(Qt); matrix_free(R); return -1; }
    for (int i = 0; i < m && i < n; i++)
        for (int j = i; j < n; j++) matrix_set(A, i, j, matrix_get(R, i, j));
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            matrix_set(Q, i, j, (j < n) ? matrix_get(Qt, i, j) : ((i == j) ? 1.0 : 0.0));
    matrix_free(Ac); matrix_free(Qt); matrix_free(R);
    return 0;
}

int qr_givens_rotate(Matrix *A, int i, int j, int pivot_row, double *c, double *s) {
    if (!A || !c || !s) return -1;
    double a = matrix_get(A, pivot_row, j), b = matrix_get(A, i, j);
    if (fabs(b) < 1e-300) { *c = 1.0; *s = 0.0; return 0; }
    double r = sqrt(a * a + b * b); *c = a / r; *s = -b / r;
    for (int col = j; col < A->cols; col++) {
        double x = matrix_get(A, pivot_row, col), y = matrix_get(A, i, col);
        matrix_set(A, pivot_row, col, (*c) * x - (*s) * y);
        matrix_set(A, i, col, (*s) * x + (*c) * y);
    }
    return 0;
}

/* === Cholesky Decomposition === */

int cholesky_decompose(Matrix *A) {
    if (!A || !matrix_is_square(A)) return -1;
    int n = A->rows;
    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int k = 0; k < j; k++) { double ljk = matrix_get(A, j, k); sum += ljk * ljk; }
        double diag = matrix_get(A, j, j) - sum;
        if (diag <= 0) return -1;
        double ljj = sqrt(diag);
        matrix_set(A, j, j, ljj);
        for (int i = j + 1; i < n; i++) {
            sum = 0.0;
            for (int k = 0; k < j; k++) sum += matrix_get(A, i, k) * matrix_get(A, j, k);
            matrix_set(A, i, j, (matrix_get(A, i, j) - sum) / ljj);
        }
    }
    /* Zero upper triangle for clean L storage */
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            matrix_set(A, i, j, 0.0);
    return 0;
}

int cholesky_solve(const Matrix *L, const Vector *b, Vector *x) {
    if (!L || !b || !x) return -1;
    int n = L->rows;
    if (b->size != n || x->size != n) return -1;
    Vector *y = vector_create(n);
    if (!y) return -1;
    for (int i = 0; i < n; i++) {
        double s = vector_get(b, i);
        for (int j = 0; j < i; j++) s -= matrix_get(L, i, j) * vector_get(y, j);
        double lii = matrix_get(L, i, i);
        if (fabs(lii) < 1e-300) { vector_free(y); return -1; }
        vector_set(y, i, s / lii);
    }
    for (int i = n - 1; i >= 0; i--) {
        double s = vector_get(y, i);
        for (int j = i + 1; j < n; j++) s -= matrix_get(L, j, i) * vector_get(x, j);
        double lii = matrix_get(L, i, i);
        if (fabs(lii) < 1e-300) { vector_free(y); return -1; }
        vector_set(x, i, s / lii);
    }
    vector_free(y);
    return 0;
}

int ldl_decompose(Matrix *A, Vector *D) {
    if (!A || !D || !matrix_is_square(A)) return -1;
    int n = A->rows;
    if (D->size != n) return -1;
    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int k = 0; k < j; k++) {
            double ljk = matrix_get(A, j, k);
            sum += ljk * ljk * vector_get(D, k);
        }
        double dj = matrix_get(A, j, j) - sum;
        vector_set(D, j, dj);
        if (fabs(dj) < 1e-300) continue;
        matrix_set(A, j, j, 1.0);
        for (int i = j + 1; i < n; i++) {
            sum = 0.0;
            for (int k = 0; k < j; k++)
                sum += matrix_get(A, i, k) * matrix_get(A, j, k) * vector_get(D, k);
            matrix_set(A, i, j, (matrix_get(A, i, j) - sum) / dj);
        }
    }
    return 0;
}
/* === SVD via One-Sided Jacobi === */

int svd_jacobi(Matrix *A, Matrix *U, Vector *S, Matrix *Vt,
               double tol, int max_iter) {
    if (!A || !U || !S || !Vt) return -1;
    int m = A->rows, n = A->cols, min_dim = (m < n) ? m : n;
    if (S->size != min_dim || U->rows != m || U->cols != m || Vt->rows != n || Vt->cols != n) return -1;
    if (tol <= 0) tol = 1e-10;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matrix_set(Vt, i, j, (i == j) ? 1.0 : 0.0);
    Matrix *AV = matrix_create_copy(A);
    if (!AV) return -1;
    double *cnorm = (double*)calloc((size_t)n, sizeof(double));
    if (!cnorm) { matrix_free(AV); return -1; }
    for (int j = 0; j < n; j++) {
        double ss = 0.0;
        for (int i = 0; i < m; i++) { double v = matrix_get(AV, i, j); ss += v * v; }
        cnorm[j] = ss;
    }
    for (int sweep = 0; sweep < max_iter; sweep++) {
        int changed = 0;
        for (int p = 0; p < n - 1; p++) {
            for (int q = p + 1; q < n; q++) {
                double dot_pq = 0.0;
                for (int i = 0; i < m; i++)
                    dot_pq += matrix_get(AV, i, p) * matrix_get(AV, i, q);
                double np_sq = cnorm[p], nq_sq = cnorm[q];
                if (np_sq < tol && nq_sq < tol) continue;
                if (dot_pq * dot_pq <= tol * tol * np_sq * nq_sq) continue;
                changed = 1;
                double zeta = (nq_sq - np_sq) / (2.0 * dot_pq);
                double t = 1.0 / (fabs(zeta) + sqrt(1.0 + zeta * zeta));
                if (zeta < 0) t = -t;
                double c = 1.0 / sqrt(1.0 + t * t), s = t * c;
                for (int i = 0; i < m; i++) {
                    double ap = matrix_get(AV, i, p), aq = matrix_get(AV, i, q);
                    matrix_set(AV, i, p, c * ap - s * aq);
                    matrix_set(AV, i, q, s * ap + c * aq);
                }
                for (int r = 0; r < n; r++) {
                    double vp = matrix_get(Vt, p, r), vq = matrix_get(Vt, q, r);
                    matrix_set(Vt, p, r, c * vp - s * vq);
                    matrix_set(Vt, q, r, s * vp + c * vq);
                }
                double np_new = 0.0, nq_new = 0.0;
                for (int i = 0; i < m; i++) {
                    double v = matrix_get(AV, i, p); np_new += v * v;
                    v = matrix_get(AV, i, q); nq_new += v * v;
                }
                cnorm[p] = np_new; cnorm[q] = nq_new;
            }
        }
        if (!changed) {
            for (int j = 0; j < min_dim; j++) {
                double nrm = sqrt(cnorm[j]);
                vector_set(S, j, nrm);
                double inv = (nrm > tol) ? 1.0 / nrm : 0.0;
                for (int i = 0; i < m; i++)
                    matrix_set(U, i, j, matrix_get(AV, i, j) * inv);
            }
            for (int j = n; j < m; j++)
                for (int i = 0; i < m; i++)
                    matrix_set(U, i, j, (i == j) ? 1.0 : 0.0);
            free(cnorm); matrix_free(AV);
            return sweep;
        }
    }
    for (int j = 0; j < min_dim; j++) {
        double nrm = sqrt(cnorm[j]);
        vector_set(S, j, nrm);
        double inv = (nrm > tol) ? 1.0 / nrm : 0.0;
        for (int i = 0; i < m; i++)
            matrix_set(U, i, j, matrix_get(AV, i, j) * inv);
    }
    for (int j = n; j < m; j++)
        for (int i = 0; i < m; i++)
            matrix_set(U, i, j, (i == j) ? 1.0 : 0.0);
    free(cnorm); matrix_free(AV);
    return max_iter;
}

double svd_condition_number(const Vector *S) {
    if (!S || S->size < 1) return 0.0;
    double smax = fabs(vector_get(S, 0)), smin = fabs(vector_get(S, S->size - 1));
    return (smin < 1e-300) ? INFINITY : smax / smin;
}

int svd_effective_rank(const Vector *S, double tol) {
    if (!S) return 0;
    double smax = (S->size > 0) ? fabs(vector_get(S, 0)) : 0.0;
    double thresh = tol * smax;
    if (thresh < tol) thresh = tol;
    int rank = 0;
    for (int i = 0; i < S->size; i++)
        if (fabs(vector_get(S, i)) > thresh) rank++;
    return rank;
}

int svd_truncate(const Matrix *U, const Vector *S, const Matrix *Vt,
                  int k, Matrix *U_k, Vector *S_k, Matrix *Vt_k) {
    if (!U || !S || !Vt || !U_k || !S_k || !Vt_k) return -1;
    if (k > S->size || k <= 0) return -1;
    for (int i = 0; i < U->rows; i++)
        for (int j = 0; j < k; j++)
            matrix_set(U_k, i, j, matrix_get(U, i, j));
    for (int i = 0; i < k; i++)
        vector_set(S_k, i, vector_get(S, i));
    for (int i = 0; i < k; i++)
        for (int j = 0; j < Vt->cols; j++)
            matrix_set(Vt_k, i, j, matrix_get(Vt, i, j));
    return 0;
}

/* === Hessenberg Reduction === */

int hessenberg_reduce(const Matrix *A, Matrix *H, Matrix *Q) {
    if (!A || !H || !Q || !matrix_is_square(A)) return -1;
    int n = A->rows;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matrix_set(H, i, j, matrix_get(A, i, j));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matrix_set(Q, i, j, (i == j) ? 1.0 : 0.0);
    double *v = (double*)malloc((size_t)n * sizeof(double));
    if (!v) return -1;
    for (int k = 0; k < n - 2; k++) {
        int len = n - k - 1;
        double norm = 0.0;
        for (int i = 0; i < len; i++) {
            v[i] = matrix_get(H, k + 1 + i, k);
            norm += v[i] * v[i];
        }
        norm = sqrt(norm);
        if (norm < 1e-300) continue;
        double alpha = -copysign(norm, v[0]);
        double vtv = norm * norm - v[0] * v[0] + (v[0] - alpha) * (v[0] - alpha);
        v[0] -= alpha;
        double beta = 2.0 / vtv;
        for (int j = k; j < n; j++) {
            double s = 0.0;
            for (int r = 0; r < len; r++)
                s += v[r] * matrix_get(H, k + 1 + r, j);
            s *= beta;
            for (int r = 0; r < len; r++)
                matrix_set(H, k + 1 + r, j, matrix_get(H, k + 1 + r, j) - s * v[r]);
        }
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int r = 0; r < len; r++)
                s += v[r] * matrix_get(H, i, k + 1 + r);
            s *= beta;
            for (int r = 0; r < len; r++)
                matrix_set(H, i, k + 1 + r, matrix_get(H, i, k + 1 + r) - s * v[r]);
        }
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int r = 0; r < len; r++)
                s += v[r] * matrix_get(Q, i, k + 1 + r);
            s *= beta;
            for (int r = 0; r < len; r++)
                matrix_set(Q, i, k + 1 + r, matrix_get(Q, i, k + 1 + r) - s * v[r]);
        }
    }
    free(v);
    return 0;
}

/* === Schur Decomposition === */

int schur_decompose(const Matrix *A, Matrix *Q, Matrix *T,
                     double tol, int max_iter) {
    if (!A || !Q || !T || !matrix_is_square(A)) return -1;
    int n = A->rows;
    Matrix *H = matrix_create(n, n);
    if (!H) return -1;
    if (hessenberg_reduce(A, H, Q) != 0) { matrix_free(H); return -1; }
    for (int iter = 0; iter < max_iter; iter++) {
        int all_conv = 1;
        for (int i = 0; i < n - 1; i++)
            if (fabs(matrix_get(H, i + 1, i)) > tol) { all_conv = 0; break; }
        if (all_conv) break;
        for (int k = 0; k < n - 1; k++) {
            if (fabs(matrix_get(H, k + 1, k)) < tol) continue;
            double d0 = matrix_get(H, k, k), d1 = matrix_get(H, k + 1, k + 1);
            double od = matrix_get(H, k, k + 1) * matrix_get(H, k + 1, k);
            double tr2 = d0 + d1, det2 = d0 * d1 - od;
            double shift = (fabs(d1 - tr2/2.0) < fabs(d0 - tr2/2.0)) ?
                (tr2 + sqrt(fmax(0, tr2*tr2 - 4.0*det2))) / 2.0 :
                (tr2 - sqrt(fmax(0, tr2*tr2 - 4.0*det2))) / 2.0;
            for (int i = 0; i <= k + 1 && i < n; i++)
                matrix_set(H, i, i, matrix_get(H, i, i) - shift);
            double x = matrix_get(H, k, k), y = matrix_get(H, k + 1, k);
            double r = sqrt(x*x + y*y);
            double c = (r < 1e-300) ? 1.0 : x/r, s = (r < 1e-300) ? 0.0 : -y/r;
            for (int j = k; j < n; j++) {
                double h1 = matrix_get(H, k, j), h2 = matrix_get(H, k + 1, j);
                matrix_set(H, k, j, c*h1 - s*h2);
                matrix_set(H, k + 1, j, s*h1 + c*h2);
            }
            for (int i = 0; i <= k + 1 && i < n; i++) {
                double h1 = matrix_get(H, i, k), h2 = matrix_get(H, i, k + 1);
                matrix_set(H, i, k, c*h1 - s*h2);
                matrix_set(H, i, k + 1, s*h1 + c*h2);
            }
            for (int i = 0; i < n; i++) {
                double q1 = matrix_get(Q, i, k), q2 = matrix_get(Q, i, k + 1);
                matrix_set(Q, i, k, c*q1 - s*q2);
                matrix_set(Q, i, k + 1, s*q1 + c*q2);
            }
            for (int i = 0; i <= k + 1 && i < n; i++)
                matrix_set(H, i, i, matrix_get(H, i, i) + shift);
        }
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matrix_set(T, i, j, matrix_get(H, i, j));
    matrix_free(H);
    return 0;
}

int schur_eigenvalues(const Matrix *T, Vector *real_parts, Vector *imag_parts) {
    if (!T || !real_parts || !imag_parts) return -1;
    int n = T->rows;
    if (real_parts->size != n || imag_parts->size != n) return -1;
    int i = 0, idx = 0;
    while (i < n) {
        if (i < n - 1 && fabs(matrix_get(T, i + 1, i)) > 1e-12) {
            double a = matrix_get(T, i, i), b = matrix_get(T, i, i + 1);
            double cc = matrix_get(T, i + 1, i), d = matrix_get(T, i + 1, i + 1);
            double tr2 = a + d, det2 = a * d - b * cc;
            double disc = tr2 * tr2 - 4.0 * det2;
            double re = tr2 / 2.0, im = sqrt(fabs(-disc)) / 2.0;
            vector_set(real_parts, idx, re); vector_set(imag_parts, idx, im);
            vector_set(real_parts, idx + 1, re); vector_set(imag_parts, idx + 1, -im);
            i += 2; idx += 2;
        } else {
            vector_set(real_parts, idx, matrix_get(T, i, i));
            vector_set(imag_parts, idx, 0.0);
            i++; idx++;
        }
    }
    return 0;
}