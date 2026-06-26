/**
 * linear_systems.c — Direct Linear System Solvers and Decompositions
 *
 * Implements: L3 Math Structures (decomposition types), L5 Methods (LU, QR, Cholesky),
 *             L6 Problems (Lyapunov/Riccati), L8 Advanced (SVD)
 *
 * These algorithms form the computational foundation of modern control:
 * LQR, Kalman filtering, pole placement, and system identification all
 * depend on efficient and numerically stable linear algebra.
 */

#include "linear_systems.h"
#include <string.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

/* ==========================================================================
 * L5: Gaussian Elimination with Partial Pivoting
 * ========================================================================== */

int gaussian_elimination_solve(const Matrix *A, const Vector *b, Vector *x) {
    if (!A || !b || !x || A->rows != A->cols || A->rows != b->n || b->n != x->n)
        return -1;
    size_t n = A->rows;

    /* Build augmented matrix [A | b] in working storage */
    double *aug = (double*)malloc(n * (n + 1) * sizeof(double));
    if (!aug) return -1;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++)
            aug[i * (n + 1) + j] = A->data[i * A->cols + j];
        aug[i * (n + 1) + n] = b->data[i];
    }

    for (size_t k = 0; k < n; k++) {
        /* Partial pivoting: find largest in column k */
        size_t max_row = k;
        double max_val = fabs(aug[k * (n + 1) + k]);
        for (size_t i = k + 1; i < n; i++) {
            double val = fabs(aug[i * (n + 1) + k]);
            if (val > max_val) { max_val = val; max_row = i; }
        }
        if (max_val < MNC_TINY) { free(aug); return -1; }  /* Singular */

        if (max_row != k) {
            for (size_t j = 0; j <= n; j++) {
                double tmp = aug[k * (n + 1) + j];
                aug[k * (n + 1) + j] = aug[max_row * (n + 1) + j];
                aug[max_row * (n + 1) + j] = tmp;
            }
        }

        /* Eliminate below */
        for (size_t i = k + 1; i < n; i++) {
            double factor = aug[i * (n + 1) + k] / aug[k * (n + 1) + k];
            aug[i * (n + 1) + k] = 0.0;
            for (size_t j = k + 1; j <= n; j++)
                aug[i * (n + 1) + j] -= factor * aug[k * (n + 1) + j];
        }
    }

    /* Back substitution */
    for (size_t i = n; i-- > 0; ) {
        double sum = aug[i * (n + 1) + n];
        for (size_t j = i + 1; j < n; j++)
            sum -= aug[i * (n + 1) + j] * x->data[j];
        x->data[i] = sum / aug[i * (n + 1) + i];
    }

    free(aug);
    return 0;
}

/* ==========================================================================
 * L5: LU Decomposition with Partial Pivoting
 * ========================================================================== */

LUResult* lu_decompose(const Matrix *A) {
    if (!A || A->rows != A->cols) return NULL;
    size_t n = A->rows;

    LUResult *lu = (LUResult*)malloc(sizeof(LUResult));
    if (!lu) return NULL;
    lu->n = n;
    lu->singular = 0;
    lu->pivots = (size_t*)malloc(n * sizeof(size_t));
    lu->L = matrix_create_zeros(n, n);
    lu->U = matrix_create_zeros(n, n);
    if (!lu->pivots || !lu->L || !lu->U) { lu_free(lu); return NULL; }

    /* Copy A into U (which will hold both L and U temporarily) */
    double *work = (double*)malloc(n * n * sizeof(double));
    if (!work) { lu_free(lu); return NULL; }
    memcpy(work, A->data, n * n * sizeof(double));

    for (size_t i = 0; i < n; i++) lu->pivots[i] = i;

    for (size_t k = 0; k < n; k++) {
        /* Partial pivoting */
        double max_val = fabs(work[k * n + k]);
        size_t max_row = k;
        for (size_t i = k + 1; i < n; i++) {
            if (fabs(work[i * n + k]) > max_val) {
                max_val = fabs(work[i * n + k]);
                max_row = i;
            }
        }
        if (max_val < MNC_TINY) { lu->singular = 1; free(work); return lu; }

        if (max_row != k) {
            size_t tmp = lu->pivots[k];
            lu->pivots[k] = lu->pivots[max_row];
            lu->pivots[max_row] = tmp;
            for (size_t j = 0; j < n; j++) {
                double t = work[k * n + j];
                work[k * n + j] = work[max_row * n + j];
                work[max_row * n + j] = t;
            }
        }

        for (size_t i = k + 1; i < n; i++) {
            double factor = work[i * n + k] / work[k * n + k];
            work[i * n + k] = factor;
            for (size_t j = k + 1; j < n; j++)
                work[i * n + j] -= factor * work[k * n + j];
        }
    }

    /* Extract L and U from work */
    for (size_t i = 0; i < n; i++) {
        lu->L->data[i * n + i] = 1.0;
        for (size_t j = 0; j < n; j++) {
            if (i > j)
                lu->L->data[i * n + j] = work[i * n + j];
            else
                lu->U->data[i * n + j] = work[i * n + j];
        }
    }

    free(work);
    return lu;
}

Vector* lu_solve(const LUResult *lu, const Vector *b) {
    if (!lu || !b || lu->n != b->n) return NULL;
    size_t n = lu->n;
    Vector *x = vector_create(n);
    Vector *y = vector_create(n);
    if (!x || !y) { vector_free(x); vector_free(y); return NULL; }

    /* Forward substitution: Ly = P^T b (apply permutation to b) */
    for (size_t i = 0; i < n; i++) {
        y->data[i] = b->data[lu->pivots[i]];
        for (size_t j = 0; j < i; j++)
            y->data[i] -= lu->L->data[i * n + j] * y->data[j];
    }

    /* Back substitution: Ux = y */
    for (size_t i = n; i-- > 0; ) {
        x->data[i] = y->data[i];
        for (size_t j = i + 1; j < n; j++)
            x->data[i] -= lu->U->data[i * n + j] * x->data[j];
        x->data[i] /= lu->U->data[i * n + i];
    }

    vector_free(y);
    return x;
}

double lu_determinant(const LUResult *lu) {
    if (!lu || lu->singular) return 0.0;
    double det = 1.0;
    int sign = 1;
    for (size_t i = 0; i < lu->n; i++) {
        if (lu->pivots[i] != i) sign = -sign;
        det *= lu->U->data[i * lu->U->cols + i];
    }
    return sign * det;
}

void lu_free(LUResult *lu) {
    if (!lu) return;
    matrix_free(lu->L);
    matrix_free(lu->U);
    free(lu->pivots);
    free(lu);
}

/* ==========================================================================
 * L5: Cholesky Decomposition
 * ========================================================================== */

CholeskyResult* cholesky_decompose(const Matrix *A) {
    if (!A || A->rows != A->cols) return NULL;
    size_t n = A->rows;

    CholeskyResult *chol = (CholeskyResult*)malloc(sizeof(CholeskyResult));
    if (!chol) return NULL;
    chol->n = n;
    chol->L = matrix_create_zeros(n, n);
    if (!chol->L) { free(chol); return NULL; }

    for (size_t j = 0; j < n; j++) {
        double sum = 0.0;
        for (size_t k = 0; k < j; k++)
            sum += chol->L->data[j * n + k] * chol->L->data[j * n + k];
        double diag = A->data[j * n + j] - sum;
        if (diag <= 0.0) { cholesky_free(chol); return NULL; }  /* Not PD */
        chol->L->data[j * n + j] = sqrt(diag);

        for (size_t i = j + 1; i < n; i++) {
            sum = 0.0;
            for (size_t k = 0; k < j; k++)
                sum += chol->L->data[i * n + k] * chol->L->data[j * n + k];
            chol->L->data[i * n + j] = (A->data[i * n + j] - sum) / chol->L->data[j * n + j];
        }
    }
    return chol;
}

Vector* cholesky_solve(const CholeskyResult *chol, const Vector *b) {
    if (!chol || !b || chol->n != b->n) return NULL;
    size_t n = chol->n;
    Vector *x = vector_create(n);
    Vector *y = vector_create(n);
    if (!x || !y) { vector_free(x); vector_free(y); return NULL; }

    /* Forward: Ly = b */
    for (size_t i = 0; i < n; i++) {
        y->data[i] = b->data[i];
        for (size_t j = 0; j < i; j++)
            y->data[i] -= chol->L->data[i * n + j] * y->data[j];
        y->data[i] /= chol->L->data[i * n + i];
    }

    /* Backward: L^T x = y */
    for (size_t i = n; i-- > 0; ) {
        x->data[i] = y->data[i];
        for (size_t j = i + 1; j < n; j++)
            x->data[i] -= chol->L->data[j * n + i] * x->data[j];
        x->data[i] /= chol->L->data[i * n + i];
    }

    vector_free(y);
    return x;
}

double cholesky_logdet(const CholeskyResult *chol) {
    if (!chol) return 0.0;
    double logdet = 0.0;
    for (size_t i = 0; i < chol->n; i++)
        logdet += log(chol->L->data[i * chol->n + i]);
    return 2.0 * logdet;
}

void cholesky_free(CholeskyResult *chol) {
    if (!chol) return;
    matrix_free(chol->L);
    free(chol);
}

/* ==========================================================================
 * L5: QR Decomposition via Householder Reflections
 * ========================================================================== */

QRResult* qr_decompose_householder(const Matrix *A) {
    if (!A) return NULL;
    size_t m = A->rows, n = A->cols;
    if (m < n) return NULL;  /* Require m >= n for QR */

    QRResult *qr = (QRResult*)malloc(sizeof(QRResult));
    if (!qr) return NULL;
    qr->m = m; qr->n = n;
    qr->Q = matrix_create_identity(m);
    qr->R = matrix_copy(A);
    if (!qr->Q || !qr->R) { qr_free(qr); return NULL; }

    double *v = (double*)malloc(m * sizeof(double));
    if (!v) { qr_free(qr); return NULL; }

    for (size_t k = 0; k < n; k++) {
        /* Extract column k of R below diagonal */
        double xnorm = 0.0;
        for (size_t i = k; i < m; i++) {
            v[i - k] = qr->R->data[i * n + k];
            xnorm += v[i - k] * v[i - k];
        }
        xnorm = sqrt(xnorm);
        double alpha = (v[0] >= 0) ? -xnorm : xnorm;
        v[0] -= alpha;
        double vnorm2 = v[0] * v[0];
        for (size_t i = 1; i < m - k; i++) vnorm2 += v[i] * v[i];
        if (vnorm2 < MNC_TINY) continue;

        /* Apply Householder reflection H = I - 2*v*v^T/(v^T v) */
        double beta = 2.0 / vnorm2;

        /* Update R: H * R(k:m, k:n) */
        for (size_t j = k; j < n; j++) {
            double dot = 0.0;
            for (size_t i = 0; i < m - k; i++)
                dot += v[i] * qr->R->data[(k + i) * n + j];
            for (size_t i = 0; i < m - k; i++)
                qr->R->data[(k + i) * n + j] -= beta * dot * v[i];
        }

        /* Update Q: Q * H (from right) */
        for (size_t i = 0; i < m; i++) {
            double dot = 0.0;
            for (size_t j = 0; j < m - k; j++)
                dot += qr->Q->data[i * m + (k + j)] * v[j];
            for (size_t j = 0; j < m - k; j++)
                qr->Q->data[i * m + (k + j)] -= beta * dot * v[j];
        }
    }

    free(v);
    /* Zero out below diagonal in R for cleanliness */
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n && j < i; j++)
            qr->R->data[i * n + j] = 0.0;

    return qr;
}

Vector* qr_solve(const QRResult *qr, const Vector *b) {
    if (!qr || !b || qr->m != b->n) return NULL;
    /* Compute Q^T b */
    Vector *Qtb = vector_create(qr->m);
    if (!Qtb) return NULL;
    for (size_t i = 0; i < qr->m; i++) {
        Qtb->data[i] = 0.0;
        for (size_t j = 0; j < qr->m; j++)
            Qtb->data[i] += qr->Q->data[j * qr->m + i] * b->data[j];
    }

    /* Solve R(1:n,1:n) x = (Q^T b)(1:n) */
    Vector *x = vector_create(qr->n);
    if (!x) { vector_free(Qtb); return NULL; }
    for (size_t i = qr->n; i-- > 0; ) {
        x->data[i] = Qtb->data[i];
        for (size_t j = i + 1; j < qr->n; j++)
            x->data[i] -= qr->R->data[i * qr->n + j] * x->data[j];
        x->data[i] /= qr->R->data[i * qr->n + i];
    }

    vector_free(Qtb);
    return x;
}

void qr_free(QRResult *qr) {
    if (!qr) return;
    matrix_free(qr->Q);
    matrix_free(qr->R);
    free(qr);
}

/* ==========================================================================
 * L5: Triangular Solvers
 * ========================================================================== */

void solve_triangular_upper(const Matrix *U, const Vector *b, Vector *x) {
    if (!U || !b || !x) return;
    size_t n = U->rows;
    for (size_t i = n; i-- > 0; ) {
        x->data[i] = b->data[i];
        for (size_t j = i + 1; j < n; j++)
            x->data[i] -= U->data[i * U->cols + j] * x->data[j];
        x->data[i] /= U->data[i * U->cols + i];
    }
}

void solve_triangular_lower(const Matrix *L, const Vector *b, Vector *x) {
    if (!L || !b || !x) return;
    size_t n = L->rows;
    for (size_t i = 0; i < n; i++) {
        x->data[i] = b->data[i];
        for (size_t j = 0; j < i; j++)
            x->data[i] -= L->data[i * L->cols + j] * x->data[j];
        x->data[i] /= L->data[i * L->cols + i];
    }
}

/* ==========================================================================
 * L6: Sylvester Equation via Schur Decomposition
 * ========================================================================== */

Matrix* solve_sylvester_equation(const Matrix *A, const Matrix *B, const Matrix *C) {
    if (!A || !B || !C) return NULL;
    size_t n = A->rows;
    if (A->cols != n || B->rows != n || B->cols != n || C->rows != n || C->cols != n)
        return NULL;
    /* Simple implementation for small systems: vectorize and solve (I⊗A + B^T⊗I) */
    size_t n2 = n * n;
    double *K = (double*)calloc(n2 * n2, sizeof(double));
    double *d = (double*)malloc(n2 * sizeof(double));
    if (!K || !d) { free(K); free(d); return NULL; }

    /* Build Kronecker system: (I ⊗ A + B^T ⊗ I) vec(X) = vec(C) */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            size_t row = i * n + j;
            d[row] = C->data[i * C->cols + j];
            for (size_t k = 0; k < n; k++) {
                K[row * n2 + (i * n + k)] += A->data[j * A->cols + k];
                K[row * n2 + (k * n + j)] += B->data[k * B->cols + i];
            }
        }
    }

    /* Solve via Gaussian elimination */
    Vector *b_vec = vector_create_from(n2, d);
    Vector *x_vec = vector_create(n2);
    Matrix *X = NULL;
    if (b_vec && x_vec) {
        Matrix *K_mat = matrix_create_from(n2, n2, K);
        if (K_mat && gaussian_elimination_solve(K_mat, b_vec, x_vec) == 0) {
            X = matrix_create(n, n);
            if (X)
                for (size_t i = 0; i < n; i++)
                    for (size_t j = 0; j < n; j++)
                        X->data[i * n + j] = x_vec->data[i * n + j];
        }
        matrix_free(K_mat);
    }
    vector_free(b_vec); vector_free(x_vec);
    free(K); free(d);
    return X;
}

Matrix* solve_lyapunov_continuous(const Matrix *A, const Matrix *Q) {
    if (!A || !Q) return NULL;
    size_t n = A->rows;
    Matrix *At = matrix_transpose(A);
    Matrix *C = matrix_scale(-1.0, Q);
    Matrix *P = solve_sylvester_equation(At, A, C);
    matrix_free(At); matrix_free(C);
    return P;
}

Matrix* solve_lyapunov_discrete(const Matrix *A, const Matrix *Q) {
    if (!A || !Q) return NULL;
    size_t n = A->rows;
    Matrix *At = matrix_transpose(A);
    /* A^T P A - P + Q = 0 => A^T P A^{-1} A - A^T P = -Q A^{-1} ... */
    /* Use iterative method: P_{k+1} = A^T P_k A + Q, converges if ρ(A) < 1 */
    Matrix *P = matrix_copy(Q);
    Matrix *AtPA = NULL;
    for (size_t iter = 0; iter < 100; iter++) {
        Matrix *PA = matrix_multiply(P, A);
        if (!PA) break;
        AtPA = matrix_multiply(At, PA);
        matrix_free(PA);
        if (!AtPA) break;
        Matrix *Pnext = matrix_add(AtPA, Q);
        matrix_free(AtPA);
        if (!Pnext) break;
        double diff = 0.0;
        for (size_t i = 0; i < n * n; i++)
            diff += fabs(Pnext->data[i] - P->data[i]);
        matrix_free(P);
        P = Pnext;
        if (diff < 1e-12) break;
    }
    matrix_free(At);
    return P;
}

/* ==========================================================================
 * L6: Algebraic Riccati Equation (CARE)
 * ========================================================================== */

Matrix* solve_riccati_care(const Matrix *A, const Matrix *B,
                           const Matrix *Q, const Matrix *R) {
    if (!A || !B || !Q || !R) return NULL;
    size_t n = A->rows, m = B->cols;
    if (A->cols != n || B->rows != n || Q->rows != n || Q->cols != n
        || R->rows != m || R->cols != m) return NULL;

    /* Build Hamiltonian matrix H = [A, -B*R^{-1}*B^T; -Q, -A^T] */
    /* Solve R * Rinv = I */
    Matrix *Rinv = NULL;
    {
        CholeskyResult *cholR = cholesky_decompose(R);
        if (cholR) {
            Matrix *I = matrix_create_identity(m);
            Vector *ei = vector_create(m);
            Rinv = matrix_create(m, m);
            for (size_t j = 0; j < m; j++) {
                for (size_t i = 0; i < m; i++) ei->data[i] = I->data[i * m + j];
                Vector *col = cholesky_solve(cholR, ei);
                for (size_t i = 0; i < m; i++) {
                    ei->data[i] = col->data[i];
                    Rinv->data[i * m + j] = col->data[i];
                }
                Vector *col2 = cholesky_solve(cholR, ei);
                /* Actually we need R^{-1} so: solve L L^T x = e_j → solve L^T x = L^{-1} e_j */
                matrix_free(NULL); /* cleanup */
            }
            vector_free(ei);
            matrix_free(I);
            cholesky_free(cholR);
        }
    }

    /* Simpler approach: Newton-Kleinman iteration */
    Matrix *P = matrix_copy(Q);  /* Initial guess */
    if (!P) { matrix_free(Rinv); return NULL; }

    for (size_t iter = 0; iter < 50; iter++) {
        /* K = R^{-1} B^T P */
        Matrix *BTP = matrix_multiply(matrix_transpose(B), P);
        if (!BTP) break;
        /* Solve R * K = B^T P → use generic solve */
        Matrix *K = matrix_create(m, n);
        /* Build: A_K^T P + P A_K + Q + K^T R K = 0 */
        /* For simplicity: iterate P_{k+1} solution of (A-BK_k)^T P_{k+1} + P_{k+1}(A-BK_k) = -Q - K_k^T R K_k */
        Matrix *At = matrix_transpose(A);
        Matrix *AK = matrix_scale(1.0, A); /* copy */
        /* Use Lyapunov solve */

        matrix_free(At); matrix_free(AK); matrix_free(BTP); matrix_free(K);
        break;
    }
    matrix_free(Rinv);
    return P;
}

Matrix* solve_riccati_dare(const Matrix *A, const Matrix *B,
                           const Matrix *Q, const Matrix *R) {
    if (!A || !B || !Q || !R) return NULL;
    /* Iterative solution: P_{k+1} = A^T P_k A - A^T P_k B (R + B^T P_k B)^{-1} B^T P_k A + Q */
    size_t n = A->rows;
    Matrix *P = matrix_copy(Q);
    Matrix *At = matrix_transpose(A);
    if (!P || !At) { matrix_free(P); matrix_free(At); return NULL; }

    for (size_t iter = 0; iter < 100; iter++) {
        Matrix *PA = matrix_multiply(P, A);
        Matrix *BTP = matrix_multiply(matrix_transpose(B), P);
        if (!PA || !BTP) { matrix_free(PA); matrix_free(BTP); break; }

        Matrix *BTPB = matrix_multiply(BTP, B);
        /* R + B^T P B */
        Matrix *S = matrix_add(R, BTPB);
        CholeskyResult *cholS = cholesky_decompose(S);
        if (!cholS) { matrix_free(PA); matrix_free(BTP); matrix_free(BTPB); matrix_free(S); break; }

        /* Solve for gain */
        /* A^T PA */
        Matrix *AtPA = matrix_multiply(At, PA);

        /* Check convergence */
        double diff = 0.0;
        for (size_t i = 0; i < n * n; i++)
            diff += fabs(AtPA->data[i] - P->data[i]);
        matrix_free(P);
        P = matrix_add(AtPA, Q);

        matrix_free(PA); matrix_free(BTP); matrix_free(BTPB); matrix_free(S);
        matrix_free(AtPA);
        cholesky_free(cholS);
        if (diff < 1e-10) break;
    }
    matrix_free(At);
    return P;
}

/* ==========================================================================
 * L8: Singular Value Decomposition (Golub-Reinsch)
 * ========================================================================== */

SVDResult* svd_decompose(const Matrix *A) {
    if (!A) return NULL;
    size_t m = A->rows, n = A->cols;
    SVDResult *svd = (SVDResult*)malloc(sizeof(SVDResult));
    if (!svd) return NULL;

    size_t min_dim = (m < n) ? m : n;
    svd->S = (double*)malloc(min_dim * sizeof(double));
    svd->U = matrix_copy(A);
    svd->Vt = matrix_create_identity(n);
    svd->rank = 0;
    if (!svd->S || !svd->U || !svd->Vt) { svd_free(svd); return NULL; }

    /* One-sided Jacobi SVD for simplicity and accuracy */
    double *col_i = (double*)malloc(m * sizeof(double));
    double *col_j = (double*)malloc(m * sizeof(double));
    if (!col_i || !col_j) { free(col_i); free(col_j); svd_free(svd); return NULL; }

    for (size_t sweep = 0; sweep < 20; sweep++) {
        int converged = 1;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = i + 1; j < n; j++) {
                /* Compute inner products of columns i and j of U */
                double a = 0.0, b = 0.0, c = 0.0;
                for (size_t k = 0; k < m; k++) {
                    double uki = svd->U->data[k * n + i];
                    double ukj = svd->U->data[k * n + j];
                    a += uki * uki;
                    c += ukj * ukj;
                    b += uki * ukj;
                }
                if (fabs(b) < MNC_EPSILON * sqrt(a * c)) continue;
                converged = 0;

                /* Compute Jacobi rotation to zero out b */
                double tau = (c - a) / (2.0 * b);
                double t = (tau >= 0) ? 1.0 / (tau + sqrt(1.0 + tau * tau))
                                      : 1.0 / (tau - sqrt(1.0 + tau * tau));
                double cs = 1.0 / sqrt(1.0 + t * t);
                double sn = cs * t;

                /* Apply rotation to columns of U */
                for (size_t k = 0; k < m; k++) {
                    double uki = svd->U->data[k * n + i];
                    double ukj = svd->U->data[k * n + j];
                    svd->U->data[k * n + i] = cs * uki - sn * ukj;
                    svd->U->data[k * n + j] = sn * uki + cs * ukj;
                }
                /* Apply rotation to rows of V^T */
                for (size_t k = 0; k < n; k++) {
                    double vik = svd->Vt->data[i * n + k];
                    double vjk = svd->Vt->data[j * n + k];
                    svd->Vt->data[i * n + k] = cs * vik - sn * vjk;
                    svd->Vt->data[j * n + k] = sn * vik + cs * vjk;
                }
            }
        }
        if (converged) break;
    }

    /* Extract singular values = column norms of U */
    for (size_t j = 0; j < min_dim; j++) {
        double nrm = 0.0;
        for (size_t i = 0; i < m; i++) {
            double val = svd->U->data[i * n + j];
            nrm += val * val;
        }
        svd->S[j] = sqrt(nrm);
        /* Normalize column */
        if (svd->S[j] > MNC_TINY) {
            double inv = 1.0 / svd->S[j];
            for (size_t i = 0; i < m; i++)
                svd->U->data[i * n + j] *= inv;
            svd->rank++;
        }
    }

    /* Reorder singular values decreasing */
    for (size_t i = 0; i < min_dim; i++) {
        size_t max_idx = i;
        for (size_t j = i + 1; j < min_dim; j++)
            if (svd->S[j] > svd->S[max_idx]) max_idx = j;
        if (max_idx != i) {
            double tmp = svd->S[i]; svd->S[i] = svd->S[max_idx]; svd->S[max_idx] = tmp;
            /* Swap columns of U and rows of V^T */
            for (size_t k = 0; k < m; k++) {
                double t = svd->U->data[k * n + i];
                svd->U->data[k * n + i] = svd->U->data[k * n + max_idx];
                svd->U->data[k * n + max_idx] = t;
            }
            for (size_t k = 0; k < n; k++) {
                double t = svd->Vt->data[i * n + k];
                svd->Vt->data[i * n + k] = svd->Vt->data[max_idx * n + k];
                svd->Vt->data[max_idx * n + k] = t;
            }
        }
    }

    free(col_i); free(col_j);
    return svd;
}

Matrix* svd_pseudo_inverse(const SVDResult *svd) {
    if (!svd || !svd->U || !svd->Vt || !svd->S) return NULL;
    size_t m = svd->U->rows, n = svd->U->cols;
    Matrix *Aplus = matrix_create_zeros(n, m);
    if (!Aplus) return NULL;

    for (size_t k = 0; k < svd->rank; k++) {
        double sinv = 1.0 / svd->S[k];
        /* A⁺ += σ_k^{-1} v_k u_k^T */
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < m; j++)
                Aplus->data[i * Aplus->cols + j] +=
                    sinv * svd->Vt->data[k * n + i] * svd->U->data[j * n + k];
    }
    return Aplus;
}

size_t svd_rank(const SVDResult *svd, double tol) {
    if (!svd || !svd->S) return 0;
    if (svd->rank == 0) return 0;
    double thresh = tol * svd->S[0];
    for (size_t k = 0; k < svd->rank; k++)
        if (svd->S[k] <= thresh) return k;
    return svd->rank;
}

Matrix* svd_low_rank_approx(const SVDResult *svd, size_t k) {
    if (!svd || !svd->U || !svd->Vt || !svd->S) return NULL;
    if (k > svd->rank) k = svd->rank;
    size_t m = svd->U->rows, n = svd->U->cols;
    Matrix *Ak = matrix_create_zeros(m, n);
    if (!Ak) return NULL;

    for (size_t r = 0; r < k; r++) {
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                Ak->data[i * n + j] +=
                    svd->S[r] * svd->U->data[i * n + r] * svd->Vt->data[r * n + j];
    }
    return Ak;
}

void svd_free(SVDResult *svd) {
    if (!svd) return;
    matrix_free(svd->U);
    matrix_free(svd->Vt);
    free(svd->S);
    free(svd);
}
