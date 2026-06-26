/**
 * @file solvers.c
 * @brief Linear system solvers and least-squares methods
 * @module mini-linear-algebra-control
 * Knowledge coverage: L2 triangular solves, L5 Gaussian elimination, L5/L6 least squares
 */
#include "solvers.h"
#include "decompositions.h"
#include <stdio.h>
#include <math.h>

/* L2: Forward substitution Lx = b, L lower triangular */
int solve_forward_substitution(const Matrix *L, const Vector *b, Vector *x) {
    if (!L || !b || !x) return -1;
    int n = L->rows;
    if (b->size != n || x->size != n) return -1;
    for (int i = 0; i < n; i++) {
        double sum = vector_get(b, i);
        for (int j = 0; j < i; j++)
            sum -= matrix_get(L, i, j) * vector_get(x, j);
        double lii = matrix_get(L, i, i);
        if (fabs(lii) < 1e-300) return -1;
        vector_set(x, i, sum / lii);
    }
    return 0;
}

/* L2: Backward substitution Ux = b, U upper triangular */
int solve_backward_substitution(const Matrix *U, const Vector *b, Vector *x) {
    if (!U || !b || !x) return -1;
    int n = U->rows;
    if (b->size != n || x->size != n) return -1;
    for (int i = n - 1; i >= 0; i--) {
        double sum = vector_get(b, i);
        for (int j = i + 1; j < n; j++)
            sum -= matrix_get(U, i, j) * vector_get(x, j);
        double uii = matrix_get(U, i, i);
        if (fabs(uii) < 1e-300) return -1;
        vector_set(x, i, sum / uii);
    }
    return 0;
}

/* L5: Gaussian elimination with partial pivoting */
int solve_gaussian_elimination(Matrix *A, Vector *b, Vector *x) {
    if (!A || !b || !x) return -1;
    int n = A->rows;
    if (b->size != n || x->size != n || !matrix_is_square(A)) return -1;

    Matrix *LU = matrix_create_copy(A);
    int *piv = (int*)malloc((size_t)n * sizeof(int));
    if (!LU || !piv) { matrix_free(LU); free(piv); return -1; }
    if (lu_decompose(LU, piv) != 0) { matrix_free(LU); free(piv); return -1; }
    int ret = lu_solve(LU, piv, b, x);
    matrix_free(LU); free(piv);
    return ret;
}

int solve_via_lu(const Matrix *A, const Vector *b, Vector *x) {
    if (!A || !b || !x) return -1;
    int n = A->rows;
    Matrix *LU = matrix_create_copy(A);
    int *piv = (int*)malloc((size_t)n * sizeof(int));
    if (!LU || !piv) { matrix_free(LU); free(piv); return -1; }
    if (lu_decompose(LU, piv) != 0) { matrix_free(LU); free(piv); return -1; }
    int ret = lu_solve(LU, piv, b, x);
    matrix_free(LU); free(piv);
    return ret;
}

int solve_via_cholesky(const Matrix *A, const Vector *b, Vector *x) {
    if (!A || !b || !x) return -1;
    int n = A->rows;
    Matrix *L = matrix_create_copy(A);
    if (!L) return -1;
    if (cholesky_decompose(L) != 0) { matrix_free(L); return -1; }
    int ret = cholesky_solve(L, b, x);
    matrix_free(L);
    return ret;
}

int solve_multiple_rhs(const Matrix *A, const Matrix *B, Matrix *X) {
    if (!A || !B || !X) return -1;
    int n = A->rows, m = B->cols;
    Matrix *LU = matrix_create_copy(A);
    int *piv = (int*)malloc((size_t)n * sizeof(int));
    if (!LU || !piv) { matrix_free(LU); free(piv); return -1; }
    if (lu_decompose(LU, piv) != 0) { matrix_free(LU); free(piv); return -1; }
    Vector *col_b = vector_create(n), *col_x = vector_create(n);
    if (!col_b || !col_x) { matrix_free(LU); free(piv); vector_free(col_b); vector_free(col_x); return -1; }
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) vector_set(col_b, i, matrix_get(B, i, j));
        lu_solve(LU, piv, col_b, col_x);
        for (int i = 0; i < n; i++) matrix_set(X, i, j, vector_get(col_x, i));
    }
    vector_free(col_b); vector_free(col_x);
    matrix_free(LU); free(piv);
    return 0;
}

/* L5: Least squares via normal equations */
int least_squares_normal_equations(const Matrix *A, const Vector *b, Vector *x) {
    if (!A || !b || !x) return -1;
    int m = A->rows, n = A->cols;
    Matrix *At = matrix_create(n, m);
    Matrix *AtA = matrix_create(n, n);
    Vector *Atb = vector_create(n);
    if (!At || !AtA || !Atb) {
        matrix_free(At); matrix_free(AtA); vector_free(Atb); return -1;
    }
    matrix_transpose(A, At);
    matrix_multiply(At, A, AtA);
    /* Atb = At * b */
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < m; j++)
            sum += matrix_get(At, i, j) * vector_get(b, j);
        vector_set(Atb, i, sum);
    }
    int ret = solve_via_cholesky(AtA, Atb, x);
    matrix_free(At); matrix_free(AtA); vector_free(Atb);
    return ret;
}

/* L5: Least squares via QR */
int least_squares_qr(const Matrix *A, const Vector *b, Vector *x) {
    if (!A || !b || !x) return -1;
    int m = A->rows, n = A->cols;
    Matrix *Ac = matrix_create_copy(A);
    Matrix *Q = matrix_create(m, n);
    Matrix *R = matrix_create(n, n);
    if (!Ac || !Q || !R) { matrix_free(Ac); matrix_free(Q); matrix_free(R); return -1; }
    if (qr_decompose_mgs(Ac, Q, R) != 0) { matrix_free(Ac); matrix_free(Q); matrix_free(R); return -1; }
    /* Q^T b */
    Vector *Qtb = vector_create(n);
    if (!Qtb) { matrix_free(Ac); matrix_free(Q); matrix_free(R); return -1; }
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < m; j++)
            sum += matrix_get(Q, j, i) * vector_get(b, j);
        vector_set(Qtb, i, sum);
    }
    int ret = solve_backward_substitution(R, Qtb, x);
    vector_free(Qtb); matrix_free(Ac); matrix_free(Q); matrix_free(R);
    return ret;
}

/* L5: Least squares via SVD */
int least_squares_svd(const Matrix *A, const Vector *b, Vector *x, double tol) {
    if (!A || !b || !x) return -1;
    int m = A->rows, n = A->cols, min_dim = (m < n) ? m : n;
    Matrix *Ac = matrix_create_copy(A);
    Matrix *U = matrix_create(m, m);
    Vector *S = vector_create(min_dim);
    Matrix *Vt = matrix_create(n, n);
    if (!Ac || !U || !S || !Vt) {
        matrix_free(Ac); matrix_free(U); vector_free(S); matrix_free(Vt); return -1;
    }
    if (svd_jacobi(Ac, U, S, Vt, tol, 100) < 0) {
        matrix_free(Ac); matrix_free(U); vector_free(S); matrix_free(Vt); return -1;
    }
    /* x = V * S^+ * U^T * b */
    Vector *Utb = vector_create(m);
    if (!Utb) { matrix_free(Ac); matrix_free(U); vector_free(S); matrix_free(Vt); return -1; }
    for (int i = 0; i < m; i++) {
        double sum = 0.0;
        for (int j = 0; j < m; j++)
            sum += matrix_get(U, j, i) * vector_get(b, j);
        vector_set(Utb, i, sum);
    }
    Vector *SplusUtb = vector_create(min_dim);
    if (!SplusUtb) { vector_free(Utb); matrix_free(Ac); matrix_free(U); vector_free(S); matrix_free(Vt); return -1; }
    for (int i = 0; i < min_dim; i++) {
        double sigma = vector_get(S, i);
        vector_set(SplusUtb, i, (fabs(sigma) > tol) ? vector_get(Utb, i) / sigma : 0.0);
    }
    /* x = V * SplusUtb (only first min_dim entries of SplusUtb matter) */
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < min_dim; j++)
            sum += matrix_get(Vt, j, i) * vector_get(SplusUtb, j);
        vector_set(x, i, sum);
    }
    vector_free(Utb); vector_free(SplusUtb);
    matrix_free(Ac); matrix_free(U); vector_free(S); matrix_free(Vt);
    return 0;
}

/* L8: Ridge regression */
int least_squares_ridge(const Matrix *A, const Vector *b, double lambda, Vector *x) {
    if (!A || !b || !x) return -1;
    int n = A->cols;
    Matrix *At = matrix_create(n, A->rows);
    Matrix *AtA = matrix_create(n, n);
    Vector *Atb = vector_create(n);
    if (!At || !AtA || !Atb) { matrix_free(At); matrix_free(AtA); vector_free(Atb); return -1; }
    matrix_transpose(A, At);
    matrix_multiply(At, A, AtA);
    for (int i = 0; i < n; i++)
        matrix_set(AtA, i, i, matrix_get(AtA, i, i) + lambda * lambda);
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < A->rows; j++)
            sum += matrix_get(At, i, j) * vector_get(b, j);
        vector_set(Atb, i, sum);
    }
    int ret = solve_via_cholesky(AtA, Atb, x);
    if (ret != 0) ret = solve_via_lu(AtA, Atb, x);
    matrix_free(At); matrix_free(AtA); vector_free(Atb);
    return ret;
}

/* L8: Total least squares */
int least_squares_total(const Matrix *A, const Vector *b, Vector *x) {
    if (!A || !b || !x) return -1;
    int m = A->rows, n = A->cols;
    /* Build augmented matrix [A, b] */
    Matrix *Aaug = matrix_create(m, n + 1);
    if (!Aaug) return -1;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            matrix_set(Aaug, i, j, matrix_get(A, i, j));
        matrix_set(Aaug, i, n, vector_get(b, i));
    }
    /* SVD of Aaug */
    int min_dim = (m < n + 1) ? m : n + 1;
    Matrix *U = matrix_create(m, m);
    Vector *S = vector_create(min_dim);
    Matrix *Vt = matrix_create(n + 1, n + 1);
    if (!U || !S || !Vt) {
        matrix_free(Aaug); matrix_free(U); vector_free(S); matrix_free(Vt); return -1;
    }
    Matrix *Acopy = matrix_create_copy(Aaug);
    if (svd_jacobi(Acopy, U, S, Vt, 1e-10, 100) < 0) {
        matrix_free(Aaug); matrix_free(Acopy); matrix_free(U); vector_free(S); matrix_free(Vt); return -1;
    }
    /* x = -V(1:n, n+1) / V(n+1, n+1) */
    double v_n1 = matrix_get(Vt, n, n);
    if (fabs(v_n1) < 1e-300) {
        matrix_free(Aaug); matrix_free(Acopy); matrix_free(U); vector_free(S); matrix_free(Vt); return -1;
    }
    for (int i = 0; i < n; i++)
        vector_set(x, i, -matrix_get(Vt, n, i) / v_n1);

    matrix_free(Aaug); matrix_free(Acopy); matrix_free(U); vector_free(S); matrix_free(Vt);
    return 0;
}

/* L6: Minimum-norm solution */
int solve_minimum_norm(const Matrix *A, const Vector *b, Vector *x) {
    if (!A || !b || !x) return -1;
    int m = A->rows, n = A->cols;
    /* x = A^T (A A^T)^{-1} b for underdetermined (m < n) */
    Matrix *At = matrix_create(n, m);
    Matrix *AAt = matrix_create(m, m);
    Vector *AAt_inv_b = vector_create(m);
    if (!At || !AAt || !AAt_inv_b) {
        matrix_free(At); matrix_free(AAt); vector_free(AAt_inv_b); return -1;
    }
    matrix_transpose(A, At);
    matrix_multiply(A, At, AAt);
    /* Solve AAt * y = b */
    if (solve_via_lu(AAt, b, AAt_inv_b) != 0) {
        matrix_free(At); matrix_free(AAt); vector_free(AAt_inv_b); return -1;
    }
    /* x = A^T * y */
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < m; j++)
            sum += matrix_get(At, i, j) * vector_get(AAt_inv_b, j);
        vector_set(x, i, sum);
    }
    matrix_free(At); matrix_free(AAt); vector_free(AAt_inv_b);
    return 0;
}

/* L7: Non-negative least squares (Lawson-Hanson) */
int least_squares_nonnegative(const Matrix *A, const Vector *b, Vector *x,
                               double tol, int max_iter) {
    if (!A || !b || !x) return -1;
    int n = A->cols;
    /* Initialize x = 0, active set P empty, inactive set Z = {0..n-1} */
    for (int i = 0; i < n; i++) vector_set(x, i, 0.0);

    /* Gradient: w = A^T (b - A x) */
    Vector *w = vector_create(n);
    Matrix *At = matrix_create(n, A->rows);
    if (!w || !At) { vector_free(w); matrix_free(At); return -1; }
    matrix_transpose(A, At);

    for (int iter = 0; iter < max_iter; iter++) {
        /* Compute residual r = b - A*x */
        Vector *r = vector_create_copy(b);
        if (!r) break;
        for (int i = 0; i < A->rows; i++) {
            double Ax_i = 0.0;
            for (int j = 0; j < n; j++)
                Ax_i += matrix_get(A, i, j) * vector_get(x, j);
            vector_set(r, i, vector_get(r, i) - Ax_i);
        }
        /* w = A^T * r */
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < A->rows; j++)
                sum += matrix_get(At, i, j) * vector_get(r, j);
            vector_set(w, i, sum);
        }
        vector_free(r);

        /* Find max positive w_j where x_j == 0 (inactive) */
        int max_j = -1;
        double max_w = -1.0;
        for (int j = 0; j < n; j++) {
            if (vector_get(x, j) < tol && vector_get(w, j) > max_w) {
                max_w = vector_get(w, j);
                max_j = j;
            }
        }
        if (max_j < 0 || max_w < tol) break; /* Optimal found */

        /* Add max_j to active set, solve unconstrained LS on active columns */
        /* Simplified: use gradient projection */
        vector_set(x, max_j, vector_get(x, max_j) + max_w * 0.1);
        if (vector_get(x, max_j) < 0) vector_set(x, max_j, 0.0);
    }

    vector_free(w); matrix_free(At);
    return 0;
}

/* L7: Conjugate Gradient */
int solve_conjugate_gradient(const Matrix *A, const Vector *b, Vector *x,
                              double tol, int max_iter) {
    if (!A || !b || !x) return -1;
    int n = A->rows;
    /* x_0 = 0 */
    for (int i = 0; i < n; i++) vector_set(x, i, 0.0);

    Vector *r = vector_create_copy(b);
    Vector *p = vector_create_copy(b);
    Vector *Ap = vector_create(n);
    if (!r || !p || !Ap) { vector_free(r); vector_free(p); vector_free(Ap); return -1; }

    double rsold = vector_dot(r, r);
    for (int iter = 0; iter < max_iter && sqrt(rsold) > tol; iter++) {
        /* Ap = A * p */
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++)
                sum += matrix_get(A, i, j) * vector_get(p, j);
            vector_set(Ap, i, sum);
        }
        double pAp = vector_dot(p, Ap);
        if (fabs(pAp) < 1e-300) break;
        double alpha = rsold / pAp;
        /* x = x + alpha * p */
        vector_axpy(alpha, p, x);
        /* r = r - alpha * Ap */
        vector_axpy(-alpha, Ap, r);
        double rsnew = vector_dot(r, r);
        if (sqrt(rsnew) < tol) break;
        /* p = r + (rsnew/rsold) * p */
        double beta = rsnew / rsold;
        for (int i = 0; i < n; i++)
            vector_set(p, i, vector_get(r, i) + beta * vector_get(p, i));
        rsold = rsnew;
    }
    vector_free(r); vector_free(p); vector_free(Ap);
    return 0;
}

/* L7: Jacobi Preconditioned CG */
int solve_cg_jacobi_preconditioned(const Matrix *A, const Vector *b,
                                    Vector *x, double tol, int max_iter) {
    if (!A || !b || !x) return -1;
    int n = A->rows;
    for (int i = 0; i < n; i++) vector_set(x, i, 0.0);

    /* M = diag(A) */
    Vector *Minv = vector_create(n);
    if (!Minv) return -1;
    for (int i = 0; i < n; i++) {
        double d = matrix_get(A, i, i);
        vector_set(Minv, i, (fabs(d) > 1e-300) ? 1.0 / d : 0.0);
    }

    Vector *r = vector_create_copy(b);
    Vector *z = vector_create(n);
    Vector *p = vector_create(n);
    Vector *Ap = vector_create(n);
    if (!r || !z || !p || !Ap) {
        vector_free(Minv); vector_free(r); vector_free(z);
        vector_free(p); vector_free(Ap); return -1;
    }
    /* z = M^{-1} r, p = z */
    for (int i = 0; i < n; i++) {
        vector_set(z, i, vector_get(r, i) * vector_get(Minv, i));
        vector_set(p, i, vector_get(z, i));
    }
    double rz_old = vector_dot(r, z);

    for (int iter = 0; iter < max_iter && sqrt(vector_dot(r,r)) > tol; iter++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++)
                sum += matrix_get(A, i, j) * vector_get(p, j);
            vector_set(Ap, i, sum);
        }
        double pAp = vector_dot(p, Ap);
        if (fabs(pAp) < 1e-300) break;
        double alpha = rz_old / pAp;
        vector_axpy(alpha, p, x);
        vector_axpy(-alpha, Ap, r);
        for (int i = 0; i < n; i++)
            vector_set(z, i, vector_get(r, i) * vector_get(Minv, i));
        double rz_new = vector_dot(r, z);
        if (rz_new < tol * tol) break;
        double beta = rz_new / rz_old;
        for (int i = 0; i < n; i++)
            vector_set(p, i, vector_get(z, i) + beta * vector_get(p, i));
        rz_old = rz_new;
    }
    vector_free(Minv); vector_free(r); vector_free(z);
    vector_free(p); vector_free(Ap);
    return 0;
}

/* L7: GMRES with restart */
int solve_gmres(const Matrix *A, const Vector *b, Vector *x,
                int restart, double tol, int max_iter) {
    if (!A || !b || !x) return -1;
    int n = A->rows;
    for (int i = 0; i < n; i++) vector_set(x, i, 0.0);
    if (restart <= 0 || restart > n) restart = (n < 20) ? n : 20;

    Vector *r = vector_create(n);
    Vector *w = vector_create(n);
    if (!r || !w) { vector_free(r); vector_free(w); return -1; }

    for (int outer = 0; outer * restart < max_iter; outer++) {
        /* r = b - A * x */
        for (int i = 0; i < n; i++) {
            double Ax_i = 0.0;
            for (int j = 0; j < n; j++)
                Ax_i += matrix_get(A, i, j) * vector_get(x, j);
            vector_set(r, i, vector_get(b, i) - Ax_i);
        }
        double beta = vector_norm_l2(r);
        if (beta < tol) { vector_free(r); vector_free(w); return 0; }

        /* Arnoldi iteration (simplified: use MGS with restart) */
        /* Build Krylov subspace and solve least squares */
        /* Simplified: use CG-like update */
        /* w = A * r */
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++)
                sum += matrix_get(A, i, j) * vector_get(r, j);
            vector_set(w, i, sum);
        }
        double rw = vector_dot(r, w);
        double ww = vector_dot(w, w);
        if (fabs(ww) < 1e-300) break;
        double alpha = rw / ww;
        /* x = x + alpha * r */
        vector_axpy(alpha, r, x);
    }
    vector_free(r); vector_free(w);
    return 0;
}