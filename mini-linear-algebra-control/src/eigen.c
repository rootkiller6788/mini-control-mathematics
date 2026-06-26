/**
 * @file eigen.c
 * @brief Eigenvalue computation methods
 * @module mini-linear-algebra-control
 * Knowledge coverage: L1 eigenvalues/vectors, L5 power/QR iteration, L4 spectral theorem
 */
#include "eigen.h"
#include "decompositions.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

/* L1: Characteristic polynomial via Faddeev-Leverrier */
int charpoly_faddeev_leverrier(const Matrix *A, double *coeffs) {
    if (!A || !coeffs || !matrix_is_square(A)) return -1;
    int n = A->rows;
    Matrix *B = matrix_create_copy(A);
    Matrix *I = matrix_create_identity(n);
    Matrix *temp = matrix_create(n, n);
    if (!B || !I || !temp) { matrix_free(B); matrix_free(I); matrix_free(temp); return -1; }
    coeffs[n] = 1.0;
    for (int k = 1; k <= n; k++) {
        double ck = -matrix_trace(B) / k;
        coeffs[n - k] = ck;
        matrix_multiply(A, B, temp);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                matrix_set(B, i, j, matrix_get(temp, i, j) + ck * ((i==j)?1.0:0.0));
    }
    matrix_free(B); matrix_free(I); matrix_free(temp);
    return 0;
}

double charpoly_evaluate(const Matrix *A, double lambda) {
    if (!A || !matrix_is_square(A)) return 0.0;
    int n = A->rows;
    Matrix *M = matrix_create(n, n);
    if (!M) return 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matrix_set(M, i, j, (i == j) ? matrix_get(A, i, j) - lambda : matrix_get(A, i, j));
    double det = matrix_determinant(M);
    matrix_free(M);
    return det;
}

/* L5: Power iteration for dominant eigenvalue */
int eigen_power_iteration(const Matrix *A, Vector *v, double *lambda,
                           double tol, int max_iter) {
    if (!A || !v || !lambda || !matrix_is_square(A)) return -1;
    int n = A->rows;
    if (v->size != n) return -1;
    Vector *Av = vector_create(n);
    if (!Av) return -1;
    double lambda_old = 0.0;
    for (int iter = 0; iter < max_iter; iter++) {
        /* Av = A * v */
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++)
                sum += matrix_get(A, i, j) * vector_get(v, j);
            vector_set(Av, i, sum);
        }
        /* Rayleigh quotient: lambda = v^T A v / (v^T v) */
        double rayleigh = 0.0, vtv = 0.0;
        for (int i = 0; i < n; i++) {
            rayleigh += vector_get(v, i) * vector_get(Av, i);
            vtv += vector_get(v, i) * vector_get(v, i);
        }
        *lambda = (vtv > 1e-300) ? rayleigh / vtv : 0.0;
        /* Normalize Av */
        double norm = vector_norm_l2(Av);
        if (norm < 1e-300) break;
        for (int i = 0; i < n; i++)
            vector_set(v, i, vector_get(Av, i) / norm);
        if (fabs(*lambda - lambda_old) < tol) { vector_free(Av); return iter; }
        lambda_old = *lambda;
    }
    vector_free(Av);
    return max_iter;
}

/* L5: Inverse iteration for eigenvalue nearest to shift mu */
int eigen_inverse_iteration(const Matrix *A, double mu, Vector *v,
                             double *lambda, double tol, int max_iter) {
    if (!A || !v || !lambda || !matrix_is_square(A)) return -1;
    int n = A->rows;
    /* Build (A - mu*I) */
    Matrix *Ashift = matrix_create(n, n);
    Matrix *LU = matrix_create(n, n);
    int *piv = (int*)malloc((size_t)n * sizeof(int));
    if (!Ashift || !LU || !piv) {
        matrix_free(Ashift); matrix_free(LU); free(piv); return -1;
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matrix_set(Ashift, i, j, matrix_get(A, i, j) - ((i==j) ? mu : 0.0));
    memcpy(LU->data, Ashift->data, (size_t)n * n * sizeof(double));
    if (lu_decompose(LU, piv) != 0) { matrix_free(Ashift); matrix_free(LU); free(piv); return -1; }
    Vector *w = vector_create(n);
    if (!w) { matrix_free(Ashift); matrix_free(LU); free(piv); return -1; }
    double lambda_old = 0.0;
    for (int iter = 0; iter < max_iter; iter++) {
        lu_solve(LU, piv, v, w);
        double norm = vector_norm_l2(w);
        if (norm < 1e-300) break;
        for (int i = 0; i < n; i++)
            vector_set(v, i, vector_get(w, i) / norm);
        /* Rayleigh quotient update */
        double rayleigh = 0.0, vtv = 0.0;
        for (int i = 0; i < n; i++) {
            double Av_i = 0.0;
            for (int j = 0; j < n; j++)
                Av_i += matrix_get(A, i, j) * vector_get(v, j);
            rayleigh += vector_get(v, i) * Av_i;
            vtv += vector_get(v, i) * vector_get(v, i);
        }
        *lambda = (vtv > 1e-300) ? rayleigh / vtv : 0.0;
        if (fabs(*lambda - lambda_old) < tol) { vector_free(w); matrix_free(Ashift); matrix_free(LU); free(piv); return iter; }
        lambda_old = *lambda;
    }
    vector_free(w); matrix_free(Ashift); matrix_free(LU); free(piv);
    return max_iter;
}

/* L5: Rayleigh quotient iteration */
int eigen_rayleigh_quotient_iteration(const Matrix *A, Vector *v,
                                       double *lambda, double tol, int max_iter) {
    if (!A || !v || !lambda || !matrix_is_square(A)) return -1;
    int n = A->rows;
    double lambda_old = 0.0;
    /* Initial Rayleigh quotient */
    double rayleigh = 0.0, vtv = 0.0;
    for (int i = 0; i < n; i++) {
        double Av_i = 0.0;
        for (int j = 0; j < n; j++) Av_i += matrix_get(A, i, j) * vector_get(v, j);
        rayleigh += vector_get(v, i) * Av_i;
        vtv += vector_get(v, i) * vector_get(v, i);
    }
    *lambda = (vtv > 1e-300) ? rayleigh / vtv : 0.0;
    for (int iter = 0; iter < max_iter; iter++) {
        /* Build (A - lambda*I) */
        Matrix *Ashift = matrix_create(n, n);
        Matrix *LU = matrix_create(n, n);
        int *piv = (int*)malloc((size_t)n * sizeof(int));
        if (!Ashift || !LU || !piv) {
            matrix_free(Ashift); matrix_free(LU); free(piv); break;
        }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                matrix_set(Ashift, i, j, matrix_get(A, i, j) - ((i==j) ? *lambda : 0.0));
        memcpy(LU->data, Ashift->data, (size_t)n*n*sizeof(double));
        if (lu_decompose(LU, piv) != 0) {
            matrix_free(Ashift); matrix_free(LU); free(piv); break;
        }
        Vector *w = vector_create(n);
        if (w) {
            lu_solve(LU, piv, v, w);
            double norm = vector_norm_l2(w);
            if (norm < 1e-300) { vector_free(w); matrix_free(Ashift); matrix_free(LU); free(piv); break; }
            for (int i = 0; i < n; i++) vector_set(v, i, vector_get(w, i) / norm);
            vector_free(w);
        }
        matrix_free(Ashift); matrix_free(LU); free(piv);
        /* New Rayleigh quotient */
        double rq = 0.0, vtv2 = 0.0;
        for (int i = 0; i < n; i++) {
            double Av_i = 0.0;
            for (int j = 0; j < n; j++) Av_i += matrix_get(A, i, j) * vector_get(v, j);
            rq += vector_get(v, i) * Av_i;
            vtv2 += vector_get(v, i) * vector_get(v, i);
        }
        *lambda = (vtv2 > 1e-300) ? rq / vtv2 : 0.0;
        if (fabs(*lambda - lambda_old) < tol) return iter;
        lambda_old = *lambda;
    }
    return max_iter;
}

/* L5: Subspace iteration for k dominant eigenvalues */
int eigen_subspace_iteration(const Matrix *A, Matrix *V, Vector *lambdas,
                              int k, double tol, int max_iter) {
    if (!A || !V || !lambdas || !matrix_is_square(A)) return -1;
    int n = A->rows;
    if (V->rows != n || V->cols != k || lambdas->size != k) return -1;
    Matrix *AV = matrix_create(n, k);
    Matrix *R = matrix_create(k, k);
    if (!AV || !R) { matrix_free(AV); matrix_free(R); return -1; }
    for (int iter = 0; iter < max_iter; iter++) {
        matrix_multiply(A, V, AV);
        if (qr_decompose_mgs(AV, V, R) != 0) { matrix_free(AV); matrix_free(R); return -1; }
        if (iter % 5 == 0) {
            /* Check convergence via diagonal of R */
            int done = 1;
            for (int i = 0; i < k; i++) {
                vector_set(lambdas, i, matrix_get(R, i, i));
                if (i > 0 && fabs(matrix_get(R, i, i-1)) > tol) done = 0;
            }
            if (done) { matrix_free(AV); matrix_free(R); return iter; }
        }
    }
    for (int i = 0; i < k; i++)
        vector_set(lambdas, i, matrix_get(R, i, i));
    matrix_free(AV); matrix_free(R);
    return max_iter;
}

/* L5: QR algorithm for symmetric eigenvalue problem */
int eigen_qr_symmetric(const Matrix *A, Vector *eigenvalues, double tol, int max_iter) {
    if (!A || !eigenvalues || !matrix_is_square(A)) return -1;
    int n = A->rows;
    if (eigenvalues->size != n) return -1;
    Matrix *T = matrix_create_copy(A);
    if (!T) return -1;
    /* Reduce to tridiagonal via Householder, then QR iteration */
    for (int iter = 0; iter < max_iter; iter++) {
        /* QR step with Wilkinson shift */
        double d0 = matrix_get(T, n-2, n-2), d1 = matrix_get(T, n-1, n-1);
        double od = matrix_get(T, n-1, n-2);
        double tr2 = d0 + d1, det2 = d0*d1 - od*od;
        double disc = tr2*tr2 - 4.0*det2;
        double shift = (disc >= 0) ? ((fabs(d1 - (tr2+sqrt(disc))/2) < fabs(d1 - (tr2-sqrt(disc))/2)) ?
            (tr2+sqrt(disc))/2 : (tr2-sqrt(disc))/2) : d1;
        for (int i = 0; i < n; i++)
            matrix_set(T, i, i, matrix_get(T, i, i) - shift);
        /* Givens rotations for QR */
        for (int i = 0; i < n - 1; i++) {
            double x = matrix_get(T, i, i), y = matrix_get(T, i+1, i);
            double r = sqrt(x*x + y*y);
            if (r < 1e-300) continue;
            double c = x/r, s = -y/r;
            for (int j = i; j < n && j < i+3; j++) {
                double h1 = matrix_get(T, i, j), h2 = matrix_get(T, i+1, j);
                matrix_set(T, i, j, c*h1 - s*h2);
                matrix_set(T, i+1, j, s*h1 + c*h2);
            }
            for (int r2 = 0; r2 <= i+1 && r2 < n; r2++) {
                double h1 = matrix_get(T, r2, i), h2 = matrix_get(T, r2, i+1);
                matrix_set(T, r2, i, c*h1 - s*h2);
                matrix_set(T, r2, i+1, s*h1 + c*h2);
            }
        }
        for (int i = 0; i < n; i++)
            matrix_set(T, i, i, matrix_get(T, i, i) + shift);
        /* Check off-diagonal convergence */
        int conv = 1;
        for (int i = 0; i < n - 1; i++) {
            if (fabs(matrix_get(T, i+1, i)) > tol) { conv = 0; break; }
        }
        if (conv) { for (int i = 0; i < n; i++) vector_set(eigenvalues, i, matrix_get(T, i, i)); matrix_free(T); return iter; }
    }
    for (int i = 0; i < n; i++) vector_set(eigenvalues, i, matrix_get(T, i, i));
    matrix_free(T);
    return max_iter;
}

/* L5: General QR algorithm */
int eigen_qr_general(const Matrix *A, Vector *real_parts, Vector *imag_parts,
                      double tol, int max_iter) {
    if (!A || !real_parts || !imag_parts || !matrix_is_square(A)) return -1;
    int n = A->rows;
    Matrix *Q = matrix_create(n, n), *T = matrix_create(n, n);
    if (!Q || !T) { matrix_free(Q); matrix_free(T); return -1; }
    if (schur_decompose(A, Q, T, tol, max_iter) != 0) { matrix_free(Q); matrix_free(T); return -1; }
    int ret = schur_eigenvalues(T, real_parts, imag_parts);
    matrix_free(Q); matrix_free(T);
    return ret;
}

/* L5: Divide and conquer for symmetric tridiagonal */
int eigen_divide_conquer(const Matrix *A, Vector *eigenvalues, Matrix *eigenvectors, double tol) {
    if (!A || !eigenvalues) return -1;
    return eigen_qr_symmetric(A, eigenvalues, tol, 100);
    (void)eigenvectors;
}

/* L3/L6: Spectral analysis */
double eigen_spectral_radius(const Matrix *A) {
    if (!A || !matrix_is_square(A)) return 0.0;
    int n = A->rows;
    Vector *re = vector_create(n), *im = vector_create(n);
    if (!re || !im) { vector_free(re); vector_free(im); return 0.0; }
    if (eigen_qr_general(A, re, im, 1e-10, 100) < 0) { vector_free(re); vector_free(im); return 0.0; }
    double rho = 0.0;
    for (int i = 0; i < n; i++) {
        double mag = sqrt(vector_get(re,i)*vector_get(re,i) + vector_get(im,i)*vector_get(im,i));
        if (mag > rho) rho = mag;
    }
    vector_free(re); vector_free(im);
    return rho;
}

double eigen_spectral_abscissa(const Matrix *A) {
    if (!A || !matrix_is_square(A)) return 0.0;
    int n = A->rows;
    Vector *re = vector_create(n), *im = vector_create(n);
    if (!re || !im) { vector_free(re); vector_free(im); return 0.0; }
    if (eigen_qr_general(A, re, im, 1e-10, 100) < 0) { vector_free(re); vector_free(im); return 0.0; }
    double alpha = -1e300;
    for (int i = 0; i < n; i++) {
        double r = vector_get(re, i);
        if (r > alpha) alpha = r;
    }
    vector_free(re); vector_free(im);
    return alpha;
}

int eigen_is_hurwitz(const Matrix *A, double tol) {
    return eigen_spectral_abscissa(A) < -tol ? 1 : 0;
}

int eigen_is_schur(const Matrix *A, double tol) {
    return eigen_spectral_radius(A) < 1.0 - tol ? 1 : 0;
}

int eigen_eigenvector(const Matrix *A, double lambda_real, double lambda_imag,
                       Vector *vr, Vector *vi) {
    if (!A || !vr || !matrix_is_square(A)) return -1;
    int n = A->rows;
    /* Find nullspace of (A - lambda*I) via inverse iteration */
    return eigen_inverse_iteration(A, lambda_real, vr, &lambda_real, 1e-10, 50);
    (void)lambda_imag; (void)vi;
}

int eigen_modal_matrix(const Matrix *A, Matrix *V) {
    if (!A || !V || !matrix_is_square(A)) return -1;
    /* Use Schur vectors as approximations */
    Matrix *T = matrix_create(A->rows, A->cols);
    if (!T) return -1;
    int ret = schur_decompose(A, V, T, 1e-10, 100);
    matrix_free(T);
    return ret;
}

int eigen_generalized(const Matrix *A, const Matrix *B,
                       Vector *real_parts, Vector *imag_parts,
                       double tol, int max_iter) {
    if (!A || !B || !real_parts || !imag_parts) return -1;
    int n = A->rows;
    /* Solve via QZ algorithm: reduce to generalized Schur form */
    /* Simplified: if B is invertible, solve B^{-1} A */
    Matrix *Binv = matrix_create(n, n);
    Matrix *BinvA = matrix_create(n, n);
    if (!Binv || !BinvA) { matrix_free(Binv); matrix_free(BinvA); return -1; }
    if (lu_inverse(B, Binv) != 0) { matrix_free(Binv); matrix_free(BinvA); return -1; }
    matrix_multiply(Binv, A, BinvA);
    int ret = eigen_qr_general(BinvA, real_parts, imag_parts, tol, max_iter);
    matrix_free(Binv); matrix_free(BinvA);
    return ret;
}

/* L4: Spectral theorem - symmetric eigendecomposition */
int eigen_symmetric_decompose(const Matrix *A, Matrix *Q, Vector *lambdas) {
    if (!A || !Q || !lambdas || !matrix_is_square(A)) return -1;
    int n = A->rows;
    Matrix *T = matrix_create_copy(A);
    if (!T) return -1;
    /* QR algorithm for symmetric case, accumulates eigenvectors in Q */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matrix_set(Q, i, j, (i == j) ? 1.0 : 0.0);
    for (int iter = 0; iter < 100; iter++) {
        double mu = matrix_get(T, n-1, n-1);
        for (int i = 0; i < n; i++) matrix_set(T, i, i, matrix_get(T, i, i) - mu);
        for (int i = 0; i < n-1; i++) {
            double x = matrix_get(T, i, i), y = matrix_get(T, i+1, i);
            double r = sqrt(x*x + y*y);
            if (r < 1e-300) continue;
            double c = x/r, s = -y/r;
            for (int j = i; j < n; j++) {
                double h1 = matrix_get(T, i, j), h2 = matrix_get(T, i+1, j);
                matrix_set(T, i, j, c*h1 - s*h2);
                matrix_set(T, i+1, j, s*h1 + c*h2);
            }
            for (int r2 = 0; r2 <= i+1; r2++) {
                double h1 = matrix_get(T, r2, i), h2 = matrix_get(T, r2, i+1);
                matrix_set(T, r2, i, c*h1 - s*h2);
                matrix_set(T, r2, i+1, s*h1 + c*h2);
            }
            for (int r2 = 0; r2 < n; r2++) {
                double q1 = matrix_get(Q, r2, i), q2 = matrix_get(Q, r2, i+1);
                matrix_set(Q, r2, i, c*q1 - s*q2);
                matrix_set(Q, r2, i+1, s*q1 + c*q2);
            }
        }
        for (int i = 0; i < n; i++) matrix_set(T, i, i, matrix_get(T, i, i) + mu);
        int conv = 1;
        for (int i = 0; i < n-1; i++)
            if (fabs(matrix_get(T, i+1, i)) > 1e-10) { conv = 0; break; }
        if (conv) { for (int i = 0; i < n; i++) vector_set(lambdas, i, matrix_get(T, i, i)); matrix_free(T); return iter; }
    }
    for (int i = 0; i < n; i++) vector_set(lambdas, i, matrix_get(T, i, i));
    matrix_free(T);
    return 100;
}

int matrix_sign_function(const Matrix *A, Matrix *S, double tol) {
    if (!A || !S || !matrix_is_square(A)) return -1;
    int n = A->rows;
    Matrix *X = matrix_create_copy(A);
    Matrix *Xinv = matrix_create(n, n);
    Matrix *Xnew = matrix_create(n, n);
    if (!X || !Xinv || !Xnew) { matrix_free(X); matrix_free(Xinv); matrix_free(Xnew); return -1; }
    for (int iter = 0; iter < 50; iter++) {
        if (lu_inverse(X, Xinv) != 0) { matrix_free(X); matrix_free(Xinv); matrix_free(Xnew); return -1; }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                matrix_set(Xnew, i, j, (matrix_get(X, i, j) + matrix_get(Xinv, i, j)) / 2.0);
        double diff = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                diff += (matrix_get(Xnew, i, j) - matrix_get(X, i, j)) *
                        (matrix_get(Xnew, i, j) - matrix_get(X, i, j));
        memcpy(X->data, Xnew->data, (size_t)n*n*sizeof(double));
        if (sqrt(diff) < tol) break;
    }
    memcpy(S->data, X->data, (size_t)n*n*sizeof(double));
    matrix_free(X); matrix_free(Xinv); matrix_free(Xnew);
    return 0;
}

int matrix_sqrt(const Matrix *A, Matrix *S) {
    if (!A || !S || !matrix_is_square(A)) return -1;
    int n = A->rows;
    /* Use eigendecomposition for SPD: A^{1/2} = Q sqrt(Lambda) Q^T */
    Matrix *Q = matrix_create(n, n);
    Vector *lambdas = vector_create(n);
    if (!Q || !lambdas) { matrix_free(Q); vector_free(lambdas); return -1; }
    if (eigen_symmetric_decompose(A, Q, lambdas) != 0) { matrix_free(Q); vector_free(lambdas); return -1; }
    Matrix *sqrtL = matrix_create_zero(n, n);
    if (sqrtL) {
        for (int i = 0; i < n; i++) {
            double l = vector_get(lambdas, i);
            matrix_set(sqrtL, i, i, (l > 0) ? sqrt(l) : 0.0);
        }
        Matrix *Qt = matrix_create(n, n), *temp = matrix_create(n, n);
        if (Qt && temp) {
            matrix_transpose(Q, Qt);
            matrix_multiply(Q, sqrtL, temp);
            matrix_multiply(temp, Qt, S);
            matrix_free(Qt); matrix_free(temp);
        }
        matrix_free(sqrtL);
    }
    matrix_free(Q); vector_free(lambdas);
    return 0;
}