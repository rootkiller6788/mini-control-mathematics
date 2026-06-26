/**
 * eigen_methods.c — Eigenvalue Computation for Control Stability
 *
 * Implements: L5 Methods (power iteration, QR algorithm), L3 Quantities
 *             (damping ratios, natural frequencies), L8 Advanced (Arnoldi)
 *
 * Eigenvalue analysis is the cornerstone of linear control theory.
 * Stability, transient response, and robustness all depend on the
 * eigenvalue spectrum of the closed-loop system matrix.
 *
 * L7 Applications: Boeing 747 flutter analysis (eigenvalue tracking),
 * Toyota engine knock detection (frequency-domain eigenvalues),
 * nuclear reactor stability (Lyapunov exponents), GPS satellite
 * attitude determination (quaternion eigenvalue methods).
 */

#include "eigen_methods.h"
#include "linear_systems.h"
#include <string.h>
#include <stdio.h>
#include <float.h>

/* ==========================================================================
 * L5: Power Iteration — Largest Eigenvalue by Magnitude
 * ========================================================================== */

int power_iteration(const Matrix *A, double *eigenvalue, Vector *eigenvector,
                    size_t max_iter, double tol) {
    if (!A || A->rows != A->cols || !eigenvalue || !eigenvector
        || eigenvector->n != A->rows) return -1;
    size_t n = A->rows;

    /* Normalize initial guess (use caller's vector) */
    double vnorm = 0.0;
    for (size_t i = 0; i < n; i++) vnorm += eigenvector->data[i] * eigenvector->data[i];
    if (vnorm < MNC_TINY) {
        for (size_t i = 0; i < n; i++) eigenvector->data[i] = 1.0;
        vnorm = (double)n;
    }
    vnorm = sqrt(vnorm);
    for (size_t i = 0; i < n; i++) eigenvector->data[i] /= vnorm;
    double lambda_old = 0.0;

    for (size_t iter = 0; iter < max_iter; iter++) {
        /* w = A * v */
        double *w = (double*)malloc(n * sizeof(double));
        if (!w) return -1;
        for (size_t i = 0; i < n; i++) {
            w[i] = 0.0;
            for (size_t j = 0; j < n; j++)
                w[i] += A->data[i * A->cols + j] * eigenvector->data[j];
        }

        /* Rayleigh quotient: λ = v^T A v / v^T v = v^T w / v^T v */
        double num = 0.0, den = 0.0;
        for (size_t i = 0; i < n; i++) {
            num += eigenvector->data[i] * w[i];
            den += eigenvector->data[i] * eigenvector->data[i];
        }
        double lambda = num / den;

        /* Normalize w */
        double wnorm = 0.0;
        for (size_t i = 0; i < n; i++) wnorm += w[i] * w[i];
        wnorm = sqrt(wnorm);
        if (wnorm < MNC_TINY) { free(w); return -1; }
        for (size_t i = 0; i < n; i++) eigenvector->data[i] = w[i] / wnorm;

        if (fabs(lambda - lambda_old) < tol * fmax(1.0, fabs(lambda))) {
            *eigenvalue = lambda;
            free(w);
            return 0;
        }
        lambda_old = lambda;
        free(w);
    }
    *eigenvalue = lambda_old;
    return 0;
}

/* ==========================================================================
 * L5: Inverse Iteration — Eigenvalue Nearest to Shift
 * ========================================================================== */

int inverse_iteration(const Matrix *A, double shift, double *eigenvalue,
                      Vector *eigenvector, size_t max_iter, double tol) {
    if (!A || A->rows != A->cols || !eigenvalue || !eigenvector
        || eigenvector->n != A->rows) return -1;
    size_t n = A->rows;

    /* Build A - shift*I */
    Matrix *Ashift = matrix_copy(A);
    if (!Ashift) return -1;
    for (size_t i = 0; i < n; i++)
        Ashift->data[i * n + i] -= shift;

    LUResult *lu = lu_decompose(Ashift);
    if (!lu || lu->singular) { lu_free(lu); matrix_free(Ashift); return -1; }

    /* Normalize initial guess */
    double nrm = 0.0;
    for (size_t i = 0; i < n; i++) nrm += eigenvector->data[i] * eigenvector->data[i];
    nrm = sqrt(nrm);
    if (nrm > 0) for (size_t i = 0; i < n; i++) eigenvector->data[i] /= nrm;

    for (size_t iter = 0; iter < max_iter; iter++) {
        Vector *v = vector_create_from(n, eigenvector->data);
        Vector *w = lu_solve(lu, v);
        vector_free(v);
        if (!w) break;

        /* Normalize */
        double wnorm = vector_norm_l2(w);
        if (wnorm < MNC_TINY) { vector_free(w); break; }
        for (size_t i = 0; i < n; i++) eigenvector->data[i] = w->data[i] / wnorm;

        /* Rayleigh quotient for refined eigenvalue */
        double *Av = (double*)malloc(n * sizeof(double));
        for (size_t i = 0; i < n; i++) {
            Av[i] = 0.0;
            for (size_t j = 0; j < n; j++)
                Av[i] += A->data[i * A->cols + j] * eigenvector->data[j];
        }
        double lambda = 0.0;
        for (size_t i = 0; i < n; i++)
            lambda += eigenvector->data[i] * Av[i];
        free(Av);

        if (fabs(lambda - shift) < tol) {
            *eigenvalue = lambda;
            vector_free(w);
            lu_free(lu); matrix_free(Ashift);
            return 0;
        }
        shift = lambda;
        /* Recompute LU with new shift */
        lu_free(lu);
        matrix_free(Ashift);
        Ashift = matrix_copy(A);
        for (size_t i = 0; i < n; i++) Ashift->data[i * n + i] -= shift;
        lu = lu_decompose(Ashift);
        if (!lu || lu->singular) { vector_free(w); lu_free(lu); matrix_free(Ashift); return -1; }

        vector_free(w);
    }
    lu_free(lu); matrix_free(Ashift);
    return 0;
}

/* ==========================================================================
 * L5/L8: Rayleigh Quotient Iteration
 * ========================================================================== */

int rayleigh_quotient_iteration(const Matrix *A, double *eigenvalue,
                                Vector *eigenvector, size_t max_iter, double tol) {
    if (!A || A->rows != A->cols || !eigenvalue || !eigenvector) return -1;
    size_t n = A->rows;

    /* Compute initial Rayleigh quotient */
    double *Av = (double*)malloc(n * sizeof(double));
    if (!Av) return -1;
    for (size_t i = 0; i < n; i++) {
        Av[i] = 0.0;
        for (size_t j = 0; j < n; j++)
            Av[i] += A->data[i * A->cols + j] * eigenvector->data[j];
    }
    double sigma = 0.0;
    for (size_t i = 0; i < n; i++)
        sigma += eigenvector->data[i] * Av[i];
    free(Av);

    for (size_t iter = 0; iter < max_iter; iter++) {
        /* Solve (A - sigma I) y = x_k */
        Matrix *Ashift = matrix_copy(A);
        if (!Ashift) return -1;
        for (size_t i = 0; i < n; i++) Ashift->data[i * n + i] -= sigma;
        LUResult *lu = lu_decompose(Ashift);
        if (!lu || lu->singular) { lu_free(lu); matrix_free(Ashift); return -1; }

        Vector *v = vector_create_from(n, eigenvector->data);
        Vector *y = lu_solve(lu, v);
        vector_free(v); lu_free(lu); matrix_free(Ashift);
        if (!y) return -1;

        /* Normalize and update */
        double ynorm = vector_norm_l2(y);
        if (ynorm < MNC_TINY) { vector_free(y); break; }
        for (size_t i = 0; i < n; i++) eigenvector->data[i] = y->data[i] / ynorm;

        /* New Rayleigh quotient */
        Av = (double*)malloc(n * sizeof(double));
        for (size_t i = 0; i < n; i++) {
            Av[i] = 0.0;
            for (size_t j = 0; j < n; j++)
                Av[i] += A->data[i * A->cols + j] * eigenvector->data[j];
        }
        double sigma_new = 0.0;
        for (size_t i = 0; i < n; i++)
            sigma_new += eigenvector->data[i] * Av[i];
        free(Av);

        if (fabs(sigma_new - sigma) < tol) {
            *eigenvalue = sigma_new;
            vector_free(y);
            return 0;
        }
        sigma = sigma_new;
        vector_free(y);
    }
    *eigenvalue = sigma;
    return 0;
}

/* ==========================================================================
 * L5: QR Algorithm for All Eigenvalues
 * ========================================================================== */

/* Reduce matrix to upper Hessenberg form via Householder reflections */
static void hessenberg_reduction(double *H, size_t n) {
    double *v = (double*)malloc(n * sizeof(double));
    for (size_t k = 0; k < n - 2; k++) {
        /* Zero out column k below subdiagonal */
        double xnorm = 0.0;
        for (size_t i = k + 1; i < n; i++) {
            v[i - (k + 1)] = H[i * n + k];
            xnorm += v[i - (k + 1)] * v[i - (k + 1)];
        }
        xnorm = sqrt(xnorm);
        if (xnorm < MNC_TINY) continue;

        double alpha = (v[0] >= 0) ? -xnorm : xnorm;
        v[0] -= alpha;
        double vnorm2 = v[0] * v[0];
        for (size_t i = 1; i < n - k - 1; i++) vnorm2 += v[i] * v[i];
        double beta = 2.0 / vnorm2;

        /* Apply H to columns k..n-1 */
        for (size_t j = k; j < n; j++) {
            double dot = 0.0;
            for (size_t i = 0; i < n - k - 1; i++)
                dot += v[i] * H[(k + 1 + i) * n + j];
            for (size_t i = 0; i < n - k - 1; i++)
                H[(k + 1 + i) * n + j] -= beta * dot * v[i];
        }
        /* Apply H to rows 0..n-1 */
        for (size_t i = 0; i < n; i++) {
            double dot = 0.0;
            for (size_t j = 0; j < n - k - 1; j++)
                dot += v[j] * H[i * n + (k + 1 + j)];
            for (size_t j = 0; j < n - k - 1; j++)
                H[i * n + (k + 1 + j)] -= beta * dot * v[j];
        }
    }
    free(v);
}

/* Implicit double-shift QR step (Francis) on upper Hessenberg matrix */
static void qr_double_shift_step(double *H, size_t n, size_t *iter_count) {
    double s = H[(n - 1) * n + (n - 1)] + H[(n - 2) * n + (n - 2)];
    double t = H[(n - 1) * n + (n - 1)] * H[(n - 2) * n + (n - 2)]
             - H[(n - 1) * n + (n - 2)] * H[(n - 2) * n + (n - 1)];

    /* First column of (H^2 - s*H + t*I) */
    double x = H[0] * H[0] + H[0 * n + 1] * H[1 * n + 0] - s * H[0] + t;
    double y = H[1 * n + 0] * (H[0] + H[1 * n + 1] - s);
    double z = H[1 * n + 0] * H[2 * n + 1];

    for (size_t k = 0; k < n - 1; k++) {
        /* Bulge-chasing: apply Householder to create zero at (k+2, k) */
        double *v = (double*)malloc(3 * sizeof(double));
        v[0] = x; v[1] = y; v[2] = (k < n - 2) ? z : 0.0;
        double vnorm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        if (vnorm < MNC_TINY) { free(v); continue; }
        v[0] /= vnorm; v[1] /= vnorm; v[2] /= vnorm;

        /* Apply from left: H(k:k+2, k:n-1) */
        for (size_t j = k; j < n; j++) {
            double dot = v[0] * H[k * n + j] + v[1] * H[(k+1) * n + j];
            if (k < n - 2) dot += v[2] * H[(k+2) * n + j];
            double tau = 2.0 * dot;
            H[k * n + j] -= tau * v[0];
            H[(k+1) * n + j] -= tau * v[1];
            if (k < n - 2) H[(k+2) * n + j] -= tau * v[2];
        }
        /* Apply from right: H(0:min(k+3,n-1), k:k+2) */
        size_t maxrow = (k + 3 < n) ? k + 3 : n;
        for (size_t i = 0; i < maxrow; i++) {
            double dot = v[0] * H[i * n + k] + v[1] * H[i * n + (k+1)];
            if (k < n - 2) dot += v[2] * H[i * n + (k+2)];
            double tau = 2.0 * dot;
            H[i * n + k] -= tau * v[0];
            H[i * n + (k+1)] -= tau * v[1];
            if (k < n - 2) H[i * n + (k+2)] -= tau * v[2];
        }
        free(v);

        /* Next bulge */
        if (k < n - 2) {
            x = H[(k+1) * n + k];
            y = H[(k+2) * n + k];
            z = (k < n - 3) ? H[(k+3) * n + k] : 0.0;
        }
    }
    (*iter_count)++;
}

/* Extract eigenvalues from quasi-triangular Hessenberg */
static void extract_eigenvalues_real(const double *H, size_t n,
                                      double *real_parts, double *imag_parts) {
    size_t i = 0;
    while (i < n) {
        if (i < n - 1 && fabs(H[(i + 1) * n + i]) > MNC_SQRT_EPSILON) {
            /* 2×2 diagonal block → complex conjugate pair */
            double a = H[i * n + i], b = H[i * n + (i+1)];
            double c = H[(i+1) * n + i], d = H[(i+1) * n + (i+1)];
            double tr = a + d, det = a * d - b * c;
            double disc = tr * tr - 4.0 * det;
            if (disc < 0) {
                real_parts[i] = real_parts[i+1] = tr / 2.0;
                imag_parts[i] = sqrt(-disc) / 2.0;
                imag_parts[i+1] = -imag_parts[i];
            } else {
                real_parts[i] = (tr + sqrt(disc)) / 2.0; imag_parts[i] = 0.0;
                real_parts[i+1] = (tr - sqrt(disc)) / 2.0; imag_parts[i+1] = 0.0;
            }
            i += 2;
        } else {
            real_parts[i] = H[i * n + i];
            imag_parts[i] = 0.0;
            i++;
        }
    }
}

EigenResult* qr_algorithm_eigenvalues(const Matrix *A) {
    if (!A || A->rows != A->cols) return NULL;
    size_t n = A->rows;

    EigenResult *eig = (EigenResult*)malloc(sizeof(EigenResult));
    if (!eig) return NULL;
    eig->n = n;
    eig->real_parts = (double*)malloc(n * sizeof(double));
    eig->imag_parts = (double*)malloc(n * sizeof(double));
    eig->eigenvectors = NULL;
    eig->converged = 1;
    if (!eig->real_parts || !eig->imag_parts) { eigen_result_free(eig); return NULL; }

    double *H = (double*)malloc(n * n * sizeof(double));
    if (!H) { eigen_result_free(eig); return NULL; }
    memcpy(H, A->data, n * n * sizeof(double));

    /* Step 1: Reduce to Hessenberg form */
    hessenberg_reduction(H, n);

    /* Step 2: QR iteration with double shift */
    size_t iter_total = 0;
    size_t p = n;
    while (p > 1 && iter_total < 100 * n) {
        /* Check for deflation: if subdiagonal element is small, deflate */
        while (p > 1 && fabs(H[(p-1) * n + (p-2)]) < MNC_SQRT_EPSILON) p--;
        if (p == 1) break;

        /* If there's a 2×2 block at bottom, deflate it */
        if (p > 2 && fabs(H[(p-2) * n + (p-3)]) < MNC_SQRT_EPSILON) p--;

        /* QR double shift step on leading p×p block */
        qr_double_shift_step(H, p, &iter_total);
        if (iter_total > 100 * n) { eig->converged = 0; break; }
    }

    extract_eigenvalues_real(H, n, eig->real_parts, eig->imag_parts);
    free(H);
    return eig;
}

/* ==========================================================================
 * L5: Symmetric Eigenvalue via Tridiagonalization + QR
 * ========================================================================== */

EigenResult* eigen_symmetric(const Matrix *A) {
    if (!A || A->rows != A->cols) return NULL;
    return qr_algorithm_eigenvalues(A); /* For now, use general QR which works for symmetric too */
}

/* ==========================================================================
 * L3/L5: Spectral Analysis for Control
 * ========================================================================== */

SpectralRadius* spectral_radius_compute(const EigenResult *eig) {
    if (!eig) return NULL;
    SpectralRadius *sr = (SpectralRadius*)malloc(sizeof(SpectralRadius));
    if (!sr) return NULL;

    sr->spectral_radius = 0.0;
    sr->spectral_abscissa = -1e308;
    sr->min_real_part = 1e308;
    sr->max_imag_part = 0.0;
    sr->is_hurwitz = 1;
    sr->is_schur = 1;

    for (size_t i = 0; i < eig->n; i++) {
        double mag = sqrt(eig->real_parts[i] * eig->real_parts[i]
                        + eig->imag_parts[i] * eig->imag_parts[i]);
        if (mag > sr->spectral_radius) sr->spectral_radius = mag;
        if (eig->real_parts[i] > sr->spectral_abscissa)
            sr->spectral_abscissa = eig->real_parts[i];
        if (eig->real_parts[i] < sr->min_real_part)
            sr->min_real_part = eig->real_parts[i];
        if (fabs(eig->imag_parts[i]) > sr->max_imag_part)
            sr->max_imag_part = fabs(eig->imag_parts[i]);

        if (eig->real_parts[i] >= 0) sr->is_hurwitz = 0;
        if (mag >= 1.0) sr->is_schur = 0;
    }
    return sr;
}

void spectral_radius_free(SpectralRadius *sr) { free(sr); }

void eigen_result_free(EigenResult *eig) {
    if (!eig) return;
    free(eig->real_parts);
    free(eig->imag_parts);
    matrix_free(eig->eigenvectors);
    free(eig);
}

int eigen_is_stable_continuous(const EigenResult *eig, double tol) {
    if (!eig) return 0;
    for (size_t i = 0; i < eig->n; i++)
        if (eig->real_parts[i] >= -tol) return 0;
    return 1;
}

int eigen_is_stable_discrete(const EigenResult *eig, double tol) {
    if (!eig) return 0;
    for (size_t i = 0; i < eig->n; i++)
        if (eig->real_parts[i] * eig->real_parts[i]
            + eig->imag_parts[i] * eig->imag_parts[i] >= 1.0 - tol) return 0;
    return 1;
}

Vector* eigen_damping_ratios(const EigenResult *eig) {
    if (!eig) return NULL;
    Vector *zeta = vector_create(eig->n);
    if (!zeta) return NULL;
    for (size_t i = 0; i < eig->n; i++) {
        double mag = sqrt(eig->real_parts[i] * eig->real_parts[i]
                        + eig->imag_parts[i] * eig->imag_parts[i]);
        zeta->data[i] = (mag > MNC_TINY) ? -eig->real_parts[i] / mag : 0.0;
    }
    return zeta;
}

Vector* eigen_natural_frequencies(const EigenResult *eig) {
    if (!eig) return NULL;
    Vector *wn = vector_create(eig->n);
    if (!wn) return NULL;
    for (size_t i = 0; i < eig->n; i++)
        wn->data[i] = sqrt(eig->real_parts[i] * eig->real_parts[i]
                         + eig->imag_parts[i] * eig->imag_parts[i]);
    return wn;
}

/* ==========================================================================
 * L8: Arnoldi Iteration for Large-Scale Problems
 * ========================================================================== */

int arnoldi_iteration(const Matrix *A, const Vector *v1, size_t m,
                      Matrix **V_out, Matrix **H_out) {
    if (!A || !v1 || A->rows != A->cols || v1->n != A->rows || m < 1) return -1;
    size_t n = A->rows;
    if (m > n) m = n;

    Matrix *V = matrix_create_zeros(n, m + 1);
    Matrix *H = matrix_create_zeros(m + 1, m);
    if (!V || !H) { matrix_free(V); matrix_free(H); return -1; }

    /* Normalize v1 */
    double beta = vector_norm_l2(v1);
    if (beta < MNC_TINY) { matrix_free(V); matrix_free(H); return -1; }
    for (size_t i = 0; i < n; i++) V->data[i * V->cols + 0] = v1->data[i] / beta;

    for (size_t j = 0; j < m; j++) {
        /* w = A * v_j */
        double *w = (double*)calloc(n, sizeof(double));
        if (!w) { matrix_free(V); matrix_free(H); return -1; }
        for (size_t i = 0; i < n; i++)
            for (size_t k = 0; k < n; k++)
                w[i] += A->data[i * A->cols + k] * V->data[k * V->cols + j];

        /* Modified Gram-Schmidt */
        for (size_t i = 0; i <= j; i++) {
            double dot = 0.0;
            for (size_t k = 0; k < n; k++)
                dot += V->data[k * V->cols + i] * w[k];
            H->data[i * H->cols + j] = dot;
            for (size_t k = 0; k < n; k++)
                w[k] -= dot * V->data[k * V->cols + i];
        }

        double h = 0.0;
        for (size_t k = 0; k < n; k++) h += w[k] * w[k];
        h = sqrt(h);
        H->data[(j + 1) * H->cols + j] = h;

        if (h > MNC_TINY && j + 1 < m + 1) {
            for (size_t k = 0; k < n; k++)
                V->data[k * V->cols + (j + 1)] = w[k] / h;
        }
        free(w);
        if (h < MNC_TINY) break;  /* Invariant subspace found */
    }

    *V_out = V;
    *H_out = H;
    return 0;
}

EigenResult* arnoldi_eigenvalues(const Matrix *A, const Vector *v1,
                                 size_t m, size_t n_eigen) {
    if (!A || !v1) return NULL;
    Matrix *V = NULL, *H = NULL;
    if (arnoldi_iteration(A, v1, m, &V, &H) != 0) return NULL;

    /* Extract leading m×m block of H and compute eigenvalues */
    Matrix *Hm = matrix_create(m, m);
    if (!Hm) { matrix_free(V); matrix_free(H); return NULL; }
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < m; j++)
            Hm->data[i * m + j] = H->data[i * H->cols + j];

    EigenResult *eig = qr_algorithm_eigenvalues(Hm);
    (void)n_eigen;
    matrix_free(Hm);
    matrix_free(V);
    matrix_free(H);
    return eig;
}
