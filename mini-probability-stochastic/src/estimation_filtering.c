/**
 * @file estimation_filtering.c
 * @brief Kalman filter (linear, extended, unscented), particle filter (SIR),
 *        Kalman smoother, log-likelihood computation, AIC/BIC,
 *        and matrix utility functions.
 *
 * Knowledge: L5 engineering methods (Kalman, EKF, UKF, particle filter),
 * L6 engineering problems (state estimation, sensor fusion),
 * L7 applications (GPS/IMU, radar tracking).
 */
#include "estimation_filtering.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double sq(double x) { return x * x; }

static double lcg_uniform(uint64_t *seed) {
    *seed = *seed * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)((*seed) >> 32) / 4294967296.0;
}

static double box_muller(uint64_t *seed) {
    double u1 = lcg_uniform(seed), u2 = lcg_uniform(seed);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* ================================================================
 * Matrix Utilities
 * ================================================================ */

/*
 * mat_mul_n — C = A * B for n-by-n square matrices. O(n^3).
 */
void mat_mul_n(const prob_real_t *A, const prob_real_t *B,
               prob_real_t *C, size_t n) {
    if (!A || !B || !C || n == 0) return;
    size_t i, j, k;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            double s = 0.0;
            for (k = 0; k < n; k++)
                s += A[i * n + k] * B[k * n + j];
            C[i * n + j] = s;
        }
}

/*
 * mat_vec_mul — y = A * x, A is rows-by-cols. O(rows*cols).
 */
void mat_vec_mul(const prob_real_t *A, const prob_real_t *x,
                 prob_real_t *y, size_t rows, size_t cols) {
    if (!A || !x || !y || rows == 0 || cols == 0) return;
    size_t i, j;
    for (i = 0; i < rows; i++) {
        double s = 0.0;
        for (j = 0; j < cols; j++)
            s += A[i * cols + j] * x[j];
        y[i] = s;
    }
}

/*
 * mat_transpose — B = A^T. O(rows*cols).
 */
void mat_transpose(const prob_real_t *A, prob_real_t *B,
                   size_t rows, size_t cols) {
    if (!A || !B || rows == 0 || cols == 0) return;
    size_t i, j;
    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            B[j * rows + i] = A[i * cols + j];
}

/*
 * mat_cholesky_solve — solve Ax = b via Cholesky decomposition.
 * A must be n-by-n symmetric positive definite.
 * O(n^3). Returns 0 on success, -1 if non-positive-definite.
 */
int mat_cholesky_solve(const prob_real_t *A, const prob_real_t *b,
                       prob_real_t *x, size_t n) {
    if (!A || !b || !x || n == 0) return -1;
    double *L = (double *)calloc(n * n, sizeof(double));
    double *y = (double *)calloc(n, sizeof(double));
    if (!L || !y) { free(L); free(y); return -1; }
    size_t i, j, k;
    /* Cholesky: A = L * L^T */
    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            double s = A[i * n + j];
            for (k = 0; k < j; k++)
                s -= L[i * n + k] * L[j * n + k];
            if (i == j) {
                if (s <= 0.0) { free(L); free(y); return -1; }
                L[i * n + i] = sqrt(s);
            } else {
                L[i * n + j] = s / L[j * n + j];
            }
        }
    }
    /* Forward substitution: L*y = b */
    for (i = 0; i < n; i++) {
        double s = b[i];
        for (j = 0; j < i; j++)
            s -= L[i * n + j] * y[j];
        y[i] = s / L[i * n + i];
    }
    /* Back substitution: L^T * x = y */
    for (i = n; i > 0; ) { i--;
        double s = y[i];
        for (j = i + 1; j < n; j++)
            s -= L[j * n + i] * x[j];
        x[i] = s / L[i * n + i];
    }
    free(L); free(y);
    return 0;
}

/*
 * mat_gaussian_solve — solve via Gaussian elimination with partial pivoting.
 * A and b are modified in-place. O(n^3).
 */
int mat_gaussian_solve(prob_real_t *A, prob_real_t *b, size_t n,
                       prob_real_t *x) {
    if (!A || !b || !x || n == 0) return -1;
    size_t i, j, k;
    for (k = 0; k < n; k++) {
        size_t max_r = k;
        double max_v = fabs(A[k * n + k]);
        for (i = k + 1; i < n; i++)
            if (fabs(A[i * n + k]) > max_v) {
                max_v = fabs(A[i * n + k]); max_r = i;
            }
        if (max_v < 1e-15) return -1;
        if (max_r != k) {
            for (j = k; j < n; j++) {
                double t = A[k*n+j]; A[k*n+j] = A[max_r*n+j];
                A[max_r*n+j] = t;
            }
            double tb = b[k]; b[k] = b[max_r]; b[max_r] = tb;
        }
        for (i = k + 1; i < n; i++) {
            double f = A[i * n + k] / A[k * n + k];
            for (j = k; j < n; j++) A[i*n+j] -= f * A[k*n+j];
            b[i] -= f * b[k];
        }
    }
    for (i = n; i > 0; ) { i--;
        double s = b[i];
        for (j = i + 1; j < n; j++) s -= A[i*n+j] * x[j];
        x[i] = s / A[i*n+i];
    }
    return 0;
}

/*
 * mat_inverse — B = A^{-1} via Gauss-Jordan on augmented matrix.
 * O(n^3). Returns 0 on success, -1 if singular.
 */
int mat_inverse(const prob_real_t *A, prob_real_t *B, size_t n) {
    if (!A || !B || n == 0) return -1;
    /* Build augmented matrix [A | I] */
    double *aug = (double *)malloc(n * (2 * n) * sizeof(double));
    if (!aug) return -1;
    size_t i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            aug[i * (2*n) + j] = A[i * n + j];
            aug[i * (2*n) + n + j] = (i == j) ? 1.0 : 0.0;
        }
    }
    /* Forward elimination */
    size_t k;
    for (k = 0; k < n; k++) {
        size_t max_r = k;
        double max_v = fabs(aug[k * (2*n) + k]);
        for (i = k + 1; i < n; i++)
            if (fabs(aug[i * (2*n) + k]) > max_v) {
                max_v = fabs(aug[i * (2*n) + k]); max_r = i;
            }
        if (max_v < 1e-15) { free(aug); return -1; }
        if (max_r != k)
            for (j = 0; j < 2*n; j++) {
                double t = aug[k*(2*n)+j];
                aug[k*(2*n)+j] = aug[max_r*(2*n)+j];
                aug[max_r*(2*n)+j] = t;
            }
        double piv = aug[k * (2*n) + k];
        for (j = 0; j < 2*n; j++) aug[k*(2*n)+j] /= piv;
        for (i = 0; i < n; i++) {
            if (i == k) continue;
            double f = aug[i * (2*n) + k];
            for (j = 0; j < 2*n; j++)
                aug[i*(2*n)+j] -= f * aug[k*(2*n)+j];
        }
    }
    /* Extract inverse */
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            B[i * n + j] = aug[i * (2*n) + n + j];
    free(aug);
    return 0;
}

/* ================================================================
 * L5: Linear Kalman Filter (Kalman, J. Basic Eng., 1960)
 * ================================================================ */

/*
 * kalman_predict — time update (prior).
 * x_pred = F * x + B * u
 * P_pred = F * P * F^T + Q
 * O(d^3 + d^2*m) where d=state_dim. Independent knowledge point:
 * the Kalman predict step propagates the state and its uncertainty.
 */
void kalman_predict(const kalman_model_t *model,
                    const prob_real_t *x, const prob_real_t *P,
                    const prob_real_t *u,
                    prob_real_t *x_pred, prob_real_t *P_pred) {
    if (!model || !x || !P || !x_pred || !P_pred) return;
    size_t d = model->state_dim;
    /* x_pred = F * x */
    mat_vec_mul(model->F, x, x_pred, d, d);
    /* Add B * u if control exists */
    if (model->B && u && model->input_dim > 0) {
        double *Bu = (double *)calloc(d, sizeof(double));
        if (Bu) {
            mat_vec_mul(model->B, u, Bu, d, model->input_dim);
            size_t i;
            for (i = 0; i < d; i++) x_pred[i] += Bu[i];
            free(Bu);
        }
    }
    /* P_pred = F * P * F^T + Q */
    double *FP = (double *)calloc(d * d, sizeof(double));
    double *FT = (double *)calloc(d * d, sizeof(double));
    if (!FP || !FT) { free(FP); free(FT); return; }
    mat_mul_n(model->F, P, FP, d);          /* FP = F * P */
    mat_transpose(model->F, FT, d, d);      /* FT = F^T */
    mat_mul_n(FP, FT, P_pred, d);           /* FP * F^T */
    size_t i;
    for (i = 0; i < d * d; i++) P_pred[i] += model->Q[i];
    free(FP); free(FT);
}

/*
 * kalman_update — measurement update (posterior).
 * K = P_pred * H^T * inv(H * P_pred * H^T + R)  (Kalman gain)
 * x = x_pred + K * (y - H * x_pred)             (state correction)
 * P = (I - K * H) * P_pred                       (covariance correction)
 * O(d^3 + d^2*m + d*m^2 + m^3). Independent knowledge point:
 * optimal fusion of prediction and observation under Gaussian noise.
 */
void kalman_update(const kalman_model_t *model,
                   const prob_real_t *x_pred, const prob_real_t *P_pred,
                   const prob_real_t *y,
                   prob_real_t *x, prob_real_t *P,
                   prob_real_t *K) {
    if (!model || !x_pred || !P_pred || !y || !x || !P) return;
    size_t d = model->state_dim, m = model->obs_dim;
    /* Innovation: z = y - H * x_pred */
    double *z = (double *)malloc(m * sizeof(double));
    if (!z) return;
    memcpy(z, y, m * sizeof(double));
    {
        double *Hx = (double *)calloc(m, sizeof(double));
        if (Hx) {
            mat_vec_mul(model->H, x_pred, Hx, m, d);
            size_t i;
            for (i = 0; i < m; i++) z[i] -= Hx[i];
            free(Hx);
        }
    }
    /* S = H * P_pred * H^T + R */
    double *PHt = (double *)calloc(d * m, sizeof(double));
    double *S   = (double *)calloc(m * m, sizeof(double));
    double *Ht  = (double *)calloc(d * m, sizeof(double));
    double *K_local = (double *)calloc(d * m, sizeof(double));
    double *Sinv = (double *)calloc(m * m, sizeof(double));
    if (!PHt || !S || !Ht || !K_local || !Sinv) {
        free(z); free(PHt); free(S); free(Ht); free(K_local); free(Sinv);
        return;
    }
    mat_transpose(model->H, Ht, m, d);            /* Ht = H^T */
    /* PHt = P_pred * H^T */
    size_t i, j, kk;
    for (i = 0; i < d; i++)
        for (j = 0; j < m; j++) {
            double s = 0.0;
            for (kk = 0; kk < d; kk++)
                s += P_pred[i * d + kk] * Ht[kk * m + j];
            PHt[i * m + j] = s;
        }
    /* S = H * PHt + R */
    for (i = 0; i < m; i++)
        for (j = 0; j < m; j++) {
            double s = 0.0;
            for (kk = 0; kk < d; kk++)
                s += model->H[i * d + kk] * PHt[kk * m + j];
            S[i * m + j] = s + model->R[i * m + j];
        }
    /* K = PHt * S^{-1} — solve S * K^T = PHt^T */
    {
        double *S_copy = (double *)malloc(m * m * sizeof(double));
        double *rhs = (double *)malloc(m * d * sizeof(double));
        if (!S_copy || !rhs) { free(S_copy); free(rhs);
            free(z); free(PHt); free(S); free(Ht); free(K_local); free(Sinv);
            return; }
        memcpy(S_copy, S, m * m * sizeof(double));
        /* For each column of K, solve S * col = PHt^T_col */
        size_t col;
        for (col = 0; col < d; col++) {
            for (i = 0; i < m; i++) rhs[i] = PHt[col * m + i];
            /* Gaussian elimination on S_copy (copy each time) */
            double *Sc = (double *)malloc(m * m * sizeof(double));
            double *rc = (double *)malloc(m * sizeof(double));
            if (!Sc || !rc) { free(Sc); free(rc); continue; }
            memcpy(Sc, S, m * m * sizeof(double));
            memcpy(rc, rhs, m * sizeof(double));
            /* Solve Sc * kgain_col = rc */
            size_t kp;
            for (kp = 0; kp < m; kp++) {
                double piv = Sc[kp * m + kp];
                if (fabs(piv) < 1e-15) break;
                for (j = kp; j < m; j++) Sc[kp*m+j] /= piv;
                rc[kp] /= piv;
                for (i = 0; i < m; i++) {
                    if (i == kp) continue;
                    double f = Sc[i*m+kp];
                    for (j = kp; j < m; j++) Sc[i*m+j] -= f * Sc[kp*m+j];
                    rc[i] -= f * rc[kp];
                }
            }
            for (i = 0; i < m; i++) K_local[i * d + col] = rc[i];
            free(Sc); free(rc);
        }
        free(S_copy); free(rhs);
    }
    if (K) memcpy(K, K_local, d * m * sizeof(double));
    /* x = x_pred + K * z */
    for (i = 0; i < d; i++) {
        double kz = 0.0;
        for (j = 0; j < m; j++) kz += K_local[j * d + i] * z[j];
        x[i] = x_pred[i] + kz;
    }
    /* P = (I - K*H) * P_pred */
    {
        double *KH = (double *)calloc(d * d, sizeof(double));
        if (KH) {
            for (i = 0; i < d; i++)
                for (j = 0; j < d; j++) {
                    double s = 0.0;
                    size_t k2;
                    for (k2 = 0; k2 < m; k2++)
                        s += K_local[k2 * d + i] * model->H[k2 * d + j];
                    KH[i * d + j] = s;
                }
            for (i = 0; i < d; i++)
                for (j = 0; j < d; j++) {
                    double s = 0.0;
                    size_t k2;
                    for (k2 = 0; k2 < d; k2++)
                        s += ((i==k2?1.0:0.0) - KH[i*d+k2]) * P_pred[k2*d+j];
                    P[i*d+j] = s;
                }
            free(KH);
        }
    }
    free(z); free(PHt); free(S); free(Ht); free(K_local); free(Sinv);
}

/*
 * kalman_filter_run — full forward pass over T observations.
 * O(T * d^3). Independent knowledge point: the complete Kalman
 * filtering pipeline for time-series state estimation.
 */
void kalman_filter_run(const kalman_model_t *model,
                       const prob_real_t *observations,
                       const prob_real_t *controls,
                       size_t T,
                       const prob_real_t *x_init,
                       const prob_real_t *P_init,
                       prob_real_t *x_estimates,
                       prob_real_t *P_estimates) {
    if (!model || !observations || !x_init || !P_init
        || !x_estimates || !P_estimates || T == 0) return;
    size_t d = model->state_dim;
    double *x = (double *)malloc(d * sizeof(double));
    double *P = (double *)malloc(d * d * sizeof(double));
    double *x_pred = (double *)malloc(d * sizeof(double));
    double *P_pred = (double *)malloc(d * d * sizeof(double));
    if (!x || !P || !x_pred || !P_pred) {
        free(x); free(P); free(x_pred); free(P_pred); return;
    }
    memcpy(x, x_init, d * sizeof(double));
    memcpy(P, P_init, d * d * sizeof(double));
    size_t t;
    for (t = 0; t < T; t++) {
        const double *u = (controls) ? &controls[t * model->input_dim] : NULL;
        kalman_predict(model, x, P, u, x_pred, P_pred);
        kalman_update(model, x_pred, P_pred,
                      &observations[t * model->obs_dim],
                      x, P, NULL);
        memcpy(&x_estimates[t * d], x, d * sizeof(double));
        memcpy(&P_estimates[t * d * d], P, d * d * sizeof(double));
    }
    free(x); free(P); free(x_pred); free(P_pred);
}

/*
 * kalman_smoother — Rauch-Tung-Striebel backward pass.
 * x_{k|T} = x_{k|k} + G_k * (x_{k+1|T} - x_{k+1|k})
 * P_{k|T} = P_{k|k} + G_k * (P_{k+1|T} - P_{k+1|k}) * G_k^T
 * where G_k = P_{k|k} * F^T * P_{k+1|k}^{-1}.
 * O(T * d^3). Independent knowledge point: optimal retrospective
 * state estimation using all available data.
 */
void kalman_smoother(const kalman_model_t *model,
                     const prob_real_t *x_filtered,
                     const prob_real_t *P_filtered,
                     const prob_real_t *x_predicted,
                     const prob_real_t *P_predicted,
                     size_t T,
                     prob_real_t *x_smoothed,
                     prob_real_t *P_smoothed) {
    if (!model || !x_filtered || !P_filtered
        || !x_predicted || !P_predicted
        || !x_smoothed || !P_smoothed || T == 0) return;
    size_t d = model->state_dim;
    /* Copy last estimate */
    memcpy(&x_smoothed[(T-1)*d], &x_filtered[(T-1)*d], d * sizeof(double));
    memcpy(&P_smoothed[(T-1)*d*d], &P_filtered[(T-1)*d*d],
           d * d * sizeof(double));
    size_t k;
    for (k = T - 1; k > 0; k--) {
        size_t km1 = k - 1;
        /* G = P_{k|k} * F^T * inv(P_{k+1|k}) */
        double *FT = (double *)malloc(d * d * sizeof(double));
        double *PFt = (double *)malloc(d * d * sizeof(double));
        double *Pkp1k_inv = (double *)malloc(d * d * sizeof(double));
        double *G = (double *)malloc(d * d * sizeof(double));
        if (!FT || !PFt || !Pkp1k_inv || !G) {
            free(FT); free(PFt); free(Pkp1k_inv); free(G); return;
        }
        mat_transpose(model->F, FT, d, d);
        mat_mul_n(&P_filtered[km1*d*d], FT, PFt, d);
        /* Invert P_predicted[k] */
        mat_inverse(&P_predicted[k*d*d], Pkp1k_inv, d);
        mat_mul_n(PFt, Pkp1k_inv, G, d);
        /* x_smoothed[km1] = x_filtered[km1] + G*(x_smoothed[k]-x_predicted[k]) */
        size_t i, j, kk;
        double *diff = (double *)malloc(d * sizeof(double));
        if (diff) {
            for (i = 0; i < d; i++)
                diff[i] = x_smoothed[k*d + i] - x_predicted[k*d + i];
            for (i = 0; i < d; i++) {
                double gd = 0.0;
                for (j = 0; j < d; j++) gd += G[i*d+j] * diff[j];
                x_smoothed[km1*d + i] = x_filtered[km1*d + i] + gd;
            }
            free(diff);
        }
        /* P_smoothed[km1] = P_filtered[km1] + G*(P_smoothed[k]-P_predicted[k])*G^T */
        double *Pd = (double *)malloc(d * d * sizeof(double));
        double *GT = (double *)malloc(d * d * sizeof(double));
        double *tmp1 = (double *)malloc(d * d * sizeof(double));
        double *tmp2 = (double *)malloc(d * d * sizeof(double));
        if (Pd && GT && tmp1 && tmp2) {
            for (i = 0; i < d*d; i++)
                Pd[i] = P_smoothed[k*d*d+i] - P_predicted[k*d*d+i];
            mat_mul_n(G, Pd, tmp1, d);
            mat_transpose(G, GT, d, d);
            mat_mul_n(tmp1, GT, tmp2, d);
            for (i = 0; i < d*d; i++)
                P_smoothed[km1*d*d+i] = P_filtered[km1*d*d+i] + tmp2[i];
        }
        free(Pd); free(GT); free(tmp1); free(tmp2);
        free(FT); free(PFt); free(Pkp1k_inv); free(G);
    }
}

/* ================================================================
 * L5: Extended Kalman Filter
 * ================================================================ */

/*
 * ekf_predict — EKF predict with finite-difference Jacobians.
 * F_ij = df_i/dx_j computed by central differences.
 * O(d^3 + d^2 + d*input_dim). Independent knowledge point:
 * first-order linearisation for mildly nonlinear systems.
 */
void ekf_predict(const nlinear_model_t *model,
                 const prob_real_t *x, const prob_real_t *P,
                 const prob_real_t *u,
                 prob_real_t *x_pred, prob_real_t *P_pred,
                 prob_real_t eps) {
    if (!model || !x || !P || !x_pred || !P_pred || eps <= 0.0) return;
    size_t d = model->state_dim;
    /* Compute f(x, u) directly */
    model->f(x, u, x_pred);
    /* Jacobian via central finite differences */
    double *F = (double *)calloc(d * d, sizeof(double));
    double *xp = (double *)malloc(d * sizeof(double));
    double *fp = (double *)malloc(d * sizeof(double));
    double *fm = (double *)malloc(d * sizeof(double));
    if (!F || !xp || !fp || !fm) {
        free(F); free(xp); free(fp); free(fm); return;
    }
    size_t i, j;
    for (j = 0; j < d; j++) {
        memcpy(xp, x, d * sizeof(double)); xp[j] += eps;
        memcpy(fm, x, d * sizeof(double));
        /* Need also fm = f(x - eps*e_j) — actually use central diff */
        double *x_tmp = (double *)malloc(d * sizeof(double));
        if (!x_tmp) continue;
        memcpy(x_tmp, x, d * sizeof(double)); x_tmp[j] -= eps;
        /* Compute f(x+eps) and f(x-eps) */
        double *f_plus = (double *)malloc(d * sizeof(double));
        double *f_minus = (double *)malloc(d * sizeof(double));
        if (f_plus && f_minus) {
            model->f(xp, u, f_plus);
            model->f(x_tmp, u, f_minus);
            for (i = 0; i < d; i++)
                F[i * d + j] = (f_plus[i] - f_minus[i]) / (2.0 * eps);
        }
        free(f_plus); free(f_minus); free(x_tmp);
    }
    /* P_pred = F * P * F^T + Q */
    double *FP = (double *)malloc(d * d * sizeof(double));
    double *FT = (double *)malloc(d * d * sizeof(double));
    if (FP && FT) {
        mat_mul_n(F, P, FP, d);
        mat_transpose(F, FT, d, d);
        mat_mul_n(FP, FT, P_pred, d);
        for (i = 0; i < d * d; i++) P_pred[i] += model->Q[i];
    }
    free(F); free(xp); free(fp); free(fm); free(FP); free(FT);
}

/*
 * ekf_update — EKF measurement update with Jacobian of h.
 */
void ekf_update(const nlinear_model_t *model,
                const prob_real_t *x_pred, const prob_real_t *P_pred,
                const prob_real_t *y,
                prob_real_t *x, prob_real_t *P,
                prob_real_t eps) {
    if (!model || !x_pred || !P_pred || !y || !x || !P || eps <= 0.0)
        return;
    size_t d = model->state_dim, m = model->obs_dim;
    /* Jacobian H = dh/dx via central finite differences */
    double *H = (double *)calloc(m * d, sizeof(double));
    if (!H) return;
    size_t i, j;
    for (j = 0; j < d; j++) {
        double *xp = (double *)malloc(d * sizeof(double));
        double *xm = (double *)malloc(d * sizeof(double));
        double *hp = (double *)malloc(m * sizeof(double));
        double *hm = (double *)malloc(m * sizeof(double));
        if (!xp || !xm || !hp || !hm) {
            free(xp); free(xm); free(hp); free(hm); continue;
        }
        memcpy(xp, x_pred, d * sizeof(double)); xp[j] += eps;
        memcpy(xm, x_pred, d * sizeof(double)); xm[j] -= eps;
        model->h(xp, hp); model->h(xm, hm);
        for (i = 0; i < m; i++)
            H[i * d + j] = (hp[i] - hm[i]) / (2.0 * eps);
        free(xp); free(xm); free(hp); free(hm);
    }
    /* Build a kalman_model_t from the linearised H for the update */
    kalman_model_t km;
    km.state_dim = d; km.input_dim = 0; km.obs_dim = m;
    km.F = NULL; km.B = NULL; km.H = H;
    km.Q = model->Q; km.R = model->R;
    kalman_update(&km, x_pred, P_pred, y, x, P, NULL);
    free(H);
}

/* ================================================================
 * L5: Unscented Kalman Filter (Julier & Uhlmann, Proc. IEEE, 2004)
 * ================================================================ */

/*
 * ukf_step — single UKF predict-update.
 * 2n+1 sigma points capture mean and covariance to 2nd order exactly
 * (3rd for Gaussian priors). No Jacobians required.
 * O(d^3). Independent knowledge point: derivative-free nonlinear
 * filtering with sigma-point deterministic sampling.
 */
void ukf_step(const nlinear_model_t *model,
              const prob_real_t *x, const prob_real_t *P,
              const prob_real_t *u, const prob_real_t *y,
              prob_real_t *x_out, prob_real_t *P_out,
              prob_real_t alpha, prob_real_t beta, prob_real_t kappa) {
    if (!model || !x || !P || !u || !y || !x_out || !P_out) return;
    size_t d = model->state_dim, m = model->obs_dim;
    size_t n_sigma = 2 * d + 1;
    prob_real_t lambda = alpha * alpha * (d + kappa) - d;
    /* Weights */
    double Wm0 = lambda / (d + lambda);
    double Wc0 = Wm0 + (1.0 - alpha * alpha + beta);
    double Wi  = 0.5 / (d + lambda);
    /* Generate sigma points */
    double *sigma_x = (double *)malloc(n_sigma * d * sizeof(double));
    double *sigma_y = (double *)malloc(n_sigma * m * sizeof(double));
    if (!sigma_x || !sigma_y) { free(sigma_x); free(sigma_y); return; }
    /* Cholesky of (d+lambda)*P */
    double *L = (double *)calloc(d * d, sizeof(double));
    if (!L) { free(sigma_x); free(sigma_y); return; }
    /* Simple Cholesky: L * L^T = (d+lambda) * P */
    size_t i, j, kk;
    for (i = 0; i < d; i++) {
        for (j = 0; j <= i; j++) {
            double s = (d + lambda) * P[i*d+j];
            size_t mm;
            for (mm = 0; mm < j; mm++)
                s -= L[i*d+mm] * L[j*d+mm];
            if (i == j) {
                L[i*d+i] = (s > 0.0) ? sqrt(s) : 0.0;
            } else {
                L[i*d+j] = (L[j*d+j] > 0.0) ? s / L[j*d+j] : 0.0;
            }
        }
    }
    /* Sigma point 0 = x */
    memcpy(&sigma_x[0], x, d * sizeof(double));
    for (i = 0; i < d; i++) {
        size_t idx_p = (i + 1) * d;
        size_t idx_m = (i + 1 + d) * d;
        for (j = 0; j < d; j++) {
            sigma_x[idx_p + j] = x[j] + L[j*d+i];
            sigma_x[idx_m + j] = x[j] - L[j*d+i];
        }
    }
    /* Predict sigma points through f */
    double *sigma_xp = (double *)malloc(n_sigma * d * sizeof(double));
    if (!sigma_xp) { free(L); free(sigma_x); free(sigma_y); return; }
    size_t s;
    for (s = 0; s < n_sigma; s++)
        model->f(&sigma_x[s*d], u, &sigma_xp[s*d]);
    /* Predicted mean */
    double *x_pred = (double *)calloc(d, sizeof(double));
    if (x_pred) {
        for (j = 0; j < d; j++)
            x_pred[j] = Wm0 * sigma_xp[j];
        for (s = 1; s < n_sigma; s++)
            for (j = 0; j < d; j++)
                x_pred[j] += Wi * sigma_xp[s*d+j];
    }
    /* Predicted covariance */
    double *P_pred = (double *)calloc(d * d, sizeof(double));
    if (P_pred && x_pred) {
        for (j = 0; j < d; j++)
            for (kk = 0; kk < d; kk++)
                P_pred[j*d+kk] = Wc0 * (sigma_xp[j]-x_pred[j])
                                * (sigma_xp[kk]-x_pred[kk]);
        for (s = 1; s < n_sigma; s++)
            for (j = 0; j < d; j++)
                for (kk = 0; kk < d; kk++)
                    P_pred[j*d+kk] += Wi * (sigma_xp[s*d+j]-x_pred[j])
                                     * (sigma_xp[s*d+kk]-x_pred[kk]);
        for (i = 0; i < d * d; i++) P_pred[i] += model->Q[i];
    }
    /* Propagate through observation function */
    for (s = 0; s < n_sigma; s++)
        model->h(&sigma_xp[s*d], &sigma_y[s*d]);
    /* Predicted observation mean */
    double *y_pred = (double *)calloc(m, sizeof(double));
    if (y_pred) {
        for (j = 0; j < m; j++)
            y_pred[j] = Wm0 * sigma_y[j];
        for (s = 1; s < n_sigma; s++)
            for (j = 0; j < m; j++)
                y_pred[j] += Wi * sigma_y[s*d+j];
    }
    /* Innovation covariance Pyy and cross covariance Pxy */
    double *Pyy = (double *)calloc(m * m, sizeof(double));
    double *Pxy = (double *)calloc(d * m, sizeof(double));
    if (Pyy && Pxy && y_pred && x_pred) {
        for (j = 0; j < m; j++)
            for (kk = 0; kk < m; kk++)
                Pyy[j*m+kk] = Wc0 * (sigma_y[j]-y_pred[j])
                             * (sigma_y[kk]-y_pred[kk]);
        for (j = 0; j < d; j++)
            for (kk = 0; kk < m; kk++)
                Pxy[j*m+kk] = Wc0 * (sigma_xp[j]-x_pred[j])
                             * (sigma_y[kk]-y_pred[kk]);
        for (s = 1; s < n_sigma; s++) {
            for (j = 0; j < m; j++)
                for (kk = 0; kk < m; kk++)
                    Pyy[j*m+kk] += Wi * (sigma_y[s*m+j]-y_pred[j])
                                  * (sigma_y[s*m+kk]-y_pred[kk]);
            for (j = 0; j < d; j++)
                for (kk = 0; kk < m; kk++)
                    Pxy[j*m+kk] += Wi * (sigma_xp[s*d+j]-x_pred[j])
                                  * (sigma_y[s*m+kk]-y_pred[kk]);
        }
        for (i = 0; i < m*m; i++) Pyy[i] += model->R[i];
    }
    /* Kalman gain K = Pxy * Pyy^{-1} */
    if (Pyy && Pxy && y_pred && x_pred) {
        double *Pyy_inv = (double *)malloc(m * m * sizeof(double));
        if (Pyy_inv) {
            mat_inverse(Pyy, Pyy_inv, m);
            /* K = Pxy * Pyy_inv */
            double *K = (double *)calloc(d * m, sizeof(double));
            if (K) {
                size_t r, c, k2;
                for (r = 0; r < d; r++)
                    for (c = 0; c < m; c++) {
                        double sval = 0.0;
                        for (k2 = 0; k2 < m; k2++)
                            sval += Pxy[r*m+k2] * Pyy_inv[k2*m+c];
                        K[r*m+c] = sval;
                    }
                /* x_out = x_pred + K*(y - y_pred) */
                double *innov = (double *)malloc(m * sizeof(double));
                if (innov) {
                    for (j = 0; j < m; j++)
                        innov[j] = y[j] - y_pred[j];
                    for (j = 0; j < d; j++) {
                        double k_innov = 0.0;
                        for (kk = 0; kk < m; kk++)
                            k_innov += K[j*m+kk] * innov[kk];
                        x_out[j] = x_pred[j] + k_innov;
                    }
                    free(innov);
                }
                /* P_out = P_pred - K * Pyy * K^T */
                double *KPy = (double *)calloc(d * m, sizeof(double));
                double *KPyyKt = (double *)calloc(d * d, sizeof(double));
                if (KPy && KPyyKt) {
                    size_t r, c, k2;
                    for (r = 0; r < d; r++)
                        for (c = 0; c < m; c++) {
                            double sval = 0.0;
                            for (k2 = 0; k2 < m; k2++)
                                sval += K[r*m+k2] * Pyy[k2*m+c];
                            KPy[r*m+c] = sval;
                        }
                    for (r = 0; r < d; r++)
                        for (c = 0; c < d; c++) {
                            double sval = 0.0;
                            for (k2 = 0; k2 < m; k2++)
                                sval += KPy[r*m+k2] * K[c*m+k2];
                            KPyyKt[r*d+c] = sval;
                        }
                    for (i = 0; i < d*d; i++)
                        P_out[i] = P_pred[i] - KPyyKt[i];
                }
                free(KPy); free(KPyyKt);
                free(K);
            }
            free(Pyy_inv);
        }
    }
    free(L); free(sigma_x); free(sigma_y); free(sigma_xp);
    free(x_pred); free(P_pred); free(y_pred); free(Pyy); free(Pxy);
}

/* ================================================================
 * L5: Particle Filter — Sequential Importance Resampling
 * ================================================================ */

/*
 * particle_filter_sir — bootstrap particle filter (SIR) single step.
 * Reference: Gordon, Salmond, Smith, IEE Proc-F (1993).
 * 1. Propagate each particle through f (with process noise).
 * 2. Compute importance weight w_i = p(y | x_i).
 * 3. Resample with replacement proportional to weights.
 * O(n_particles * state_dim). Independent knowledge point:
 * Monte Carlo approximation of the Bayes filter for arbitrary
 * nonlinear/non-Gaussian systems.
 */
void particle_filter_sir(size_t n_particles, size_t state_dim,
                         size_t obs_dim,
                         prob_real_t *particles,
                         prob_real_t *weights,
                         void (*f)(const prob_real_t*, const prob_real_t*,
                                   prob_real_t*, uint64_t*),
                         prob_real_t (*obs_likelihood)(const prob_real_t*,
                                                       const prob_real_t*,
                                                       size_t),
                         const prob_real_t *u, const prob_real_t *y,
                         uint64_t *seed) {
    if (!particles || !weights || !f || !obs_likelihood
        || n_particles == 0) return;
    /* Propagate particles */
    double *new_p = (double *)malloc(n_particles * state_dim
                                     * sizeof(double));
    double *new_w = (double *)malloc(n_particles * sizeof(double));
    if (!new_p || !new_w) { free(new_p); free(new_w); return; }
    double w_sum = 0.0;
    size_t i;
    for (i = 0; i < n_particles; i++) {
        f(&particles[i * state_dim], u,
          &new_p[i * state_dim], seed);
        new_w[i] = obs_likelihood(&new_p[i * state_dim], y, obs_dim);
        w_sum += new_w[i];
    }
    /* Normalise weights */
    if (w_sum > 0.0)
        for (i = 0; i < n_particles; i++) new_w[i] /= w_sum;
    else {
        /* Degenerate: uniform weights */
        double uni = 1.0 / (double)n_particles;
        for (i = 0; i < n_particles; i++) new_w[i] = uni;
    }
    /* Systematic resampling */
    double *cum_w = (double *)malloc(n_particles * sizeof(double));
    if (!cum_w) { free(new_p); free(new_w); return; }
    cum_w[0] = new_w[0];
    for (i = 1; i < n_particles; i++)
        cum_w[i] = cum_w[i - 1] + new_w[i];
    double u0 = lcg_uniform(seed) / (double)n_particles;
    size_t idx = 0;
    for (i = 0; i < n_particles; i++) {
        double u = u0 + (double)i / (double)n_particles;
        while (idx < n_particles && cum_w[idx] < u) idx++;
        if (idx >= n_particles) idx = n_particles - 1;
        memcpy(&particles[i * state_dim],
               &new_p[idx * state_dim],
               state_dim * sizeof(double));
        weights[i] = new_w[idx];
    }
    free(new_p); free(new_w); free(cum_w);
}

/*
 * particle_filter_estimate — weighted mean of particles.
 * O(n_particles * state_dim).
 */
void particle_filter_estimate(const prob_real_t *particles,
                              const prob_real_t *weights,
                              size_t n_particles, size_t state_dim,
                              prob_real_t *x_est) {
    if (!particles || !weights || !x_est || n_particles == 0) return;
    size_t j;
    for (j = 0; j < state_dim; j++) {
        double s = 0.0;
        size_t i;
        for (i = 0; i < n_particles; i++)
            s += weights[i] * particles[i * state_dim + j];
        x_est[j] = s;
    }
}

/* ================================================================
 * L6: Likelihood and Model Selection
 * ================================================================ */

/*
 * kalman_log_likelihood — log-likelihood of observations under Kalman model.
 * log L = -0.5 * sum_t [ m*log(2*pi) + log|S_t| + z_t^T * S_t^{-1} * z_t ]
 * Useful for parameter estimation (e.g., EM algorithm for unknown Q, R).
 * O(T * d^3). Independent knowledge point: model evaluation via
 * prediction-error decomposition.
 */
prob_real_t kalman_log_likelihood(const kalman_model_t *model,
                                  const prob_real_t *observations,
                                  size_t T,
                                  const prob_real_t *x_init,
                                  const prob_real_t *P_init) {
    if (!model || !observations || !x_init || !P_init || T == 0) return NAN;
    size_t d = model->state_dim, m = model->obs_dim;
    double *x = (double *)malloc(d * sizeof(double));
    double *P = (double *)malloc(d * d * sizeof(double));
    if (!x || !P) { free(x); free(P); return NAN; }
    memcpy(x, x_init, d * sizeof(double));
    memcpy(P, P_init, d * d * sizeof(double));
    double logL = 0.0;
    double log_2pi = log(2.0 * M_PI);
    size_t t;
    double *x_pred = (double *)malloc(d * sizeof(double));
    double *P_pred = (double *)malloc(d * d * sizeof(double));
    if (!x_pred || !P_pred) {
        free(x); free(P); free(x_pred); free(P_pred); return NAN;
    }
    for (t = 0; t < T; t++) {
        kalman_predict(model, x, P, NULL, x_pred, P_pred);
        /* Innovation cov S = H*P_pred*H^T + R */
        double *Ht = (double *)malloc(d * m * sizeof(double));
        double *PHt = (double *)malloc(d * m * sizeof(double));
        double *S = (double *)malloc(m * m * sizeof(double));
        if (!Ht || !PHt || !S) {
            free(Ht); free(PHt); free(S); break;
        }
        mat_transpose(model->H, Ht, m, d);
        size_t i, j, kk;
        for (i = 0; i < d; i++)
            for (j = 0; j < m; j++) {
                double sval = 0.0;
                for (kk = 0; kk < d; kk++)
                    sval += P_pred[i*d+kk] * Ht[kk*m+j];
                PHt[i*m+j] = sval;
            }
        for (i = 0; i < m; i++)
            for (j = 0; j < m; j++) {
                double sval = 0.0;
                for (kk = 0; kk < d; kk++)
                    sval += model->H[i*d+kk] * PHt[kk*m+j];
                S[i*m+j] = sval + model->R[i*m+j];
            }
        /* log det |S| */
        double log_det_S = 0.0;
        for (i = 0; i < m; i++) log_det_S += log(fabs(S[i*m+i]));
        /* z = y - H*x_pred */
        double *z = (double *)malloc(m * sizeof(double));
        if (z) {
            memcpy(z, &observations[t*m], m * sizeof(double));
            double *Hx = (double *)calloc(m, sizeof(double));
            if (Hx) {
                mat_vec_mul(model->H, x_pred, Hx, m, d);
                for (i = 0; i < m; i++) z[i] -= Hx[i];
                free(Hx);
            }
            /* z^T * S^{-1} * z — solve S*sol = z */
            double *Sc = (double *)malloc(m * m * sizeof(double));
            double *zc = (double *)malloc(m * sizeof(double));
            if (Sc && zc) {
                memcpy(Sc, S, m * m * sizeof(double));
                memcpy(zc, z, m * sizeof(double));
                double *sol = (double *)malloc(m * sizeof(double));
                if (sol) {
                    int ret = mat_gaussian_solve(Sc, zc, m, sol);
                    if (ret == 0) {
                        double quad = 0.0;
                        for (i = 0; i < m; i++)
                            quad += z[i] * sol[i];
                        logL += -0.5 * ((double)m * log_2pi
                                        + log_det_S + quad);
                    }
                    free(sol);
                }
            }
            free(Sc); free(zc);
            free(z);
        }
        /* Update for next iteration */
        kalman_update(model, x_pred, P_pred,
                      &observations[t*m], x, P, NULL);
        free(Ht); free(PHt); free(S);
    }
    free(x); free(P); free(x_pred); free(P_pred);
    return logL;
}

/*
 * info_criterion_aic — AIC = 2*k - 2*log(L).
 * Lower is better. Penalty 2*k conservatively estimates model
 * complexity (Akaike, 1974).
 * O(1). Independent knowledge point: AIC for model selection.
 */
prob_real_t info_criterion_aic(prob_real_t log_likelihood, size_t n_params) {
    return 2.0 * (double)n_params - 2.0 * log_likelihood;
}

/*
 * info_criterion_bic — BIC = k*ln(n) - 2*ln(L).
 * Heavier penalty than AIC; asymptotically consistent (Schwarz, 1978).
 * O(1). Independent knowledge point: BIC / Schwarz criterion.
 */
prob_real_t info_criterion_bic(prob_real_t log_likelihood, size_t n_params,
                               size_t n_observations) {
    return (double)n_params * log((double)n_observations)
           - 2.0 * log_likelihood;
}
