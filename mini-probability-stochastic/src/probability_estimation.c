/**
 * @file probability_estimation.c
 * @brief Maximum likelihood estimation, KS goodness-of-fit test,
 *        empirical distribution with Welford-Terriberry online algorithm,
 *        quantile estimation (Hyndman-Fan Type 7), and bootstrap CI.
 *
 * Knowledge: L6 engineering problems (MLE, KS test, empirical distribution,
 * bootstrap), L3 engineering quantities (online moments).
 */
#include "probability_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int cmp_double(const void *a, const void *b) {
    double d = *(const double *)a - *(const double *)b;
    return (d > 0.0) ? 1 : ((d < 0.0) ? -1 : 0);
}

/*
 * prob_mle_gaussian — MLE for N(mu, sigma^2).
 * mu_hat = (1/n) * sum x_i.  sigma2_hat = (1/n) * sum (x_i - mu_hat)^2.
 * Two-pass algorithm for numerical stability (avoids catastrophic
 * cancellation in single-pass variance).
 * Knowledge: Gaussian MLE is the sample mean and biased sample variance.
 * The bias E[sigma2_hat] = ((n-1)/n) * sigma^2 motivates Bessel correction.
 * O(n). Independent knowledge point: MLE for the most important distribution.
 */
void prob_mle_gaussian(const prob_real_t *data, size_t n,
                       prob_real_t *mu_hat, prob_real_t *sigma2_hat) {
    if (!data || n == 0 || (!mu_hat && !sigma2_hat)) return;
    double sum = 0.0;
    size_t i;
    for (i = 0; i < n; i++) sum += data[i];
    double mu = sum / (double)n;
    double sum_sq = 0.0;
    for (i = 0; i < n; i++) {
        double d = data[i] - mu;
        sum_sq += d * d;
    }
    if (mu_hat) *mu_hat = mu;
    if (sigma2_hat) *sigma2_hat = sum_sq / (double)n;
}

/*
 * prob_mle_exponential — MLE for Exp(lambda): lambda_hat = n / sum(x_i).
 * Knowledge: exponential MLE is the reciprocal sample mean.
 * The exponential family has sufficient statistic T = sum x_i.
 * Returns NaN if any data is negative or sum is zero.
 * O(n). Independent knowledge point: MLE for lifetime/failure-time models.
 */
prob_real_t prob_mle_exponential(const prob_real_t *data, size_t n) {
    if (!data || n == 0) return NAN;
    double sum = 0.0;
    size_t i;
    for (i = 0; i < n; i++) {
        if (data[i] < 0.0) return NAN;
        sum += data[i];
    }
    return (sum > 0.0) ? (double)n / sum : NAN;
}

/*
 * prob_mle_poisson — MLE for Poisson(lambda): lambda_hat = sample mean.
 * Knowledge: Poisson MLE is the sample mean. It is unbiased and efficient
 * (achieves Cramer-Rao lower bound). Sufficient statistic: sum x_i.
 * O(n). Independent knowledge point: MLE for count data models.
 */
prob_real_t prob_mle_poisson(const int32_t *data, size_t n) {
    if (!data || n == 0) return NAN;
    double sum = 0.0;
    size_t i;
    for (i = 0; i < n; i++) sum += (double)data[i];
    return sum / (double)n;
}

/*
 * prob_ks_test — Kolmogorov-Smirnov statistic D_n = sup_x |F_n(x) - F(x)|.
 * Tests null hypothesis that data follow the specified continuous distribution.
 * Distribution-free under the null: the distribution of D_n does not depend on F.
 * Reference: Massey, JASA (1951).
 * O(n log n) due to sorting. Independent knowledge point: nonparametric
 * goodness-of-fit testing.
 */
prob_real_t prob_ks_test(const prob_real_t *data, size_t n,
                         const prob_dist_cont_t *dist) {
    if (!data || n == 0 || !dist) return NAN;
    double *sorted = (double *)malloc(n * sizeof(double));
    if (!sorted) return NAN;
    memcpy(sorted, data, n * sizeof(double));
    qsort(sorted, n, sizeof(double), cmp_double);
    double sup_diff = 0.0;
    size_t i;
    for (i = 0; i < n; i++) {
        double Fn_plus  = (double)(i + 1) / (double)n;
        double Fn_minus = (double)i / (double)n;
        double F0 = prob_cdf_cont(dist, sorted[i]);
        double d1 = fabs(Fn_plus - F0);
        double d2 = fabs(Fn_minus - F0);
        if (d1 > sup_diff) sup_diff = d1;
        if (d2 > sup_diff) sup_diff = d2;
    }
    free(sorted);
    return sup_diff;
}

/* ================================================================
 * L6: Empirical Distribution with Welford's Online Algorithm
 * ================================================================ */

/*
 * prob_empirical_create — allocate empirical distribution buffer.
 * O(1). Independent knowledge point: incremental statistics data structure.
 */
prob_empirical_t *prob_empirical_create(size_t capacity) {
    if (capacity == 0) return NULL;
    prob_empirical_t *emp = (prob_empirical_t *)malloc(
        sizeof(prob_empirical_t));
    if (!emp) return NULL;
    emp->samples = (double *)malloc(capacity * sizeof(double));
    if (!emp->samples) { free(emp); return NULL; }
    emp->n = 0; emp->cap = capacity; emp->sum = 0.0;
    emp->mean = 0.0; emp->var = 0.0; emp->skew = 0.0; emp->kurt = 0.0;
    emp->m2 = 0.0; emp->m3 = 0.0; emp->m4 = 0.0;
    return emp;
}

/*
 * prob_empirical_insert — insert one sample with Welford-Terriberry update.
 * Knowledge: online moment computation (Welford, 1962; Terriberry, 2007).
 * Avoids catastrophic cancellation inherent in the textbook two-pass formula.
 * Stores only O(1) accumulators regardless of n: the central moments m_2,
 * m_3, m_4 (sums of powers of deviations).
 * O(1) amortised (realloc on capacity doubling). Space: O(capacity).
 * Independent knowledge point: numerically stable streaming statistics.
 */
void prob_empirical_insert(prob_empirical_t *emp, prob_real_t x) {
    if (!emp || !emp->samples) return;
    if (emp->n >= emp->cap) {
        size_t new_cap = emp->cap * 2;
        double *ns = (double *)realloc(emp->samples,
                                        new_cap * sizeof(double));
        if (!ns) return;
        emp->samples = ns; emp->cap = new_cap;
    }
    emp->samples[emp->n] = x;
    size_t n1 = emp->n;
    emp->n++;
    double n2 = (double)emp->n;
    double delta = x - emp->mean;
    double delta_n = delta / n2;
    double term1 = delta * delta_n * (double)n1;
    emp->mean += delta_n;
    /* Terriberry (2007) extension of Welford to m3 and m4 */
    emp->m4 += term1 * delta_n * delta_n * (n2 * n2 - 3.0 * n2 + 3.0)
             + 6.0 * delta_n * delta_n * emp->m2
             - 4.0 * delta_n * emp->m3;
    emp->m3 += term1 * delta_n * (n2 - 2.0) - 3.0 * delta_n * emp->m2;
    emp->m2 += term1;
    if (emp->n > 1) emp->var = emp->m2 / (n2 - 1.0);
    if (emp->n > 2) emp->skew = sqrt(n2) * emp->m3 / pow(emp->m2, 1.5);
    if (emp->n > 3) emp->kurt = n2 * emp->m4 / (emp->m2 * emp->m2) - 3.0;
    emp->sum += x;
}

/*
 * prob_empirical_free — release empirical distribution. O(1).
 */
void prob_empirical_free(prob_empirical_t *emp) {
    if (emp) { free(emp->samples); free(emp); }
}

/*
 * prob_empirical_quantile — empirical quantile via Hyndman-Fan Type 7
 * (linear interpolation between order statistics).
 * This is the default quantile definition in R, S-PLUS, NumPy, and Julia.
 * Given n sorted values x_(1) <= ... <= x_(n), the p-quantile is:
 *   Q(p) = (1 - gamma) * x_(j) + gamma * x_(j+1)
 * where j = floor((n-1)*p + 1), gamma = (n-1)*p + 1 - j.
 * Knowledge: Type 7 is optimal for continuous distributions under
 * the L2 loss criterion (Hyndman & Fan, 1996, The American Statistician).
 * O(n log n) for the sort; O(1) lookup.
 * Independent knowledge point: nonparametric quantile estimation.
 */
prob_real_t prob_empirical_quantile(const prob_empirical_t *emp, prob_t p) {
    if (!emp || emp->n == 0 || p < 0.0 || p > 1.0) return NAN;
    double *sorted = (double *)malloc(emp->n * sizeof(double));
    if (!sorted) return NAN;
    memcpy(sorted, emp->samples, emp->n * sizeof(double));
    qsort(sorted, emp->n, sizeof(double), cmp_double);
    double h = (emp->n - 1) * p;
    size_t lo = (size_t)h;
    double gamma = h - (double)lo;
    double q;
    if (lo + 1 >= emp->n)
        q = sorted[emp->n - 1];
    else
        q = (1.0 - gamma) * sorted[lo] + gamma * sorted[lo + 1];
    free(sorted);
    return q;
}

/* ================================================================
 * L7: Bootstrap Confidence Interval (Efron, 1979)
 * ================================================================ */

/*
 * prob_bootstrap_ci — percentile bootstrap CI for the mean.
 * Resamples with replacement B times; returns (lower, upper) bounds.
 * Knowledge: bootstrap is a nonparametric method for estimating sampling
 * distributions. Does not require normality or closed-form variance.
 * The percentile method is the simplest bootstrap CI.
 * Reference: Efron & Tibshirani, "An Introduction to the Bootstrap" (1993).
 * O(B * n). Independent knowledge point: resampling-based inference.
 */
void prob_bootstrap_ci(const prob_real_t *data, size_t n, size_t B,
                       prob_t confidence,
                       prob_real_t *lower, prob_real_t *upper) {
    if (!data || n == 0 || B == 0 || !lower || !upper) return;
    double *boot_means = (double *)malloc(B * sizeof(double));
    if (!boot_means) return;
    /* Simple 64-bit LCG for deterministic reproducible resampling */
    uint64_t seed = 12345ULL;
    size_t b;
    for (b = 0; b < B; b++) {
        double sum_b = 0.0;
        size_t i;
        for (i = 0; i < n; i++) {
            seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
            size_t idx = (size_t)(((double)(seed >> 32)
                          / 4294967296.0) * (double)n);
            if (idx >= n) idx = n - 1;
            sum_b += data[idx];
        }
        boot_means[b] = sum_b / (double)n;
    }
    qsort(boot_means, B, sizeof(double), cmp_double);
    double alpha = (1.0 - confidence) / 2.0;
    size_t idx_lo = (size_t)(alpha * (double)B);
    size_t idx_hi = (size_t)((1.0 - alpha) * (double)B);
    if (idx_lo >= B) idx_lo = 0;
    if (idx_hi >= B) idx_hi = B - 1;
    *lower = boot_means[idx_lo];
    *upper = boot_means[idx_hi];
    free(boot_means);
}
