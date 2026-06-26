/**
 * @file probability_stats.c
 * @brief Statistical summaries: moments, covariance, correlation,
 *        LLN verification, CLT test, total probability/expectation checks.
 * Knowledge: L3 engineering quantities, L4 conservation laws.
 */
#include "probability_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double sq(double x) { return x * x; }
static double cu(double x) { return x * x * x; }

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

static double prob_erf(double x) {
    const double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
    const double a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
    double s = (x >= 0.0) ? 1.0 : -1.0;
    x = fabs(x);
    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - ((((a5 * t + a4) * t + a3) * t + a2) * t + a1)
                    * t * exp(-x * x);
    return s * y;
}

static double std_norm_cdf(double x) {
    return 0.5 * (1.0 + prob_erf(x * M_SQRT1_2));
}

static int cmp_double(const void *a, const void *b) {
    double d = *(const double *)a - *(const double *)b;
    return (d > 0.0) ? 1 : ((d < 0.0) ? -1 : 0);
}

/* ================================================================
 * L3: Engineering Quantities — Moments and Dependence
 * ================================================================ */

/*
 * prob_mean_disc — E[X] from cached value. O(1).
 * Knowledge: first raw moment / centre of mass of distribution.
 */
prob_real_t prob_mean_disc(const prob_dist_disc_t *dist) {
    return dist ? dist->mean : NAN;
}

/*
 * prob_var_disc — Var[X] = E[(X-mu)^2] from cached value. O(1).
 * Knowledge: second central moment; measures dispersion.
 */
prob_real_t prob_var_disc(const prob_dist_disc_t *dist) {
    return dist ? dist->var : NAN;
}

/*
 * prob_std_disc — sigma = sqrt(Var[X]). O(1).
 * Knowledge: standard deviation has same units as X.
 */
prob_real_t prob_std_disc(const prob_dist_disc_t *dist) {
    return (dist && dist->var >= 0.0) ? sqrt(dist->var) : NAN;
}

/*
 * prob_skewness_disc — gamma_1 = E[((X-mu)/sigma)^3].
 * gamma_1 > 0: right-skewed (e.g. exponential). =0: symmetric. <0: left-skewed.
 * O(n). Independent knowledge: Fisher-Pearson standardised moment coefficient.
 */
prob_real_t prob_skewness_disc(const prob_dist_disc_t *dist) {
    if (!dist || dist->n == 0 || dist->var <= 0.0) return NAN;
    double sigma = sqrt(dist->var), sum_z3 = 0.0;
    size_t i;
    for (i = 0; i < dist->n; i++) {
        double z = (dist->points[i].value - dist->mean) / sigma;
        sum_z3 += cu(z) * dist->points[i].pmass;
    }
    return sum_z3;
}

/*
 * prob_kurtosis_disc — gamma_2 = E[((X-mu)/sigma)^4] - 3.
 * gamma_2 > 0: leptokurtic (heavy tails vs Gaussian). =0: mesokurtic. <0: platykurtic.
 * O(n). Independent knowledge: fourth standardised moment minus 3 (excess kurtosis).
 */
prob_real_t prob_kurtosis_disc(const prob_dist_disc_t *dist) {
    if (!dist || dist->n == 0 || dist->var <= 0.0) return NAN;
    double sigma = sqrt(dist->var), sum_z4 = 0.0;
    size_t i;
    for (i = 0; i < dist->n; i++) {
        double z = (dist->points[i].value - dist->mean) / sigma;
        sum_z4 += z * z * z * z * dist->points[i].pmass;
    }
    return sum_z4 - 3.0;
}

/*
 * prob_covariance — sample covariance Cov[X,Y] with Bessel correction (n-1).
 * Cov[aX+b, cY+d] = ac * Cov[X,Y].  Cov[X,X] = Var[X].
 * O(n). Knowledge: bilinear form measuring joint variability.
 */
prob_real_t prob_covariance(const prob_real_t *x, const prob_real_t *y,
                            size_t n) {
    if (!x || !y || n < 2) return NAN;
    double mx = 0.0, my = 0.0;
    size_t i;
    for (i = 0; i < n; i++) { mx += x[i]; my += y[i]; }
    mx /= (double)n; my /= (double)n;
    double cov = 0.0;
    for (i = 0; i < n; i++)
        cov += (x[i] - mx) * (y[i] - my);
    return cov / (double)(n - 1);
}

/*
 * prob_correlation — Pearson rho = Cov[X,Y] / (sigmaX * sigmaY) in [-1, 1].
 * |rho|=1 iff perfect linear relationship. rho=0 means no linear correlation.
 * O(n). Knowledge: dimensionless measure of linear dependence.
 */
prob_real_t prob_correlation(const prob_real_t *x, const prob_real_t *y,
                             size_t n) {
    if (!x || !y || n < 2) return NAN;
    double mx = 0.0, my = 0.0, vx = 0.0, vy = 0.0, cxy = 0.0;
    size_t i;
    for (i = 0; i < n; i++) { mx += x[i]; my += y[i]; }
    mx /= (double)n; my /= (double)n;
    for (i = 0; i < n; i++) {
        double dx = x[i] - mx, dy = y[i] - my;
        cxy += dx * dy; vx += dx * dx; vy += dy * dy;
    }
    if (vx <= 0.0 || vy <= 0.0) return NAN;
    return cxy / sqrt(vx * vy);
}

/*
 * prob_covariance_matrix — k-by-k covariance matrix from k-by-n data.
 * Output cov[r*k + c] = Cov[X_r, X_c] (Bessel-corrected). Diagonal = variances.
 * O(k^2 * n). Knowledge: multivariate generalisation of variance.
 */
void prob_covariance_matrix(const prob_real_t *data, size_t k, size_t n,
                            prob_real_t *cov) {
    if (!data || !cov || k == 0 || n < 2) return;
    double *means = (double *)malloc(k * sizeof(double));
    if (!means) return;
    size_t r, c, i;
    for (r = 0; r < k; r++) {
        means[r] = 0.0;
        for (i = 0; i < n; i++) means[r] += data[r * n + i];
        means[r] /= (double)n;
    }
    for (r = 0; r < k; r++)
        for (c = 0; c < k; c++) {
            double sum = 0.0;
            for (i = 0; i < n; i++)
                sum += (data[r * n + i] - means[r])
                     * (data[c * n + i] - means[c]);
            cov[r * k + c] = sum / (double)(n - 1);
        }
    free(means);
}

/* ================================================================
 * L4: Conservation Law / Theorem Verification
 * ================================================================ */

/*
 * prob_total_probability_check — verify sum P(A_i) = 1 for a partition.
 * Returns 1.0 - sum. Should be ~0 (within fp epsilon) for valid partition.
 * Knowledge: Kolmogorov axiom — probability measure has total mass 1.
 * O(n). Independent knowledge point: partitioning exhausts the sample space.
 */
prob_real_t prob_total_probability_check(const prob_t *partition, size_t n) {
    if (!partition || n == 0) return 1.0;
    double total = 0.0;
    size_t i;
    for (i = 0; i < n; i++) total += partition[i];
    return 1.0 - total;
}

/*
 * prob_total_expectation_check — verify E[X] = sum_j E[X|B_j] * P(B_j).
 * Returns absolute deviation from the identity.
 * Knowledge: law of total expectation / tower property.
 * E[X] = E[E[X|Y]] for any conditioning variable Y.
 * O(n * m). Independent knowledge point: smoothing property of conditional
 * expectation.
 */
prob_real_t prob_total_expectation_check(const prob_dist_disc_t *dist,
                                         const prob_t *event_probs,
                                         prob_t (**indicators)(prob_real_t),
                                         size_t m) {
    if (!dist || !event_probs || !indicators || m == 0) return NAN;
    double true_mean = prob_mean_disc(dist);
    double total_expect = 0.0;
    size_t j;
    for (j = 0; j < m; j++) {
        double cond_exp = prob_cond_expect(dist, indicators[j],
                                           event_probs[j]);
        total_expect += cond_exp * event_probs[j];
    }
    return fabs(true_mean - total_expect);
}

/*
 * prob_lln_check — verify |sample_mean - true_mean| -> 0 as n grows.
 * Knowledge: weak law of large numbers. Sample mean converges in probability
 * to population mean for i.i.d. data with finite first moment.
 * O(n). Independent knowledge point: convergence of empirical average.
 */
prob_real_t prob_lln_check(const prob_real_t *samples, size_t n,
                           prob_real_t true_mean) {
    if (!samples || n == 0) return NAN;
    double sample_mean = 0.0;
    size_t i;
    for (i = 0; i < n; i++) sample_mean += samples[i];
    sample_mean /= (double)n;
    return fabs(sample_mean - true_mean);
}

/*
 * prob_clt_ks_test — KS statistic comparing standardised sample means to N(0,1).
 * Knowledge: central limit theorem. sqrt(n)*(Xbar_n - mu)/sigma -> N(0,1)
 * in distribution as n -> infinity, under finite variance.
 * The KS statistic D_n = sup_x |F_n(x) - Phi(x)| quantifies the approximation.
 * O(n log n) due to sorting. Independent knowledge point: rate and quality
 * of CLT convergence.
 */
prob_real_t prob_clt_ks_test(const prob_real_t *z, size_t n) {
    if (!z || n == 0) return NAN;
    double *sorted = (double *)malloc(n * sizeof(double));
    if (!sorted) return NAN;
    memcpy(sorted, z, n * sizeof(double));
    qsort(sorted, n, sizeof(double), cmp_double);
    double sup_diff = 0.0;
    size_t i;
    for (i = 0; i < n; i++) {
        double ecdf = (double)(i + 1) / (double)n;
        double tcdf = std_norm_cdf(sorted[i]);
        double d1 = fabs(ecdf - tcdf);
        double d2 = (i > 0) ? fabs((double)i / (double)n - tcdf) : 0.0;
        double diff = (d1 > d2) ? d1 : d2;
        if (diff > sup_diff) sup_diff = diff;
    }
    free(sorted);
    return sup_diff;
}
