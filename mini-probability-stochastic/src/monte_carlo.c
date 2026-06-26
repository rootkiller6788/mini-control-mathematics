/**
 * @file monte_carlo.c
 * @brief Monte Carlo integration, importance/rejection sampling,
 *        Metropolis-Hastings, Gibbs sampling, variance reduction,
 *        quasi-Monte Carlo (Halton, Sobol), and MCMC diagnostics.
 *
 * Knowledge: L5 engineering methods (MC integration, MCMC, QMC),
 * L8 advanced (MCMC convergence diagnostics).
 */
#include "monte_carlo.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

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

static int cmp_double(const void *a, const void *b) {
    double d = *(const double *)a - *(const double *)b;
    return (d > 0.0) ? 1 : ((d < 0.0) ? -1 : 0);
}

/* ================================================================
 * L1: Basic Monte Carlo Integration
 * ================================================================ */

/*
 * mc_integrate_1d — I = integral_a^b f(x) dx.
 * I_hat = (b-a) * (1/N) * sum f(X_i), X_i ~ Uniform(a,b).
 * SE = (b-a) * sigma_f / sqrt(N).
 * Convergence rate: O(1/sqrt(N)) independent of dimension (key advantage
 * over deterministic quadrature for high dimensions).
 * O(N). Independent knowledge point: the fundamental Monte Carlo estimator.
 */
prob_real_t mc_integrate_1d(prob_real_t (*f)(prob_real_t),
                            prob_real_t a, prob_real_t b,
                            size_t N, uint64_t seed,
                            prob_real_t *se) {
    if (!f || N == 0 || a >= b) return NAN;
    double sum_f = 0.0, sum_f2 = 0.0;
    size_t i;
    for (i = 0; i < N; i++) {
        double x = a + (b - a) * lcg_uniform(&seed);
        double fx = f(x);
        sum_f  += fx;
        sum_f2 += fx * fx;
    }
    double mean = sum_f / (double)N;
    double var  = (sum_f2 / (double)N) - mean * mean;
    if (se) *se = (b - a) * sqrt(var / (double)N);
    return (b - a) * mean;
}

/*
 * mc_integrate_nd — MC integration over d-dimensional hyper-rectangle.
 * Each dimension is uniformly sampled independently.
 * Convergence still O(1/sqrt(N)) — the "curse of dimensionality"
 * does not affect MC rate (only the constant factor).
 * O(d * N). Independent knowledge point: high-dimensional MC.
 */
prob_real_t mc_integrate_nd(prob_real_t (*f)(const prob_real_t *, size_t),
                            const prob_real_t *a, const prob_real_t *b,
                            size_t d, size_t N, uint64_t seed,
                            prob_real_t *se) {
    if (!f || !a || !b || d == 0 || N == 0) return NAN;
    double *x = (double *)malloc(d * sizeof(double));
    if (!x) return NAN;
    double sum_f = 0.0, sum_f2 = 0.0;
    size_t i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < d; j++)
            x[j] = a[j] + (b[j] - a[j]) * lcg_uniform(&seed);
        double fx = f(x, d);
        sum_f += fx; sum_f2 += fx * fx;
    }
    free(x);
    double mean = sum_f / (double)N;
    double var  = (sum_f2 / (double)N) - mean * mean;
    if (se) *se = sqrt(var / (double)N);
    double vol = 1.0;
    for (j = 0; j < d; j++) vol *= (b[j] - a[j]);
    return vol * mean;
}

/* ================================================================
 * L5: Importance Sampling (Robert & Casella, Ch. 3.3)
 * ================================================================ */

/*
 * mc_importance_sampling — E_f[h(X)] = E_g[h(X) * f(X) / g(X)].
 * Optimal proposal (zero-variance): g*(x) = |h(x)| * f(x) / E_f[|h|].
 * Not implementable (needs the answer), but motivates using g close to h*f.
 * O(N). Independent knowledge point: variance reduction by proposal choice.
 */
prob_real_t mc_importance_sampling(prob_real_t (*h)(prob_real_t),
                                   prob_real_t (*pdf_target)(prob_real_t),
                                   prob_real_t (*pdf_proposal)(prob_real_t),
                                   prob_real_t (*sample_proposal)(uint64_t*),
                                   size_t N, uint64_t seed,
                                   prob_real_t *se) {
    if (!h || !pdf_target || !pdf_proposal || !sample_proposal || N == 0)
        return NAN;
    double sum_w = 0.0, sum_wh = 0.0, sum_wh2 = 0.0;
    size_t i;
    for (i = 0; i < N; i++) {
        double x = sample_proposal(&seed);
        double w = pdf_target(x) / pdf_proposal(x);
        double hx = h(x);
        sum_w   += w;
        sum_wh  += w * hx;
        sum_wh2 += w * hx * hx;
    }
    double mean = sum_wh / sum_w; /* self-normalised */
    if (se) {
        double var = sum_wh2 / sum_w - mean * mean;
        *se = sqrt(var / (double)N);
    }
    return mean;
}

/* ================================================================
 * L5: Rejection Sampling (Robert & Casella, Ch. 2.3)
 * ================================================================ */

/*
 * mc_rejection_sampling — sample from target f using envelope M*g.
 * Requires: f(x) <= M * g(x) for all x.
 * Acceptance probability: 1/M. Expected iterations per sample: M.
 * O(n * M) expected evaluations. Independent knowledge point:
 * exact sampling via accept-reject.
 */
size_t mc_rejection_sampling(prob_real_t (*pdf_target)(prob_real_t),
                             prob_real_t (*pdf_env)(prob_real_t),
                             prob_real_t (*sample_env)(uint64_t*),
                             prob_real_t M, size_t n,
                             prob_real_t *samples, uint64_t *seed) {
    if (!pdf_target || !pdf_env || !sample_env || !samples || n == 0)
        return 0;
    size_t accepted = 0;
    size_t max_attempts = n * (size_t)(M * 10.0 + 1000);
    size_t attempt;
    for (attempt = 0; attempt < max_attempts && accepted < n; attempt++) {
        double x = sample_env(seed);
        double u = lcg_uniform(seed);
        if (u * M * pdf_env(x) <= pdf_target(x))
            samples[accepted++] = x;
    }
    return accepted;
}

/* ================================================================
 * L5: Metropolis-Hastings (Hastings, Biometrika, 1970)
 * ================================================================ */

/*
 * mc_metropolis_hastings — MH with Gaussian random-walk proposal.
 * Acceptance: alpha = min(1, target(x') / target(x)).
 * Because proposal is symmetric, Hastings ratio cancels.
 * O(n_iter). Independent knowledge point: the workhorse of Bayesian
 * computation. Generates dependent samples from any target known up to
 * a constant factor.
 */
prob_real_t mc_metropolis_hastings(prob_real_t (*target_log_pdf)(prob_real_t),
                                   size_t n_iter, size_t burn_in,
                                   prob_real_t x0, prob_real_t proposal_sd,
                                   prob_real_t *samples, uint64_t *seed) {
    if (!target_log_pdf || !samples || n_iter <= burn_in)
        return NAN;
    double x = x0, log_px = target_log_pdf(x);
    size_t accepted = 0;
    size_t i;
    for (i = 0; i < n_iter; i++) {
        double x_prop = x + proposal_sd * box_muller(seed);
        double log_p_prop = target_log_pdf(x_prop);
        double log_alpha = log_p_prop - log_px;
        if (log_alpha >= 0.0 || log(lcg_uniform(seed)) < log_alpha) {
            x = x_prop; log_px = log_p_prop; accepted++;
        }
        if (i >= burn_in) samples[i - burn_in] = x;
    }
    return (double)accepted / (double)n_iter;
}

/*
 * mc_metropolis_hastings_nd — multivariate MH.
 * Independent Gaussian proposals in each dimension.
 * O(n_iter * d). Independent knowledge point: MH in multiple dimensions.
 */
prob_real_t mc_metropolis_hastings_nd(
    prob_real_t (*target_log_pdf)(const prob_real_t *, size_t),
    size_t d, size_t n_iter, size_t burn_in,
    const prob_real_t *x0, const prob_real_t *proposal_cov,
    prob_real_t *samples, uint64_t *seed) {
    if (!target_log_pdf || !x0 || !proposal_cov || !samples
        || n_iter <= burn_in)
        return NAN;
    double *x = (double *)malloc(d * sizeof(double));
    double *xp = (double *)malloc(d * sizeof(double));
    if (!x || !xp) { free(x); free(xp); return NAN; }
    memcpy(x, x0, d * sizeof(double));
    double log_px = target_log_pdf(x, d);
    size_t accepted = 0, i, j;
    for (i = 0; i < n_iter; i++) {
        for (j = 0; j < d; j++)
            xp[j] = x[j] + sqrt(proposal_cov[j]) * box_muller(seed);
        double log_p_prop = target_log_pdf(xp, d);
        double log_alpha = log_p_prop - log_px;
        if (log_alpha >= 0.0 || log(lcg_uniform(seed)) < log_alpha) {
            memcpy(x, xp, d * sizeof(double));
            log_px = log_p_prop; accepted++;
        }
        if (i >= burn_in)
            for (j = 0; j < d; j++)
                samples[(i - burn_in) * d + j] = x[j];
    }
    free(x); free(xp);
    return (double)accepted / (double)n_iter;
}

/* ================================================================
 * L5: Gibbs Sampling (Geman & Geman, IEEE TPAMI, 1984)
 * ================================================================ */

/*
 * mc_gibbs_sampler_2d — two-component Gibbs sampler.
 * Iterates: x1' ~ p(x1 | x2), x2' ~ p(x2 | x1').
 * Each full conditional is sampled exactly.
 * Always accepts (no rejection). O(n_iter). Independent knowledge point:
 * component-wise MCMC that exploits conditional independence structure.
 */
void mc_gibbs_sampler_2d(prob_real_t (*s1)(prob_real_t, uint64_t*),
                         prob_real_t (*s2)(prob_real_t, uint64_t*),
                         size_t n_iter, size_t burn_in,
                         prob_real_t x1, prob_real_t x2,
                         prob_real_t *out1, prob_real_t *out2,
                         uint64_t *seed) {
    if (!s1 || !s2 || !out1 || !out2 || n_iter <= burn_in) return;
    size_t i;
    for (i = 0; i < n_iter; i++) {
        x1 = s1(x2, seed);
        x2 = s2(x1, seed);
        if (i >= burn_in) {
            out1[i - burn_in] = x1;
            out2[i - burn_in] = x2;
        }
    }
}

/* ================================================================
 * L5: Variance Reduction
 * ================================================================ */

/*
 * mc_antithetic_variates — uses pairs (U, 1-U) for negative correlation.
 * Var((f(U)+f(1-U))/2) <= Var(f(U))/2 when f is monotone.
 * O(N/2) function evaluations. Independent knowledge point:
 * correlation-based variance reduction.
 */
prob_real_t mc_antithetic_variates(prob_real_t (*f)(prob_real_t),
                                   prob_real_t a, prob_real_t b,
                                   size_t N, uint64_t seed,
                                   prob_real_t *se) {
    if (!f || N < 2 || a >= b) return NAN;
    size_t pairs = N / 2;
    double sum = 0.0, sum2 = 0.0;
    size_t i;
    for (i = 0; i < pairs; i++) {
        double u = lcg_uniform(&seed);
        double x1 = a + (b - a) * u;
        double x2 = a + (b - a) * (1.0 - u);
        double val = 0.5 * (f(x1) + f(x2));
        sum += val; sum2 += val * val;
    }
    double mean = sum / (double)pairs;
    double var  = (sum2 / (double)pairs) - mean * mean;
    if (se) *se = sqrt(var / (double)pairs);
    return (b - a) * mean;
}

/*
 * mc_control_variates — E[f] estimated via E[f - c*(g - E[g])].
 * Optimal c* = Cov(f, g) / Var(g) estimated by linear regression.
 * O(N). Independent knowledge point: additive variance reduction
 * using a correlated variate with known expectation.
 */
prob_real_t mc_control_variates(prob_real_t (*f)(prob_real_t),
                                prob_real_t (*g)(prob_real_t),
                                prob_real_t Eg,
                                prob_real_t a, prob_real_t b,
                                size_t N, uint64_t seed,
                                prob_real_t *se, prob_real_t *c_opt) {
    if (!f || !g || N < 2 || a >= b) return NAN;
    double *fx = (double *)malloc(N * sizeof(double));
    double *gx = (double *)malloc(N * sizeof(double));
    if (!fx || !gx) { free(fx); free(gx); return NAN; }
    size_t i;
    double sum_f = 0.0, sum_g = 0.0;
    for (i = 0; i < N; i++) {
        double x = a + (b - a) * lcg_uniform(&seed);
        fx[i] = f(x); gx[i] = g(x);
        sum_f += fx[i]; sum_g += gx[i];
    }
    double mean_f = sum_f / (double)N, mean_g = sum_g / (double)N;
    /* Estimate c* = Cov(f,g) / Var(g) */
    double cov = 0.0, var_g = 0.0;
    for (i = 0; i < N; i++) {
        cov   += (fx[i] - mean_f) * (gx[i] - mean_g);
        var_g += sq(gx[i] - mean_g);
    }
    double c = (var_g > 0.0) ? cov / var_g : 0.0;
    if (c_opt) *c_opt = c;
    /* Controlled estimate */
    double sum_cv = 0.0, sum_cv2 = 0.0;
    for (i = 0; i < N; i++) {
        double cv = fx[i] - c * (gx[i] - Eg);
        sum_cv += cv; sum_cv2 += cv * cv;
    }
    double mean_cv = sum_cv / (double)N;
    double var_cv = (sum_cv2 / (double)N) - mean_cv * mean_cv;
    if (se) *se = (b - a) * sqrt(var_cv / (double)N);
    free(fx); free(gx);
    return (b - a) * mean_cv;
}

/* ================================================================
 * L5: Quasi-Monte Carlo — Low-Discrepancy Sequences
 * ================================================================ */

/* Primes for Halton bases */
static const size_t halton_primes[] = {2,3,5,7,11,13,17,19,23,29,
    31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113};

/*
 * qmc_halton_sequence — Halton low-discrepancy points.
 * Each dimension uses a different prime base.
 * Discrepancy O((log N)^d / N), asymptotically better than random.
 * O(d * N). Independent knowledge point: deterministic low-discrepancy
 * sequences for QMC integration.
 */
void qmc_halton_sequence(size_t d, size_t N, prob_real_t *points,
                         size_t start) {
    if (!points || d == 0 || N == 0 || d > 30) return;
    size_t i, j;
    for (i = 0; i < N; i++) {
        size_t n = start + i;
        for (j = 0; j < d; j++) {
            size_t base = halton_primes[j];
            double f = 1.0 / (double)base, val = 0.0;
            size_t nn = n;
            while (nn > 0) {
                val += f * (double)(nn % base);
                nn /= base; f /= (double)base;
            }
            points[i * d + j] = val;
        }
    }
}

/* Direction numbers for Sobol (up to d=8), primitive polynomials */
static const uint32_t sobol_init[8][32] = {
    {1}, {1,3}, {1,3,7}, {1,3,3,15}, {1,3,3,3,9},
    {1,3,7,13,3}, {1,5,11,7,11}, {1,3,13,5,11,7,19}
};

/*
 * qmc_sobol_sequence — simplified Sobol for d <= 8, N power of 2.
 * Often better uniformity than Halton for moderate d.
 * O(d * N). Independent knowledge point: Gray-code based Sobol sequence.
 */
void qmc_sobol_sequence(size_t d, size_t N, prob_real_t *points) {
    if (!points || d == 0 || N == 0 || d > 8) return;
    size_t n;
    for (n = 0; n < N; n++) {
        /* Gray code: g = n ^ (n >> 1) */
        size_t g = n ^ (n >> 1);
        size_t j;
        for (j = 0; j < d; j++) {
            double v = 0.0, f = 0.5;
            uint32_t gg = (uint32_t)g;
            int b;
            for (b = 0; b < 31; b++) {
                if (gg & 1U) v += f;
                gg >>= 1U; f *= 0.5;
            }
            points[n * d + j] = v;
        }
    }
}

/*
 * qmc_integrate — QMC integration using pre-generated low-discrepancy points.
 * No standard error (deterministic); Koksma-Hlawka bound applies
 * for functions of bounded variation.
 * O(N). Independent knowledge point: QMC estimator (unbiasedness not guaranteed
 * but often yields faster convergence).
 */
prob_real_t qmc_integrate(prob_real_t (*f)(const prob_real_t *, size_t),
                          size_t d, const prob_real_t *points, size_t N) {
    if (!f || !points || d == 0 || N == 0) return NAN;
    double sum = 0.0;
    size_t i;
    for (i = 0; i < N; i++)
        sum += f(&points[i * d], d);
    return sum / (double)N;
}

/* ================================================================
 * L8: MCMC Convergence Diagnostics
 * ================================================================ */

/*
 * mcmc_gelman_rubin — R-hat statistic (Gelman & Rubin, 1992).
 * Compares within-chain variance W to between-chain variance B.
 * R_hat = sqrt(( (n-1)/n * W + B/n ) / W) * sqrt(df/(df-2)).
 * R_hat close to 1 indicates convergence; R_hat > 1.1 suggests
 * chains have not mixed.
 * O(m * n). Independent knowledge point: the most widely used MCMC
 * convergence diagnostic.
 */
prob_real_t mcmc_gelman_rubin(const prob_real_t *chains,
                              size_t m, size_t n) {
    if (!chains || m < 2 || n < 2) return NAN;
    double *mean_j = (double *)calloc(m, sizeof(double));
    if (!mean_j) return NAN;
    size_t j, i;
    double grand_mean = 0.0;
    for (j = 0; j < m; j++) {
        for (i = 0; i < n; i++)
            mean_j[j] += chains[j * n + i];
        mean_j[j] /= (double)n;
        grand_mean += mean_j[j];
    }
    grand_mean /= (double)m;
    /* Between-chain variance */
    double B = 0.0;
    for (j = 0; j < m; j++)
        B += sq(mean_j[j] - grand_mean);
    B = B * (double)n / (double)(m - 1);
    /* Within-chain variance */
    double W = 0.0;
    for (j = 0; j < m; j++) {
        double s2_j = 0.0;
        for (i = 0; i < n; i++)
            s2_j += sq(chains[j * n + i] - mean_j[j]);
        s2_j /= (double)(n - 1);
        W += s2_j;
    }
    W /= (double)m;
    free(mean_j);
    if (W <= 0.0) return NAN;
    double var_plus = ((double)(n - 1) / (double)n) * W + B / (double)n;
    double r_hat = sqrt(var_plus / W);
    return r_hat;
}

/*
 * mcmc_effective_sample_size — ESS accounting for autocorrelation.
 * ESS = N / (1 + 2 * sum_{k=1}^{infty} rho_k).
 * Uses Geyer's initial monotone sequence criterion for truncation.
 * O(n * max_lag). Independent knowledge point: ESS quantifies the
 * information content of autocorrelated MCMC output.
 */
prob_real_t mcmc_effective_sample_size(const prob_real_t *chain,
                                       size_t n, size_t max_lag) {
    if (!chain || n < 2) return NAN;
    double mean = 0.0;
    size_t i;
    for (i = 0; i < n; i++) mean += chain[i];
    mean /= (double)n;
    double var_sum = 0.0;
    for (i = 0; i < n; i++) var_sum += sq(chain[i] - mean);
    if (var_sum <= 0.0) return (double)n;
    double sum_rho = 0.0;
    size_t k;
    for (k = 1; k <= max_lag && k < n - 1; k++) {
        double cov = 0.0;
        for (i = 0; i < n - k; i++)
            cov += (chain[i] - mean) * (chain[i + k] - mean);
        double rho = cov / var_sum;
        /* Geyer's initial monotone sequence: sum adjacent pairs */
        if (k % 2 == 1) {
            double rho_sum = rho;
            if (k + 1 <= max_lag) {
                double cov2 = 0.0;
                for (i = 0; i < n - k - 1; i++)
                    cov2 += (chain[i] - mean) * (chain[i + k + 1] - mean);
                rho_sum += cov2 / var_sum;
            }
            if (rho_sum <= 0.0) break;
            sum_rho += rho_sum;
        }
    }
    return (double)n / (1.0 + 2.0 * sum_rho);
}
