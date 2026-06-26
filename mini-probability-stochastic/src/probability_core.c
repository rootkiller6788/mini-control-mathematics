/**
 * @file probability_core.c
 * @brief Shared numeric infrastructure: PRNG (LCG, xorshift, Box-Muller,
 *        Marsaglia polar), direct distribution sampling (inverse transform,
 *        Marsaglia-Tsang gamma, beta via gamma ratio), and special
 *        mathematical functions (erf, lgamma, regularised incomplete
 *        gamma and beta).
 *
 * Knowledge: L1 definitions (sampling from standard distributions),
 * L5 engineering methods (PRNG algorithms, inverse transform,
 * accept-reject, special functions).
 */
#include "probability_core.h"
#include <math.h>
#include <stdlib.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

static double sq(double x) { return x * x; }

/* ================================================================
 * L5: Pseudo-Random Number Generators
 * ================================================================ */

/*
 * prob_lcg64_uniform — 64-bit LCG (MMIX parameters, Knuth).
 * Returns uniform double in (0,1). Period ~2^64. O(1).
 * Independent knowledge point: LCG is the simplest widely-used PRNG.
 */
double prob_lcg64_uniform(uint64_t *state) {
    *state = *state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)((*state) >> 32) / 4294967296.0;
}

/*
 * prob_xorshift64_uniform — Marsaglia's xorshift64 PRNG.
 * Returns uniform double in (0,1). Period 2^64-1. O(1).
 * Independent knowledge point: nonlinear shift-register PRNG; faster
 * than LCG on many architectures.
 */
double prob_xorshift64_uniform(uint64_t *state) {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    return (double)(*state >> 11) * (1.0 / 9007199254740992.0);
}

/*
 * prob_box_muller — Box-Muller transform for N(0,1).
 * sqrt(-2*ln(U1)) * cos(2*pi*U2). O(1).
 * Independent knowledge point: fundamental algorithm for Gaussian
 * random variate generation from uniform inputs (Box & Muller, 1958).
 */
double prob_box_muller(uint64_t *state) {
    double u1 = prob_lcg64_uniform(state);
    double u2 = prob_lcg64_uniform(state);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/*
 * prob_marsaglia_polar — Marsaglia polar method for N(0,1).
 * Avoids trigonometric function calls. O(1) expected.
 * Independent knowledge point: polar rejection method generates
 * Gaussian pairs without sin/cos (Marsaglia & Bray, 1964).
 */
double prob_marsaglia_polar(uint64_t *state) {
    double u1, u2, s;
    do {
        u1 = 2.0 * prob_lcg64_uniform(state) - 1.0;
        u2 = 2.0 * prob_lcg64_uniform(state) - 1.0;
        s = u1 * u1 + u2 * u2;
    } while (s >= 1.0 || s == 0.0);
    return u1 * sqrt(-2.0 * log(s) / s);
}

/* ================================================================
 * L2: Direct Sampling from Standard Distributions
 * ================================================================ */

/*
 * prob_sample_uniform — X ~ Uniform(a, b). Inverse CDF. O(1).
 * Knowledge: the most basic sampling scheme; maps (0,1) to (a,b).
 */
double prob_sample_uniform(double a, double b, uint64_t *state) {
    return a + (b - a) * prob_lcg64_uniform(state);
}

/*
 * prob_sample_exponential — X ~ Exp(lambda) via inverse CDF.
 * F^{-1}(p) = -ln(1-p)/lambda. O(1).
 * Knowledge: inverse transform sampling (Devroye, 1986, Ch. 2).
 * Any continuous distribution can be sampled as F^{-1}(U).
 * Exponential is the canonical example with closed-form inverse CDF.
 */
double prob_sample_exponential(double lambda, uint64_t *state) {
    if (lambda <= 0.0) return NAN;
    return -log(1.0 - prob_lcg64_uniform(state)) / lambda;
}

/*
 * prob_sample_gaussian — X ~ N(mu, sigma^2) via Box-Muller. O(1).
 * Knowledge: location-scale transformation: if Z ~ N(0,1),
 * then mu + sigma*Z ~ N(mu, sigma^2).
 */
double prob_sample_gaussian(double mu, double sigma, uint64_t *state) {
    return mu + sigma * prob_box_muller(state);
}

/*
 * prob_sample_bernoulli — X ~ Bernoulli(p). O(1).
 * Knowledge: fundamental binary outcome model; building block
 * for binomial, geometric, and negative binomial.
 */
int prob_sample_bernoulli(double p, uint64_t *state) {
    return (prob_lcg64_uniform(state) < p) ? 1 : 0;
}

/*
 * prob_sample_binomial — X ~ Binomial(n, p) via sum of n Bernoulli trials.
 * O(n). For large n, use normal approximation or BTPE algorithm (Kachitvichyanukul).
 * Knowledge: sum of n i.i.d. Bernoulli(p) is Binomial(n,p) —
 * reproductive property of the Bernoulli family.
 */
int prob_sample_binomial(int n, double p, uint64_t *state) {
    if (n <= 0 || p <= 0.0) return 0;
    if (p >= 1.0) return n;
    int count = 0, i;
    for (i = 0; i < n; i++)
        if (prob_lcg64_uniform(state) < p) count++;
    return count;
}

/*
 * prob_sample_poisson — X ~ Poisson(lambda) via Knuth's method.
 * Uses the Poisson-exponential relationship: the number of
 * Exp(1) events needed to exceed lambda follows Poisson(lambda).
 * O(lambda) expected. Knowledge: duality between Poisson counts
 * and exponential inter-arrival times.
 */
int prob_sample_poisson(double lambda, uint64_t *state) {
    if (lambda <= 0.0) return 0;
    double L = exp(-lambda), p = 1.0;
    int k = 0;
    do { k++; p *= prob_lcg64_uniform(state); } while (p > L);
    return k - 1;
}

/*
 * prob_sample_gamma_marsaglia — Marsaglia-Tsang (ACM TOMS, 2000).
 * For alpha >= 1: acceptance-rejection with transformed Gaussian.
 * For alpha < 1: Gamma(alpha) ~ Gamma(alpha+1) * U^{1/alpha}.
 * O(1) expected. Knowledge: efficient gamma sampling without
 * computing log or exp in the inner loop.
 */
double prob_sample_gamma_marsaglia(double alpha, uint64_t *state) {
    if (alpha < 1.0) {
        double g = prob_sample_gamma_marsaglia(alpha + 1.0, state);
        return g * pow(prob_lcg64_uniform(state), 1.0 / alpha);
    }
    double d = alpha - 1.0 / 3.0;
    double c = 1.0 / sqrt(9.0 * d);
    for (;;) {
        double x, v;
        do { x = prob_box_muller(state); v = 1.0 + c * x; }
        while (v <= 0.0);
        v = v * v * v;
        double u = prob_lcg64_uniform(state);
        if (u < 1.0 - 0.0331 * (x * x) * (x * x)) return d * v;
        if (log(u) < 0.5 * x * x + d * (1.0 - v + log(v))) return d * v;
    }
}

/*
 * prob_sample_gamma — X ~ Gamma(alpha, beta) with scale beta.
 * Uses Marsaglia-Tsang. O(1).
 */
double prob_sample_gamma(double alpha, double beta, uint64_t *state) {
    return beta * prob_sample_gamma_marsaglia(alpha, state);
}

/*
 * prob_sample_beta — X ~ Beta(a, b) via gamma ratio.
 * If X ~ Gamma(a,1) and Y ~ Gamma(b,1) independent, X/(X+Y) ~ Beta(a,b).
 * O(1). Knowledge: the gamma-beta relationship is the foundation
 * of Dirichlet processes, Bayesian Beta-binomial conjugacy, and
 * order-statistic distributions.
 */
double prob_sample_beta(double a, double b, uint64_t *state) {
    double x = prob_sample_gamma_marsaglia(a, state);
    double y = prob_sample_gamma_marsaglia(b, state);
    return x / (x + y);
}

/*
 * prob_sample_discrete — sample from discrete distribution via CDF. O(n).
 * Knowledge: general discrete sampling; for repeated draws use
 * the alias method (Walker, 1977) which gives O(1) after O(n) setup.
 */
size_t prob_sample_discrete(const double *probs, size_t n, uint64_t *state) {
    if (!probs || n == 0) return 0;
    double r = prob_lcg64_uniform(state), cum = 0.0;
    size_t i;
    for (i = 0; i < n; i++) {
        cum += probs[i];
        if (r <= cum) return i;
    }
    return n - 1;
}

/* ================================================================
 * L1/L3: Special Mathematical Functions
 * ================================================================ */

/*
 * prob_erf_approx — Gauss error function (Abramowitz & Stegun 7.1.26).
 * Rational approximation with max absolute error 1.5e-7. O(1).
 * Knowledge: erf is the CDF of the half-normal distribution scaled
 * by sqrt(2). It appears in diffusion solutions, Gaussian tail
 * probabilities, and heat equation problems.
 */
double prob_erf_approx(double x) {
    const double a1 = 0.254829592, a2 = -0.284496736;
    const double a3 = 1.421413741, a4 = -1.453152027;
    const double a5 = 1.061405429, p  = 0.3275911;
    double sign = (x >= 0.0) ? 1.0 : -1.0;
    x = fabs(x);
    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - ((((a5 * t + a4) * t + a3) * t + a2) * t + a1)
                    * t * exp(-x * x);
    return sign * y;
}

/*
 * prob_std_normal_cdf — Phi(x) = 0.5 * (1 + erf(x/sqrt(2))). O(1).
 * Knowledge: the standard normal CDF is the cumulative probability
 * function for N(0,1). Not expressible in elementary functions.
 */
double prob_std_normal_cdf(double x) {
    return 0.5 * (1.0 + prob_erf_approx(x * M_SQRT1_2));
}

/*
 * prob_lgamma — log Gamma function via Lanczos approximation.
 * Precision ~1e-14 for all x > 0. Uses reflection for x < 0.5.
 * Reference: Press et al., Numerical Recipes, Ch. 6.1. O(1).
 * Knowledge: lgamma is numerically stable for large arguments
 * compared to tgamma(). It is central to computing binomial
 * coefficients, beta function values, and gamma PDF normalisation.
 */
double prob_lgamma(double x) {
    if (x < 0.5)
        return log(M_PI / sin(M_PI * x)) - prob_lgamma(1.0 - x);
    double c[6] = {76.18009172947146, -86.50532032941677,
                   24.01409824083091, -1.231739572450155,
                   0.1208650973866179e-2, -0.5395239384953e-5};
    double y = x, tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    double ser = 1.000000000190015;
    int j;
    for (j = 0; j < 6; j++) { y += 1.0; ser += c[j] / y; }
    return -tmp + log(2.5066282746310005 * ser / x);
}

/*
 * prob_gamma_cdf_reg — P(a,x), regularised lower incomplete gamma.
 * Series expansion for x < a+1; continued fraction for x >= a+1.
 * Reference: Press et al., NR Ch. 6.2. O(iterations <= 200).
 * Knowledge: foundation for chi-squared, gamma, and Erlang CDFs.
 */
double prob_gamma_cdf_reg(double a, double x) {
    if (x <= 0.0) return 0.0;
    if (x < a + 1.0) {
        double ap = a, sum = 1.0 / a, del = sum;
        int i;
        for (i = 1; i < 200; i++) {
            ap += 1.0; del *= x / ap; sum += del;
            if (fabs(del) < fabs(sum) * 1e-15) break;
        }
        return sum * exp(-x + a * log(x) - prob_lgamma(a));
    } else {
        double b = x + 1.0 - a, c = 1.0 / 1e-30;
        double d = 1.0 / b, h = d;
        int i;
        for (i = 1; i <= 200; i++) {
            double an = -i * (i - a); b += 2.0;
            d = an * d + b;
            if (fabs(d) < 1e-30) d = 1e-30;
            c = b + an / c;
            if (fabs(c) < 1e-30) c = 1e-30;
            d = 1.0 / d; h *= d * c;
            if (fabs(d * c - 1.0) < 1e-15) break;
        }
        return 1.0 - exp(-x + a * log(x) - prob_lgamma(a)) * h;
    }
}

/*
 * prob_beta_cdf_reg — I_x(a,b), regularised incomplete beta function.
 * Continued fraction via modified Lentz method.
 * Reference: Press et al., NR Ch. 6.4. O(iterations <= 200).
 * Knowledge: foundation for beta, binomial tail, Student-t,
 * and Fisher-F distribution CDFs.
 */
double prob_beta_cdf_reg(double a, double b, double x) {
    if (x <= 0.0) return 0.0; if (x >= 1.0) return 1.0;
    const int MAXIT = 200;
    const double EPS = 3.0e-14;
    double qab = a + b, qap = a + 1.0, qam = a - 1.0;
    double c = 1.0, d = 1.0 - qab * x / qap;
    if (fabs(d) < 1e-30) d = 1e-30;
    d = 1.0 / d; double h = d;
    int m;
    for (m = 1; m <= MAXIT; m++) {
        int m2 = 2 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d; if (fabs(d) < 1e-30) d = 1e-30;
        c = 1.0 + aa / c; if (fabs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d; h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d; if (fabs(d) < 1e-30) d = 1e-30;
        c = 1.0 + aa / c; if (fabs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        double del = d * c; h *= del;
        if (fabs(del - 1.0) < EPS) break;
    }
    double bt = exp(prob_lgamma(a + b) - prob_lgamma(a) - prob_lgamma(b)
                    + a * log(x) + b * log(1.0 - x));
    if (x < (a + 1.0) / (a + b + 2.0))
        return bt * h / a;
    else {
        /* I_x(a,b) = 1 - I_{1-x}(b,a) — symmetric property */
        double x2 = 1.0 - x;
        double qab2 = b + a, qap2 = b + 1.0, qam2 = b - 1.0;
        double c2 = 1.0, d2 = 1.0 - qab2 * x2 / qap2;
        if (fabs(d2) < 1e-30) d2 = 1e-30;
        d2 = 1.0 / d2; double h2 = d2;
        int m2;
        for (m2 = 1; m2 <= MAXIT; m2++) {
            int m22 = 2 * m2;
            double aa2 = m2 * (a - m2) * x2 / ((qam2 + m22) * (b + m22));
            d2 = 1.0 + aa2 * d2; if (fabs(d2) < 1e-30) d2 = 1e-30;
            c2 = 1.0 + aa2 / c2; if (fabs(c2) < 1e-30) c2 = 1e-30;
            d2 = 1.0 / d2; h2 *= d2 * c2;
            aa2 = -(b + m2) * (qab2 + m2) * x2 / ((b + m22) * (qap2 + m22));
            d2 = 1.0 + aa2 * d2; if (fabs(d2) < 1e-30) d2 = 1e-30;
            c2 = 1.0 + aa2 / c2; if (fabs(c2) < 1e-30) c2 = 1e-30;
            d2 = 1.0 / d2;
            double del2 = d2 * c2; h2 *= del2;
            if (fabs(del2 - 1.0) < EPS) break;
        }
        double bt2 = exp(prob_lgamma(b + a) - prob_lgamma(b) - prob_lgamma(a)
                         + b * log(x2) + a * log(1.0 - x2));
        return 1.0 - bt2 * h2 / b;
    }
}