/**
 * @file probability_distributions.c
 * @brief PDF, CDF, inverse CDF for 11 continuous distribution families
 *        plus discrete distribution construction and moment computation.
 *
 * Knowledge: L1 definitions (11 parametric families), L2 core concepts
 * (PMF/PDF/CDF/quantile functions), L6 engineering problems (distribution
 * construction and parameterisation).
 */
#include "probability_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

/* ---- Internal numeric helpers ---- */
static double sq(double x) { return x * x; }
static double cu(double x) { return x * x * x; }

/*
 * erf_approx: Gauss error function.
 * Abramowitz & Stegun 7.1.26 rational approximation.
 * Max absolute error ~ 1.5e-7.  O(1).
 */
static double erf_approx(double x) {
    const double a1 =  0.254829592, a2 = -0.284496736;
    const double a3 =  1.421413741, a4 = -1.453152027;
    const double a5 =  1.061405429, p  =  0.3275911;
    double s = (x >= 0.0) ? 1.0 : -1.0;
    x = fabs(x);
    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - ((((a5 * t + a4) * t + a3) * t + a2) * t + a1)
                    * t * exp(-x * x);
    return s * y;
}

/* Standard normal CDF: Phi(x) = 0.5 * [1 + erf(x/sqrt(2))]. O(1). */
static double std_norm_cdf(double x) {
    return 0.5 * (1.0 + erf_approx(x * M_SQRT1_2));
}

/*
 * lngamma: log Gamma function via Lanczos (Press et al. NR Ch.6.1).
 * Precision ~1e-14.  Uses reflection for x < 0.5.  O(1).
 */
static double lngamma(double x) {
    if (x < 0.5) return log(M_PI / sin(M_PI * x)) - lngamma(1.0 - x);
    double c[6] = { 76.18009172947146, -86.50532032941677,
                    24.01409824083091, -1.231739572450155,
                    0.1208650973866179e-2, -0.5395239384953e-5 };
    double y = x, tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    double ser = 1.000000000190015;
    int j;
    for (j = 0; j < 6; j++) { y += 1.0; ser += c[j] / y; }
    return -tmp + log(2.5066282746310005 * ser / x);
}

/*
 * Regularised incomplete gamma P(a,x).
 * Series for x < a+1, continued fraction otherwise (NR Ch.6.2).
 * O(iterations), typically < 200.
 */
static double gamma_series(double a, double x) {
    if (x <= 0.0) return 0.0;
    double ap = a, sum = 1.0 / a, del = sum;
    int i;
    for (i = 1; i < 200; i++) {
        ap += 1.0; del *= x / ap; sum += del;
        if (fabs(del) < fabs(sum) * 1e-15) break;
    }
    return sum * exp(-x + a * log(x) - lngamma(a));
}

static double gamma_cf(double a, double x) {
    double b = x + 1.0 - a, c = 1.0 / 1e-30, d = 1.0 / b, h = d;
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
    return exp(-x + a * log(x) - lngamma(a)) * h;
}

static double gamma_cdf_reg(double a, double x) {
    if (x <= 0.0) return 0.0;
    return (x < a + 1.0) ? gamma_series(a, x) : 1.0 - gamma_cf(a, x);
}

/*
 * Regularised incomplete beta I_x(a,b) by continued fraction (NR Ch.6.4).
 * O(iterations), typically < 200.
 */
static double betacf(double a, double b, double x) {
    const int MAXIT = 200;
    const double EPS = 3.0e-14;
    double qab = a + b, qap = a + 1.0, qam = a - 1.0;
    double c = 1.0, d = 1.0 - qab * x / qap;
    if (fabs(d) < 1e-30) d = 1e-30;
    d = 1.0 / d;
    double h = d;
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
        double del = d * c;
        h *= del;
        if (fabs(del - 1.0) < EPS) break;
    }
    return h;
}

static double beta_cdf_reg(double a, double b, double x) {
    if (x <= 0.0) return 0.0; if (x >= 1.0) return 1.0;
    double bt = exp(lngamma(a + b) - lngamma(a) - lngamma(b)
                    + a * log(x) + b * log(1.0 - x));
    if (x < (a + 1.0) / (a + b + 2.0))
        return bt * betacf(a, b, x) / a;
    else
        return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

static int cmp_double(const void *a, const void *b) {
    double d = *(const double *)a - *(const double *)b;
    return (d > 0.0) ? 1 : ((d < 0.0) ? -1 : 0);
}

/* ================================================================
 * L2: PMF, CDF, PDF, Inverse CDF
 * ================================================================ */

/*
 * prob_pmf: P(X = x) for discrete distribution.
 * Knowledge: PMF definition — fundamental building block of
 * discrete probability.  O(n) linear scan.
 */
prob_t prob_pmf(const prob_dist_disc_t *dist, prob_real_t x) {
    if (!dist || !dist->points || dist->n == 0) return 0.0;
    size_t i;
    for (i = 0; i < dist->n; i++)
        if (fabs(dist->points[i].value - x) < 1e-12)
            return dist->points[i].pmass;
    return 0.0;
}

/*
 * prob_cdf_disc: F(x) = P(X <= x) = sum_{x_i <= x} p_i.
 * Knowledge: CDF is right-continuous, monotone non-decreasing,
 * F(-inf)=0, F(inf)=1.  O(n).
 */
prob_t prob_cdf_disc(const prob_dist_disc_t *dist, prob_real_t x) {
    if (!dist || !dist->points || dist->n == 0) return 0.0;
    prob_t cdf = 0.0;
    size_t i;
    for (i = 0; i < dist->n; i++)
        if (dist->points[i].value <= x) cdf += dist->points[i].pmass;
    return cdf;
}

/*
 * prob_cdf_cont: closed-form CDF for 11 continuous families.
 * Knowledge: each parametric family has a unique analytic CDF.
 *   Gaussian:  Phi((x-mu)/sigma)
 *   Exponential: 1 - exp(-lambda*x)
 *   Uniform: (x-a)/(b-a)
 *   Gamma: P(alpha, x/beta) — regularised lower incomplete gamma
 *   Beta: I_x(alpha, beta) — regularised incomplete beta
 *   Weibull: 1 - exp(-(x/lambda)^k)
 *   Lognormal: Phi((ln x - mu)/sigma)
 *   Cauchy: 0.5 + arctan((x-x0)/gamma)/pi
 *   Chi-squared: P(k/2, x/2)
 *   Student-t: regularised incomplete beta function
 *   Fisher-F: regularised incomplete beta function
 * O(1) per call for all families.
 */
prob_t prob_cdf_cont(const prob_dist_cont_t *dist, prob_real_t x) {
    if (!dist) return NAN;
    double m = dist->param1, s = dist->param2;
    switch (dist->type) {
    case PROB_DIST_GAUSSIAN:
        if (s <= 0.0) return (x < m) ? 0.0 : 1.0;
        return std_norm_cdf((x - m) / s);
    case PROB_DIST_EXPONENTIAL:
        if (x < 0.0) return 0.0; return 1.0 - exp(-m * x);
    case PROB_DIST_UNIFORM:
        if (x < m) return 0.0; if (x >= s) return 1.0;
        return (x - m) / (s - m);
    case PROB_DIST_GAMMA:
        if (x <= 0.0) return 0.0; return gamma_cdf_reg(m, x / s);
    case PROB_DIST_BETA:
        if (x <= 0.0) return 0.0; if (x >= 1.0) return 1.0;
        return beta_cdf_reg(m, s, x);
    case PROB_DIST_WEIBULL:
        if (x <= 0.0) return 0.0; return 1.0 - exp(-pow(x / m, s));
    case PROB_DIST_LOGNORMAL:
        if (x <= 0.0) return 0.0; return std_norm_cdf((log(x) - m) / s);
    case PROB_DIST_CAUCHY:
        return 0.5 + atan((x - m) / s) / M_PI;
    case PROB_DIST_CHISQ:
        if (x <= 0.0) return 0.0; return gamma_cdf_reg(m / 2.0, x / 2.0);
    case PROB_DIST_STUDENT_T: {
        double nu = m, tr = x * x / (nu + x * x);
        double ib = beta_cdf_reg(nu / 2.0, 0.5, 1.0 - tr);
        return (x >= 0.0) ? 1.0 - 0.5 * ib : 0.5 * ib;
    }
    case PROB_DIST_FISHER_F: {
        if (x <= 0.0) return 0.0;
        double d1 = m, d2 = s, xx = d1 * x / (d1 * x + d2);
        return beta_cdf_reg(d1 / 2.0, d2 / 2.0, xx);
    }
    default: return NAN;
    }
}

/*
 * prob_pdf_cont: density f(x) for 11 continuous families.
 * Knowledge: PDF is the derivative of CDF wherever it exists.
 * Each family has a closed-form density (Papoulis, Ch. 5).
 * O(1) per call.
 */
prob_real_t prob_pdf_cont(const prob_dist_cont_t *dist, prob_real_t x) {
    if (!dist) return NAN;
    double m = dist->param1, s = dist->param2;
    switch (dist->type) {
    case PROB_DIST_GAUSSIAN:
        if (s <= 0.0) return 0.0;
        return exp(-0.5 * sq((x - m) / s)) / (s * sqrt(2.0 * M_PI));
    case PROB_DIST_EXPONENTIAL:
        if (x < 0.0) return 0.0; return m * exp(-m * x);
    case PROB_DIST_UNIFORM:
        if (x < m || x > s) return 0.0; return 1.0 / (s - m);
    case PROB_DIST_GAMMA:
        if (x <= 0.0) return 0.0;
        return pow(x, m - 1.0) * exp(-x / s)
               / (pow(s, m) * tgamma(m));
    case PROB_DIST_BETA:
        if (x <= 0.0 || x >= 1.0) return 0.0;
        return pow(x, m - 1.0) * pow(1.0 - x, s - 1.0)
               / exp(lngamma(m) + lngamma(s) - lngamma(m + s));
    case PROB_DIST_WEIBULL:
        if (x <= 0.0) return 0.0;
        return (s / m) * pow(x / m, s - 1.0) * exp(-pow(x / m, s));
    case PROB_DIST_LOGNORMAL:
        if (x <= 0.0) return 0.0;
        return exp(-0.5 * sq((log(x) - m) / s))
               / (x * s * sqrt(2.0 * M_PI));
    case PROB_DIST_CAUCHY:
        return 1.0 / (M_PI * s * (1.0 + sq((x - m) / s)));
    case PROB_DIST_CHISQ:
        if (x <= 0.0) return 0.0;
        return pow(x, m / 2.0 - 1.0) * exp(-x / 2.0)
               / (pow(2.0, m / 2.0) * tgamma(m / 2.0));
    case PROB_DIST_STUDENT_T: {
        double nu = m;
        return exp(lngamma((nu + 1.0) / 2.0) - lngamma(nu / 2.0))
               / (sqrt(nu * M_PI))
               * pow(1.0 + x * x / nu, -(nu + 1.0) / 2.0);
    }
    case PROB_DIST_FISHER_F: {
        if (x <= 0.0) return 0.0;
        double d1 = m, d2 = s;
        return exp(lngamma((d1 + d2) / 2.0)
                   - lngamma(d1 / 2.0) - lngamma(d2 / 2.0))
               * pow(d1 / d2, d1 / 2.0) * pow(x, d1 / 2.0 - 1.0)
               * pow(1.0 + d1 * x / d2, -(d1 + d2) / 2.0);
    }
    default: return NAN;
    }
}

/*
 * prob_inverse_cdf: quantile function F^{-1}(p).
 * Knowledge: mapping [0,1] -> support.  Core to inverse-transform
 * sampling (Devroye, 1986).  Uses:
 *   - Abramowitz 26.2.23 rational approx for Gaussian
 *   - Analytic inverses for exponential, uniform, Cauchy, Weibull
 *   - Newton-Raphson (up to 60 iterations) for gamma, beta, t, F
 * O(1) analytic; O(iter) ~ O(1) for Newton.
 */
prob_real_t prob_inverse_cdf(const prob_dist_cont_t *dist, prob_t p) {
    if (!dist) return NAN;
    if (p <= 0.0 || p >= 1.0) return NAN;
    double m = dist->param1, s = dist->param2;
    switch (dist->type) {
    case PROB_DIST_GAUSSIAN: {
        double pp = (p < 0.5) ? p : 1.0 - p;
        double t = sqrt(-2.0 * log(pp));
        double num = 2.515517 + 0.802853 * t + 0.010328 * t * t;
        double den = 1.0 + 1.432788 * t + 0.189269 * t * t
                        + 0.001308 * t * t * t;
        double zp = t - num / den;
        if (p < 0.5) zp = -zp;
        return m + s * zp;
    }
    case PROB_DIST_EXPONENTIAL:
        return -log(1.0 - p) / m;
    case PROB_DIST_UNIFORM:
        return m + p * (s - m);
    case PROB_DIST_CAUCHY:
        return m + s * tan(M_PI * (p - 0.5));
    case PROB_DIST_WEIBULL:
        return m * pow(-log(1.0 - p), 1.0 / s);
    case PROB_DIST_LOGNORMAL: {
        prob_dist_cont_t nrm;
        prob_dist_cont_init(&nrm, PROB_DIST_GAUSSIAN, m, s);
        return exp(prob_inverse_cdf(&nrm, p));
    }
    case PROB_DIST_CHISQ: {
        prob_dist_cont_t gam;
        prob_dist_cont_init(&gam, PROB_DIST_GAMMA, m / 2.0, 2.0);
        return prob_inverse_cdf(&gam, p);
    }
    case PROB_DIST_GAMMA:
    case PROB_DIST_BETA:
    case PROB_DIST_STUDENT_T:
    case PROB_DIST_FISHER_F: {
        double q = isnan(dist->mean) ? 1.0 : dist->mean;
        int iter;
        for (iter = 0; iter < 60; iter++) {
            double F = prob_cdf_cont(dist, q) - p;
            double f = prob_pdf_cont(dist, q);
            if (f > 1e-300) {
                double step = F / f;
                q -= step;
                if (fabs(step) < 1e-10) break;
            } else {
                q += (p > prob_cdf_cont(dist, q))
                     ? 0.1 * s : -0.1 * s;
            }
        }
        return q;
    }
    default: return NAN;
    }
}

/* ================================================================
 * L2: Expectation operator
 * ================================================================ */

/*
 * prob_expect_disc: E[g(X)] = sum g(x_i) * p_i.
 * Knowledge: expectation is the Lebesgue integral w.r.t. probability
 * measure.  Linearity holds for any g.  O(n).
 */
prob_real_t prob_expect_disc(const prob_dist_disc_t *dist,
                             prob_real_t (*g)(prob_real_t)) {
    if (!dist || !dist->points || dist->n == 0 || !g) return NAN;
    double e = 0.0;
    size_t i;
    for (i = 0; i < dist->n; i++)
        e += g(dist->points[i].value) * dist->points[i].pmass;
    return e;
}

/*
 * prob_expect_cont: E[g(X)] under continuous distribution.
 * Gauss-Legendre 32-node for Gaussian (exact to polynomial degree 63).
 * Trapezoidal with 256 subdivisions for others.
 * O(1) for Gaussian; O(256) for others.
 */
prob_real_t prob_expect_cont(const prob_dist_cont_t *dist,
                             prob_real_t (*g)(prob_real_t)) {
    if (!dist || !g) return NAN;
    if (dist->type == PROB_DIST_GAUSSIAN) {
        static const double n32[32] = {
            -0.9972638618494816, -0.9856115115452684, -0.9647622555875064,
            -0.9349060759377397, -0.8963211557660521, -0.8493676137325700,
            -0.7944837959679424, -0.7321821187402897, -0.6630442669302152,
            -0.5877157572407623, -0.5068999089322294, -0.4213512761306353,
            -0.3318686022821277, -0.2392873622521371, -0.1444719615827965,
            -0.0483076656877383,  0.0483076656877383,  0.1444719615827965,
             0.2392873622521371,  0.3318686022821277,  0.4213512761306353,
             0.5068999089322294,  0.5877157572407623,  0.6630442669302152,
             0.7321821187402897,  0.7944837959679424,  0.8493676137325700,
             0.8963211557660521,  0.9349060759377397,  0.9647622555875064,
             0.9856115115452684,  0.9972638618494816
        };
        static const double w32[32] = {
            0.0070186100094701, 0.0162743947309057, 0.0253920653092621,
            0.0342738629130214, 0.0428358980222267, 0.0509980592623762,
            0.0586840934785355, 0.0658222227763618, 0.0723457941088485,
            0.0781938957870703, 0.0833119242269467, 0.0876520930044038,
            0.0911738786957639, 0.0938443990808046, 0.0956387200792749,
            0.0965400885147278, 0.0965400885147278, 0.0956387200792749,
            0.0938443990808046, 0.0911738786957639, 0.0876520930044038,
            0.0833119242269467, 0.0781938957870703, 0.0723457941088485,
            0.0658222227763618, 0.0586840934785355, 0.0509980592623762,
            0.0428358980222267, 0.0342738629130214, 0.0253920653092621,
            0.0162743947309057, 0.0070186100094701
        };
        double sc = M_SQRT2 * dist->param2, e = 0.0;
        int i;
        for (i = 0; i < 32; i++) {
            double x = dist->param1 + sc * n32[i];
            e += w32[i] * g(x) * prob_pdf_cont(dist, x);
        }
        return e * sc;
    }
    /* Fallback trapezoidal integration */
    double lo = prob_inverse_cdf(dist, 0.001);
    double hi = prob_inverse_cdf(dist, 0.999);
    int ns = 256;
    double h = (hi - lo) / ns, e = 0.0;
    int i;
    for (i = 0; i < ns; i++) {
        double x0 = lo + i * h, x1 = lo + (i + 1) * h;
        e += 0.5 * (g(x0) * prob_pdf_cont(dist, x0)
                  + g(x1) * prob_pdf_cont(dist, x1)) * h;
    }
    return e;
}

/*
 * prob_cond_expect: E[X | A] = E[X * I_A] / P(A).
 * Knowledge: conditional expectation is the Radon-Nikodym derivative
 * of the restricted measure.  O(n).
 */
prob_real_t prob_cond_expect(const prob_dist_disc_t *dist,
                             prob_t (*indicator)(prob_real_t),
                             prob_t prob_event) {
    if (!dist || !indicator || prob_event <= 0.0) return NAN;
    double s = 0.0;
    size_t i;
    for (i = 0; i < dist->n; i++)
        s += dist->points[i].value
           * indicator(dist->points[i].value)
           * dist->points[i].pmass;
    return s / prob_event;
}

/*
 * prob_bayes: P(H|E) = P(E|H) * P(H) / P(E).
 * Knowledge: Bayes' theorem is the foundation of all Bayesian
 * inference and machine learning.  O(1).
 */
prob_t prob_bayes(prob_t prior_H, prob_t likelihood_E_given_H,
                  prob_t prob_E) {
    if (prob_E <= 0.0) return NAN;
    return prior_H * likelihood_E_given_H / prob_E;
}

/* ================================================================
 * L6: Distribution Construction
 * ================================================================ */

/*
 * prob_dist_disc_create: build discrete distribution from (value, pmass)
 * pairs.  Normalises masses to sum 1, caches mean and variance.
 * Returns NULL on invalid input.  O(n), O(n) space.
 */
prob_dist_disc_t *prob_dist_disc_create(const prob_dpoint_t *points,
                                        size_t n) {
    if (!points || n == 0) return NULL;
    prob_dist_disc_t *d = (prob_dist_disc_t *)malloc(
        sizeof(prob_dist_disc_t));
    if (!d) return NULL;
    d->n = n;
    d->points = (prob_dpoint_t *)malloc(n * sizeof(prob_dpoint_t));
    if (!d->points) { free(d); return NULL; }
    memcpy(d->points, points, n * sizeof(prob_dpoint_t));
    double tot = 0.0;
    size_t i;
    for (i = 0; i < n; i++) tot += d->points[i].pmass;
    if (tot > 0.0)
        for (i = 0; i < n; i++) d->points[i].pmass /= tot;
    d->mean = 0.0;
    for (i = 0; i < n; i++)
        d->mean += d->points[i].value * d->points[i].pmass;
    d->var = 0.0;
    for (i = 0; i < n; i++) {
        double delta = d->points[i].value - d->mean;
        d->var += delta * delta * d->points[i].pmass;
    }
    return d;
}

void prob_dist_disc_free(prob_dist_disc_t *dist) {
    if (dist) { free(dist->points); free(dist); }
}

/*
 * prob_dist_cont_init: initialise continuous distribution with
 * analytical moments (mean, var, skew, kurt) for the given family.
 * Knowledge: moments derivable from characteristic function or MGF
 * (Papoulis, Ch. 5).  O(1).
 */
void prob_dist_cont_init(prob_dist_cont_t *dist, prob_dist_type_t type,
                         prob_real_t p1, prob_real_t p2) {
    if (!dist) return;
    dist->type = type; dist->param1 = p1; dist->param2 = p2;
    switch (type) {
    case PROB_DIST_GAUSSIAN:
        dist->mean = p1; dist->var = p2 * p2;
        dist->skew = 0.0; dist->kurt = 0.0; break;
    case PROB_DIST_EXPONENTIAL:
        dist->mean = 1.0 / p1; dist->var = 1.0 / (p1 * p1);
        dist->skew = 2.0; dist->kurt = 6.0; break;
    case PROB_DIST_UNIFORM:
        dist->mean = (p1 + p2) / 2.0;
        dist->var  = (p2 - p1) * (p2 - p1) / 12.0;
        dist->skew = 0.0; dist->kurt = -1.2; break;
    case PROB_DIST_GAMMA:
        dist->mean = p1 * p2; dist->var = p1 * p2 * p2;
        dist->skew = 2.0 / sqrt(p1); dist->kurt = 6.0 / p1; break;
    case PROB_DIST_BETA: {
        double a = p1, b = p2, ab = a + b;
        dist->mean = a / ab;
        dist->var  = a * b / (ab * ab * (ab + 1.0));
        dist->skew = 2.0 * (b - a) * sqrt(ab + 1.0)
                     / ((ab + 2.0) * sqrt(a * b));
        dist->kurt = 6.0 * (cu(a - b) * (ab + 1.0)
                     - a * b * (ab + 2.0))
                     / (a * b * (ab + 2.0) * (ab + 3.0));
        break;
    }
    case PROB_DIST_WEIBULL: {
        double k = p2;
        double g1 = tgamma(1.0 + 1.0 / k);
        double g2 = tgamma(1.0 + 2.0 / k);
        double g3 = tgamma(1.0 + 3.0 / k);
        double g4 = tgamma(1.0 + 4.0 / k);
        dist->mean = p1 * g1;
        dist->var  = p1 * p1 * (g2 - g1 * g1);
        double s3 = pow(dist->var, 1.5);
        dist->skew = (g3 - 3.0 * g2 * g1 + 2.0 * cu(g1))
                     * p1 * p1 * p1 / s3;
        dist->kurt = (g4 - 4.0 * g3 * g1 + 6.0 * g2 * g1 * g1
                     - 3.0 * sq(sq(g1)))
                     * pow(p1, 4.0) / sq(dist->var) - 3.0;
        break;
    }
    case PROB_DIST_LOGNORMAL:
        dist->mean = exp(p1 + p2 * p2 / 2.0);
        dist->var  = (exp(p2 * p2) - 1.0)
                     * exp(2.0 * p1 + p2 * p2);
        dist->skew = (exp(p2 * p2) + 2.0)
                     * sqrt(exp(p2 * p2) - 1.0);
        dist->kurt = exp(4.0 * p2 * p2)
                     + 2.0 * exp(3.0 * p2 * p2)
                     + 3.0 * exp(2.0 * p2 * p2) - 6.0;
        break;
    case PROB_DIST_CAUCHY:
        dist->mean = NAN; dist->var = NAN;
        dist->skew = NAN; dist->kurt = NAN; break;
    case PROB_DIST_CHISQ:
        dist->mean = p1; dist->var = 2.0 * p1;
        dist->skew = sqrt(8.0 / p1);
        dist->kurt = 12.0 / p1; break;
    case PROB_DIST_STUDENT_T: {
        double nu = p1;
        dist->mean = (nu > 1.0) ? 0.0 : NAN;
        dist->var  = (nu > 2.0) ? nu / (nu - 2.0) : NAN;
        dist->skew = (nu > 3.0) ? 0.0 : NAN;
        dist->kurt = (nu > 4.0) ? 6.0 / (nu - 4.0) : NAN;
        break;
    }
    case PROB_DIST_FISHER_F: {
        double d1 = p1, d2 = p2;
        dist->mean = (d2 > 2.0) ? d2 / (d2 - 2.0) : NAN;
        if (d2 > 4.0)
            dist->var = 2.0 * d2 * d2 * (d1 + d2 - 2.0)
                / (d1 * (d2 - 2.0) * (d2 - 2.0) * (d2 - 4.0));
        else dist->var = NAN;
        dist->skew = NAN; dist->kurt = NAN; break;
    }
    default:
        dist->mean = 0.0; dist->var = 0.0;
        dist->skew = 0.0; dist->kurt = 0.0;
    }
}
