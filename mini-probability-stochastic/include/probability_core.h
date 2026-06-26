/**
 * @file probability_core.h
 * @brief Core probability definitions: random variables, distributions,
 *        expectation, variance, covariance, conditional probability.
 *
 * Knowledge coverage:
 *   L1 - Definitions:  random variable, distribution, expectation, variance
 *   L2 - Core Concepts: law of large numbers, central limit theorem,
 *                       Bayes theorem, conditional probability
 *   L3 - Engineering Quantities: mean, variance, skewness, kurtosis,
 *                                covariance matrix
 *   L4 - Conservation Laws: total probability theorem, Bayes rule identity
 *
 * Reference texts:
 *   - Papoulis & Pillai, "Probability, Random Variables and Stochastic
 *     Processes", 4th ed. (2002)
 *   - Ross, "A First Course in Probability", 9th ed. (2014)
 *   - Feller, "An Introduction to Probability Theory and Its Applications"
 *
 * Course mapping:
 *   MIT 6.041 (Probabilistic Systems Analysis)
 *   Stanford EE 178 (Probabilistic Systems Analysis)
 *   Berkeley EECS 126 (Probability and Random Processes)
 *   Cambridge Part II (Probability and Measure)
 */

#ifndef PROBABILITY_CORE_H
#define PROBABILITY_CORE_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L1: Core Definitions
 * ================================================================ */

typedef double prob_real_t;
typedef double prob_t;

typedef struct {
    prob_real_t value;
    prob_t      pmass;
} prob_dpoint_t;

typedef struct {
    prob_dpoint_t *points;
    size_t         n;
    prob_real_t    mean;
    prob_real_t    var;
} prob_dist_disc_t;

typedef enum {
    PROB_DIST_GAUSSIAN    = 0,
    PROB_DIST_EXPONENTIAL = 1,
    PROB_DIST_UNIFORM     = 2,
    PROB_DIST_GAMMA       = 3,
    PROB_DIST_BETA        = 4,
    PROB_DIST_WEIBULL     = 5,
    PROB_DIST_LOGNORMAL   = 6,
    PROB_DIST_CAUCHY      = 7,
    PROB_DIST_CHISQ       = 8,
    PROB_DIST_STUDENT_T   = 9,
    PROB_DIST_FISHER_F    = 10,
    PROB_DIST_CUSTOM      = 99
} prob_dist_type_t;

typedef struct {
    prob_dist_type_t type;
    prob_real_t      param1;
    prob_real_t      param2;
    prob_real_t      mean;
    prob_real_t      var;
    prob_real_t      skew;
    prob_real_t      kurt;
} prob_dist_cont_t;

typedef struct {
    prob_real_t *samples;
    size_t       n;
    size_t       cap;
    prob_real_t  sum;
    prob_real_t  mean;
    prob_real_t  var;
    prob_real_t  skew;
    prob_real_t  kurt;
    prob_real_t  m2;
    prob_real_t  m3;
    prob_real_t  m4;
} prob_empirical_t;

/* L2: Distribution queries */
prob_t       prob_pmf         (const prob_dist_disc_t *dist, prob_real_t x);
prob_t       prob_cdf_disc    (const prob_dist_disc_t *dist, prob_real_t x);
prob_t       prob_cdf_cont    (const prob_dist_cont_t *dist, prob_real_t x);
prob_real_t  prob_pdf_cont    (const prob_dist_cont_t *dist, prob_real_t x);
prob_real_t  prob_inverse_cdf (const prob_dist_cont_t *dist, prob_t p);

prob_real_t  prob_expect_disc (const prob_dist_disc_t *dist,
                               prob_real_t (*g)(prob_real_t));
prob_real_t  prob_expect_cont (const prob_dist_cont_t *dist,
                               prob_real_t (*g)(prob_real_t));
prob_real_t  prob_cond_expect (const prob_dist_disc_t *dist,
                               prob_t (*indicator)(prob_real_t),
                               prob_t prob_event);
prob_t       prob_bayes       (prob_t prior_H,
                               prob_t likelihood_E_given_H,
                               prob_t prob_E);

/* L3: Engineering Quantities */
prob_real_t  prob_mean_disc      (const prob_dist_disc_t *dist);
prob_real_t  prob_var_disc       (const prob_dist_disc_t *dist);
prob_real_t  prob_std_disc       (const prob_dist_disc_t *dist);
prob_real_t  prob_skewness_disc  (const prob_dist_disc_t *dist);
prob_real_t  prob_kurtosis_disc  (const prob_dist_disc_t *dist);
prob_real_t  prob_covariance     (const prob_real_t *x, const prob_real_t *y,
                                  size_t n);
prob_real_t  prob_correlation    (const prob_real_t *x, const prob_real_t *y,
                                  size_t n);
void         prob_covariance_matrix(const prob_real_t *data,
                                    size_t k, size_t n, prob_real_t *cov);

/* L4: Conservation Theorems */
prob_real_t  prob_total_probability_check (const prob_t *partition, size_t n);
prob_real_t  prob_total_expectation_check (const prob_dist_disc_t *dist,
                                          const prob_t *event_probs,
                                          prob_t (**indicators)(prob_real_t),
                                          size_t m);
prob_real_t  prob_lln_check     (const prob_real_t *samples, size_t n,
                                 prob_real_t true_mean);
prob_real_t  prob_clt_ks_test   (const prob_real_t *standardised_means,
                                 size_t batch_size);

/* L6: Construction and Fitting */
prob_dist_disc_t *prob_dist_disc_create (const prob_dpoint_t *points, size_t n);
void              prob_dist_disc_free   (prob_dist_disc_t *dist);
void              prob_dist_cont_init   (prob_dist_cont_t *dist,
                                        prob_dist_type_t type,
                                        prob_real_t param1, prob_real_t param2);
void              prob_mle_gaussian     (const prob_real_t *data, size_t n,
                                        prob_real_t *mu_hat,
                                        prob_real_t *sigma2_hat);
prob_real_t       prob_mle_exponential  (const prob_real_t *data, size_t n);
prob_real_t       prob_mle_poisson      (const int32_t *data, size_t n);
prob_real_t       prob_ks_test          (const prob_real_t *data, size_t n,
                                        const prob_dist_cont_t *dist);

prob_empirical_t *prob_empirical_create  (size_t capacity);
void              prob_empirical_insert  (prob_empirical_t *emp, prob_real_t x);
void              prob_empirical_free    (prob_empirical_t *emp);
prob_real_t       prob_empirical_quantile(const prob_empirical_t *emp, prob_t p);
void              prob_bootstrap_ci      (const prob_real_t *data, size_t n,
                                         size_t B, prob_t confidence,
                                         prob_real_t *lower, prob_real_t *upper);

/* L5: Information Theory */
prob_real_t  prob_entropy_disc          (const prob_dist_disc_t *dist);
prob_real_t  prob_differential_entropy  (const prob_dist_cont_t *dist);
prob_real_t  prob_kl_divergence         (const prob_dist_disc_t *p,
                                         const prob_dist_disc_t *q);
prob_real_t  prob_js_divergence         (const prob_dist_disc_t *p,
                                         const prob_dist_disc_t *q);
prob_real_t  prob_mutual_information    (const prob_real_t *joint,
                                         const prob_real_t *px, size_t nx,
                                         const prob_real_t *py, size_t ny);

/* ================================================================
 * L5: PRNG and Distribution Sampling (src/probability_core.c)
 * ================================================================ */

double        prob_lcg64_uniform       (uint64_t *state);
double        prob_xorshift64_uniform  (uint64_t *state);
double        prob_box_muller          (uint64_t *state);
double        prob_marsaglia_polar     (uint64_t *state);

double        prob_sample_uniform      (double a, double b, uint64_t *state);
double        prob_sample_exponential  (double lambda, uint64_t *state);
double        prob_sample_gaussian     (double mu, double sigma, uint64_t *state);
int           prob_sample_bernoulli    (double p, uint64_t *state);
int           prob_sample_binomial     (int n, double p, uint64_t *state);
int           prob_sample_poisson      (double lambda, uint64_t *state);
double        prob_sample_gamma        (double alpha, double beta, uint64_t *state);
double        prob_sample_beta         (double a, double b, uint64_t *state);
size_t        prob_sample_discrete     (const double *probs, size_t n,
                                        uint64_t *state);

/* ================================================================
 * L1/L3: Special Mathematical Functions (src/probability_core.c)
 * ================================================================ */

double        prob_erf_approx          (double x);
double        prob_std_normal_cdf      (double x);
double        prob_lgamma              (double x);
double        prob_gamma_cdf_reg       (double a, double x);
double        prob_beta_cdf_reg        (double a, double b, double x);

#ifdef __cplusplus
}
#endif

#endif /* PROBABILITY_CORE_H */
