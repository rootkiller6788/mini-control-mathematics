/**
 * @file ex01_gaussian_mle.c
 * @brief Example: Maximum likelihood estimation for Gaussian distribution
 *        with confidence interval via bootstrap.
 *
 * Demonstrates: prob_mle_gaussian(), prob_dist_cont_init(),
 * prob_bootstrap_ci(), prob_ks_test().
 *
 * Scenario: Quality control — estimate mean and variance of a
 * manufacturing process from measurements, test normality,
 * and compute a 95% bootstrap CI for the mean.
 *
 * Knowledge: L6 engineering problems (parameter estimation,
 * goodness-of-fit, uncertainty quantification).
 */
#include "probability_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(void) {
    printf("=== Example 1: Gaussian MLE and Bootstrap CI ===\n\n");

    /* Simulate 50 measurements from N(10.0, 0.5^2) */
    uint64_t seed = 42;
    double data[50];
    int i;
    printf("Simulated measurements (N(10, 0.5^2)):\n  ");
    for (i = 0; i < 50; i++) {
        data[i] = prob_sample_gaussian(10.0, 0.5, &seed);
        if (i < 8) printf("%.3f ", data[i]);
    }
    printf("...\n\n");

    /* MLE estimation */
    double mu_hat, sigma2_hat;
    prob_mle_gaussian(data, 50, &mu_hat, &sigma2_hat);
    double sigma_hat = sqrt(sigma2_hat);

    printf("True parameters:     mu=10.000, sigma=0.500\n");
    printf("MLE estimates:       mu=%.3f, sigma=%.3f\n", mu_hat, sigma_hat);
    printf("Bias in sigma (MLE): sigma_mle uses n (not n-1), biased low\n\n");

    /* KS test for normality */
    prob_dist_cont_t gauss_fitted;
    prob_dist_cont_init(&gauss_fitted, PROB_DIST_GAUSSIAN, mu_hat, sigma_hat);
    double ks_stat = prob_ks_test(data, 50, &gauss_fitted);

    printf("Kolmogorov-Smirnov normality test:\n");
    printf("  KS statistic D_50 = %.4f\n", ks_stat);
    /* Critical value for alpha=0.05, n=50 is approx 0.188 */
    printf("  Critical value (a=0.05): 0.188\n");
    printf("  Result: %s normality at 5%% level\n\n",
           ks_stat < 0.188 ? "Cannot reject" : "Reject");

    /* Bootstrap 95% CI for the mean */
    double ci_lo, ci_hi;
    prob_bootstrap_ci(data, 50, 5000, 0.95, &ci_lo, &ci_hi);
    printf("95%% Bootstrap CI for mean (B=5000):\n");
    printf("  [%.3f, %.3f]\n", ci_lo, ci_hi);
    printf("  Width: %.3f\n", ci_hi - ci_lo);
    printf("  Contains true mean (10.0): %s\n",
           (ci_lo <= 10.0 && ci_hi >= 10.0) ? "YES" : "NO");

    /* Information-theoretic summary */
    printf("\nFitted Gaussian properties:\n");
    printf("  Differential entropy: %.4f nats\n",
           prob_differential_entropy(&gauss_fitted));

    return 0;
}
