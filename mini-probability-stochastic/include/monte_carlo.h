/**
 * @file monte_carlo.h
 * @brief Monte Carlo methods: integration, sampling, variance reduction,
 *        Markov chain Monte Carlo (MCMC), and quasi-Monte Carlo.
 *
 * Knowledge coverage:
 *   L1 - Definitions:  Monte Carlo estimator, importance weight, proposal
 *   L2 - Core Concepts: law of large numbers for MC, convergence rate O(1/sqrt(N))
 *   L5 - Engineering Methods: importance sampling, rejection sampling,
 *                             Metropolis-Hastings, Gibbs sampling,
 *                             antithetic variates, control variates,
 *                             quasi-MC (Halton, Sobol)
 *   L7 - Applications: Bayesian inference, rare-event simulation,
 *                      high-dimensional integration
 *   L8 - Advanced Methods: MCMC convergence diagnostics, Hamiltonian MC
 *
 * Reference texts:
 *   - Robert & Casella, "Monte Carlo Statistical Methods", 2nd ed. (2004)
 *   - Liu, "Monte Carlo Strategies in Scientific Computing" (2001)
 *   - Glasserman, "Monte Carlo Methods in Financial Engineering" (2003)
 *   - Owen, "Monte Carlo Theory, Methods and Examples" (2013)
 *
 * Course mapping:
 *   MIT 6.437 (Inference and Information)
 *   Stanford STATS 362 (Monte Carlo)
 *   Berkeley STAT 240A (Statistical Computing)
 *   Cambridge Part III (Monte Carlo Inference)
 */

#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include "probability_core.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L1: Basic Monte Carlo integration
 * ================================================================ */

/**
 * @brief Monte Carlo integration: estimate integral_a^b f(x) dx.
 *
 * I_hat = (b - a) * (1/N) * sum_{i=1}^N f(X_i), X_i ~ Uniform(a,b).
 *
 * Standard error: SE = (b - a) * sigma_f / sqrt(N).
 * Convergence rate: O(1/sqrt(N)) ? independent of dimension.
 *
 * Reference: Robert & Casella, Ch. 3.
 * Complexity: O(N).
 *
 * @param f      Integrand function.
 * @param a      Lower bound.
 * @param b      Upper bound.
 * @param N      Number of samples.
 * @param seed   RNG seed.
 * @param se     Output: standard error estimate.
 * @return Monte Carlo estimate of the integral.
 */
prob_real_t mc_integrate_1d(prob_real_t (*f)(prob_real_t),
                            prob_real_t a, prob_real_t b,
                            size_t N, uint64_t seed,
                            prob_real_t *se);

/**
 * @brief Monte Carlo integration over a hyper-rectangle in d dimensions.
 *
 * @param f      d-dimensional function f(x_0, ..., x_{d-1}).
 * @param a      Lower bounds array (d elements).
 * @param b      Upper bounds array (d elements).
 * @param d      Dimension.
 * @param N      Number of samples.
 * @param seed   RNG seed.
 * @param se     Standard error estimate.
 * @return Estimate of the integral.
 */
prob_real_t mc_integrate_nd(prob_real_t (*f)(const prob_real_t *, size_t),
                            const prob_real_t *a, const prob_real_t *b,
                            size_t d, size_t N, uint64_t seed,
                            prob_real_t *se);

/* ================================================================
 * L5: Importance sampling
 * ================================================================ */

/**
 * @brief Importance sampling integration:
 *        E_f[h(X)] = E_g[h(X) * f(X) / g(X)].
 *
 * Samples from proposal distribution g, reweights by likelihood ratio.
 * Optimal g (zero-variance): g*(x) propto |h(x)| * f(x).
 *
 * Reference: Robert & Casella, Ch. 3.3.
 * Complexity: O(N).
 *
 * @param h           Target function.
 * @param pdf_target  Target density f(x).
 * @param pdf_proposal Proposal density g(x).
 * @param sample_proposal Function to sample from g.
 * @param N           Number of samples.
 * @param seed        RNG seed.
 * @param se          Standard error estimate.
 * @return Importance sampling estimate.
 */
prob_real_t mc_importance_sampling(prob_real_t (*h)(prob_real_t),
                                   prob_real_t (*pdf_target)(prob_real_t),
                                   prob_real_t (*pdf_proposal)(prob_real_t),
                                   prob_real_t (*sample_proposal)(uint64_t*),
                                   size_t N, uint64_t seed,
                                   prob_real_t *se);

/* ================================================================
 * L5: Rejection sampling
 * ================================================================ */

/**
 * @brief Rejection sampling: sample from target f using envelope M * g.
 *
 * Requires f(x) <= M * g(x) for all x.
 * Acceptance probability = 1 / M.
 *
 * Reference: Robert & Casella, Ch. 2.3.
 * Complexity: O(n * M) expected evaluations.
 *
 * @param pdf_target  Target density f(x) (unnormalised OK).
 * @param pdf_env     Envelope density g(x).
 * @param sample_env  Sampler from envelope.
 * @param M           Envelope constant (must satisfy f(x) <= M*g(x)).
 * @param n           Number of samples desired.
 * @param samples     Output array (n, pre-allocated).
 * @param seed        RNG seed pointer (updated).
 * @return Actual number accepted (may be < n if max attempts exceeded).
 */
size_t mc_rejection_sampling(prob_real_t (*pdf_target)(prob_real_t),
                             prob_real_t (*pdf_env)(prob_real_t),
                             prob_real_t (*sample_env)(uint64_t*),
                             prob_real_t M, size_t n,
                             prob_real_t *samples, uint64_t *seed);

/* ================================================================
 * L5: Metropolis-Hastings algorithm
 * ================================================================ */

/**
 * @brief Metropolis-Hastings MCMC sampler.
 *
 * Generates a Markov chain with stationary distribution proportional
 * to target_pdf. Uses symmetric random-walk proposal by default.
 *
 * Reference: Hastings, "Monte Carlo Sampling Methods Using Markov
 *            Chains and Their Applications", Biometrika (1970).
 * Complexity: O(n_iterations).
 *
 * @param target_log_pdf  Log of (unnormalised) target density.
 * @param n_iter          Number of iterations (including burn-in).
 * @param burn_in         Number of initial samples to discard.
 * @param x0              Starting value.
 * @param proposal_sd     Standard deviation of Gaussian proposal.
 * @param samples         Output array (n_iter - burn_in samples).
 * @param seed            RNG seed pointer (updated).
 * @return Acceptance rate in [0, 1].
 */
prob_real_t mc_metropolis_hastings(prob_real_t (*target_log_pdf)(prob_real_t),
                                   size_t n_iter, size_t burn_in,
                                   prob_real_t x0, prob_real_t proposal_sd,
                                   prob_real_t *samples, uint64_t *seed);

/**
 * @brief Multivariate Metropolis-Hastings with Gaussian proposal.
 *
 * @param target_log_pdf  Log target density of d-dimensional vector.
 * @param d               Dimension.
 * @param n_iter          Total iterations.
 * @param burn_in         Burn-in period.
 * @param x0              Starting point (d elements).
 * @param proposal_cov    Proposal covariance diagonal (d elements).
 * @param samples         Output flat array (d * (n_iter-burn_in)),
 *                        column-major: sample j, dimension i at i+j*d.
 * @param seed            RNG seed pointer.
 * @return Acceptance rate.
 */
prob_real_t mc_metropolis_hastings_nd(
    prob_real_t (*target_log_pdf)(const prob_real_t *, size_t),
    size_t d, size_t n_iter, size_t burn_in,
    const prob_real_t *x0, const prob_real_t *proposal_cov,
    prob_real_t *samples, uint64_t *seed);

/* ================================================================
 * L5: Gibbs sampling
 * ================================================================ */

/**
 * @brief Gibbs sampler for a two-component distribution.
 *
 * Iteratively samples from conditional distributions:
 *   X_1^{(t+1)} ~ p(x_1 | x_2^{(t)})
 *   X_2^{(t+1)} ~ p(x_2 | x_1^{(t+1)})
 *
 * Reference: Geman & Geman, "Stochastic Relaxation, Gibbs
 *            Distributions, and Bayesian Restoration of Images"
 *            IEEE TPAMI (1984).
 * Complexity: O(n_iter) per chain.
 *
 * @param sample_x1_given_x2  Sampler for p(x1 | x2).
 * @param sample_x2_given_x1  Sampler for p(x2 | x1).
 * @param n_iter              Iterations.
 * @param burn_in             Burn-in.
 * @param x1_init, x2_init    Initial values.
 * @param samples_x1          Output x1 samples.
 * @param samples_x2          Output x2 samples.
 * @param seed                RNG seed pointer.
 */
void mc_gibbs_sampler_2d(prob_real_t (*sample_x1_given_x2)(prob_real_t, uint64_t*),
                         prob_real_t (*sample_x2_given_x1)(prob_real_t, uint64_t*),
                         size_t n_iter, size_t burn_in,
                         prob_real_t x1_init, prob_real_t x2_init,
                         prob_real_t *samples_x1, prob_real_t *samples_x2,
                         uint64_t *seed);

/* ================================================================
 * L5: Variance reduction techniques
 * ================================================================ */

/**
 * @brief Antithetic variates for MC integration.
 *
 * Uses pairs (X_i, X_i') where X_i' is negatively correlated with X_i.
 * For uniform: X_i' = 1 - X_i.
 * Variance reduction: Var((X+X')/2) <= Var(X)/2 when rho < 0.
 *
 * Reference: Glasserman, Ch. 4.2.
 * Complexity: O(N).
 */
prob_real_t mc_antithetic_variates(prob_real_t (*f)(prob_real_t),
                                   prob_real_t a, prob_real_t b,
                                   size_t N, uint64_t seed,
                                   prob_real_t *se);

/**
 * @brief Control variates for MC integration.
 *
 * Estimates E[f(X)] using E[f(X) - c*(g(X) - E[g])] where E[g] is known.
 * Optimal c* = Cov(f, g) / Var(g).
 *
 * Reference: Glasserman, Ch. 4.1.
 * Complexity: O(N) with linear regression.
 *
 * @param f      Target function.
 * @param g      Control variate function.
 * @param Eg     Known expectation E[g(X)].
 * @param a, b   Uniform bounds.
 * @param N      Sample size.
 * @param seed   RNG seed.
 * @param se     Standard error.
 * @param c_opt  Output: optimal control coefficient.
 * @return Controlled estimate of E[f(X)].
 */
prob_real_t mc_control_variates(prob_real_t (*f)(prob_real_t),
                                prob_real_t (*g)(prob_real_t),
                                prob_real_t Eg,
                                prob_real_t a, prob_real_t b,
                                size_t N, uint64_t seed,
                                prob_real_t *se, prob_real_t *c_opt);

/* ================================================================
 * L5: Quasi-Monte Carlo (low-discrepancy sequences)
 * ================================================================ */

/**
 * @brief Generate Halton sequence points.
 *
 * Low-discrepancy sequence with asymptotic discrepancy O((log N)^d / N).
 * Better than random for moderate dimensions (d <= 20).
 *
 * Reference: Halton, "On the Efficiency of Certain Quasi-Random
 *            Sequences", Numerische Mathematik (1960).
 * Complexity: O(d * N).
 *
 * @param d       Dimension.
 * @param N       Number of points.
 * @param points  Output array (N * d), row-major.
 * @param start   Starting index in Halton sequence (>= 1).
 */
void qmc_halton_sequence(size_t d, size_t N, prob_real_t *points, size_t start);

/**
 * @brief Generate Sobol sequence points (simplified, up to d=8).
 *
 * Often better uniformity than Halton for moderate dimensions.
 * Reference: Sobol, "On the Distribution of Points in a Cube",
 *            USSR Comp. Math. Math. Phys. (1967).
 *
 * @param d       Dimension (<= 8).
 * @param N       Number of points (must be power of 2 for full uniformity).
 * @param points  Output array (N * d), row-major.
 */
void qmc_sobol_sequence(size_t d, size_t N, prob_real_t *points);

/**
 * @brief Quasi-Monte Carlo integration using given low-discrepancy points.
 *
 * @return Estimate of the integral.
 */
prob_real_t qmc_integrate(prob_real_t (*f)(const prob_real_t *, size_t),
                          size_t d, const prob_real_t *points, size_t N);

/* ================================================================
 * L8: MCMC convergence diagnostics
 * ================================================================ */

/**
 * @brief Gelman-Rubin R-hat diagnostic for MCMC convergence.
 *
 * Compares multiple chains; R_hat ~ 1 indicates convergence.
 * Reference: Gelman & Rubin, "Inference from Iterative Simulation
 *            Using Multiple Sequences", Statistical Science (1992).
 *
 * @param chains  Array of m chains, each of length n. Flat: chain c,
 *                iteration i at chains[c*n + i].
 * @param m       Number of chains.
 * @param n       Length of each chain.
 * @return R-hat statistic (>= 1). R_hat < 1.1 suggests convergence.
 */
prob_real_t mcmc_gelman_rubin(const prob_real_t *chains,
                              size_t m, size_t n);

/**
 * @brief Effective sample size (ESS) for autocorrelated MCMC output.
 *
 * ESS = N / (1 + 2 * sum_{k=1}^\infty rho_k).
 * Estimates using initial monotone sequence criterion (Geyer, 1992).
 *
 * Reference: Kass, Carlin, Gelman, Neal, "Markov Chain Monte Carlo
 *            in Practice", 1998.
 *
 * @param chain    Single-chain samples.
 * @param n        Chain length.
 * @param max_lag  Maximum autocorrelation lag to consider.
 * @return Effective sample size.
 */
prob_real_t mcmc_effective_sample_size(const prob_real_t *chain,
                                       size_t n, size_t max_lag);

#ifdef __cplusplus
}
#endif

#endif /* MONTE_CARLO_H */
