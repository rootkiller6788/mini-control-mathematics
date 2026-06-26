# Knowledge Graph — mini-probability-stochastic

## L1: Definitions (Complete)

| # | Concept | C Implementation | Location |
|---|---------|-----------------|----------|
| 1 | Random variable (prob_real_t) | `typedef double prob_real_t` | `include/probability_core.h` |
| 2 | Probability value (prob_t) | `typedef double prob_t` | `include/probability_core.h` |
| 3 | Discrete probability mass point | `typedef struct prob_dpoint_t` | `include/probability_core.h` |
| 4 | Discrete distribution | `typedef struct prob_dist_disc_t` | `include/probability_core.h` |
| 5 | Continuous distribution family enum (11 families) | `typedef enum prob_dist_type_t` | `include/probability_core.h` |
| 6 | Continuous distribution descriptor | `typedef struct prob_dist_cont_t` | `include/probability_core.h` |
| 7 | Empirical distribution | `typedef struct prob_empirical_t` | `include/probability_core.h` |
| 8 | Markov chain (DTMC) | `typedef struct markov_chain_dtmc_t` | `include/stochastic_process.h` |
| 9 | Continuous-time Markov chain (CTMC) | `typedef struct markov_chain_ctmc_t` | `include/stochastic_process.h` |
| 10 | Poisson process | `poisson_process_simulate()` | `include/stochastic_process.h` |
| 11 | Brownian motion / Wiener process | `brownian_motion_simulate()` | `include/stochastic_process.h` |
| 12 | Geometric Brownian motion | `geometric_brownian_motion()` | `include/stochastic_process.h` |
| 13 | Ornstein-Uhlenbeck process | `ornstein_uhlenbeck_simulate()` | `include/stochastic_process.h` |
| 14 | Kalman state-space model | `typedef struct kalman_model_t` | `include/estimation_filtering.h` |
| 15 | Nonlinear state-space model | `typedef struct nlinear_model_t` | `include/estimation_filtering.h` |

## L2: Core Concepts (Complete)

| # | Concept | C Implementation | Location |
|---|---------|-----------------|----------|
| 1 | Probability mass function PMF | `prob_pmf()` | `src/probability_distributions.c` |
| 2 | Cumulative distribution function CDF | `prob_cdf_disc()`, `prob_cdf_cont()` | `src/probability_distributions.c` |
| 3 | Probability density function PDF | `prob_pdf_cont()` | `src/probability_distributions.c` |
| 4 | Inverse CDF / Quantile function | `prob_inverse_cdf()` | `src/probability_distributions.c` |
| 5 | Expectation operator E[g(X)] | `prob_expect_disc()`, `prob_expect_cont()` | `src/probability_distributions.c` |
| 6 | Conditional expectation E[X|A] | `prob_cond_expect()` | `src/probability_distributions.c` |
| 7 | Bayes theorem | `prob_bayes()` | `src/probability_distributions.c` |
| 8 | Markov property (DTMC/CTMC) | `markov_dtmc_*`, `markov_ctmc_*` | `src/stochastic_process.c` |
| 9 | Stationarity (stationary distribution) | `markov_dtmc_stationary()` | `src/stochastic_process.c` |
| 10 | Ergodicity (simulation convergence) | `markov_dtmc_simulate()` | `src/stochastic_process.c` |
| 11 | Detailed balance | `markov_dtmc_detailed_balance()` | `src/stochastic_process.c` |
| 12 | Memoryless property (exponential) | `prob_sample_exponential()` | `src/probability_core.c` |
| 13 | Direct sampling (inverse transform) | `prob_sample_*` functions | `src/probability_core.c` |

## L3: Engineering Quantities (Complete)

| # | Concept | C Implementation | Location |
|---|---------|-----------------|----------|
| 1 | Mean (first moment) | `prob_mean_disc()` | `src/probability_stats.c` |
| 2 | Variance (second central moment) | `prob_var_disc()` | `src/probability_stats.c` |
| 3 | Standard deviation | `prob_std_disc()` | `src/probability_stats.c` |
| 4 | Skewness (standardised third moment) | `prob_skewness_disc()` | `src/probability_stats.c` |
| 5 | Excess kurtosis (standardised fourth moment - 3) | `prob_kurtosis_disc()` | `src/probability_stats.c` |
| 6 | Covariance Cov[X,Y] | `prob_covariance()` | `src/probability_stats.c` |
| 7 | Pearson correlation rho | `prob_correlation()` | `src/probability_stats.c` |
| 8 | Covariance matrix (multivariate) | `prob_covariance_matrix()` | `src/probability_stats.c` |
| 9 | Sample autocorrelation function | `time_series_autocorr()` | `src/stochastic_process.c` |
| 10 | Periodogram PSD | `periodogram_psd()` | `src/stochastic_process.c` |
| 11 | Analytical moments for 11 families | `prob_dist_cont_init()` | `src/probability_distributions.c` |

## L4: Conservation Laws / Theorems (Complete)

| # | Concept | C Implementation | Location |
|---|---------|-----------------|----------|
| 1 | Total probability theorem | `prob_total_probability_check()` | `src/probability_stats.c` |
| 2 | Law of total expectation | `prob_total_expectation_check()` | `src/probability_stats.c` |
| 3 | Law of large numbers (LLN) | `prob_lln_check()` | `src/probability_stats.c` |
| 4 | Central limit theorem (CLT) | `prob_clt_ks_test()` | `src/probability_stats.c` |
| 5 | Chapman-Kolmogorov equation | `markov_dtmc_nstep()` | `src/stochastic_process.c` |
| 6 | Global balance (CTMC stationarity) | `markov_ctmc_stationary()` | `src/stochastic_process.c` |
| 7 | Gibbs inequality (KL >= 0) | `prob_kl_divergence()` | `src/probability_information.c` |

## L5: Engineering Methods (Complete)

| # | Concept | C Implementation | Location |
|---|---------|-----------------|----------|
| 1 | Maximum likelihood estimation (Gaussian) | `prob_mle_gaussian()` | `src/probability_estimation.c` |
| 2 | MLE (exponential) | `prob_mle_exponential()` | `src/probability_estimation.c` |
| 3 | MLE (Poisson) | `prob_mle_poisson()` | `src/probability_estimation.c` |
| 4 | Kolmogorov-Smirnov test | `prob_ks_test()` | `src/probability_estimation.c` |
| 5 | Welford's online algorithm (mean/var) | `prob_empirical_insert()` | `src/probability_estimation.c` |
| 6 | Terriberry extension (skew/kurt online) | `prob_empirical_insert()` | `src/probability_estimation.c` |
| 7 | Hyndman-Fan Type 7 quantile | `prob_empirical_quantile()` | `src/probability_estimation.c` |
| 8 | Linear Kalman filter | `kalman_predict()`, `kalman_update()` | `src/estimation_filtering.c` |
| 9 | Kalman smoother (RTS) | `kalman_smoother()` | `src/estimation_filtering.c` |
| 10 | Extended Kalman filter (EKF) | `ekf_predict()`, `ekf_update()` | `src/estimation_filtering.c` |
| 11 | Unscented Kalman filter (UKF) | `ukf_step()` | `src/estimation_filtering.c` |
| 12 | Particle filter (SIR) | `particle_filter_sir()` | `src/estimation_filtering.c` |
| 13 | Monte Carlo integration (1D, ND) | `mc_integrate_1d()`, `mc_integrate_nd()` | `src/monte_carlo.c` |
| 14 | Importance sampling | `mc_importance_sampling()` | `src/monte_carlo.c` |
| 15 | Rejection sampling | `mc_rejection_sampling()` | `src/monte_carlo.c` |
| 16 | Metropolis-Hastings (1D, ND) | `mc_metropolis_hastings()`, `mc_metropolis_hastings_nd()` | `src/monte_carlo.c` |
| 17 | Gibbs sampling (2-component) | `mc_gibbs_sampler_2d()` | `src/monte_carlo.c` |
| 18 | Antithetic variates | `mc_antithetic_variates()` | `src/monte_carlo.c` |
| 19 | Control variates | `mc_control_variates()` | `src/monte_carlo.c` |
| 20 | Quasi-Monte Carlo (Halton, Sobol) | `qmc_halton_sequence()`, `qmc_sobol_sequence()` | `src/monte_carlo.c` |
| 21 | Shannon entropy | `prob_entropy_disc()` | `src/probability_information.c` |
| 22 | Differential entropy | `prob_differential_entropy()` | `src/probability_information.c` |
| 23 | KL divergence | `prob_kl_divergence()` | `src/probability_information.c` |
| 24 | Jensen-Shannon divergence | `prob_js_divergence()` | `src/probability_information.c` |
| 25 | Mutual information | `prob_mutual_information()` | `src/probability_information.c` |
| 26 | PRNG (LCG64, xorshift64, Box-Muller) | `prob_lcg64_uniform()`, etc. | `src/probability_core.c` |
| 27 | Gamma sampling (Marsaglia-Tsang) | `prob_sample_gamma_marsaglia()` | `src/probability_core.c` |
| 28 | AIC / BIC model selection | `info_criterion_aic()`, `info_criterion_bic()` | `src/estimation_filtering.c` |

## L6: Engineering Problems (Complete)

| # | Problem | Solution | Location |
|---|---------|----------|----------|
| 1 | Distribution fitting (MLE) | `prob_mle_gaussian()` etc. | `src/probability_estimation.c` |
| 2 | Goodness-of-fit testing | `prob_ks_test()` | `src/probability_estimation.c` |
| 3 | Online streaming statistics | `prob_empirical_*()` | `src/probability_estimation.c` |
| 4 | State estimation from noisy observations | `kalman_filter_run()` | `src/estimation_filtering.c` |
| 5 | Retrospective state smoothing | `kalman_smoother()` | `src/estimation_filtering.c` |
| 6 | Parameter estimation via likelihood | `kalman_log_likelihood()` | `src/estimation_filtering.c` |
| 7 | Weather forecasting (Markov chain) | `markov_dtmc_*()` | Example `ex02_markov_chain.c` |
| 8 | Radar target tracking | Kalman filter | Example `ex03_kalman_tracking.c` |
| 9 | Quality control (Gaussian MLE + bootstrap) | `prob_bootstrap_ci()` + MLE | Example `ex01_gaussian_mle.c` |
| 10 | Markov chain absorption analysis | `markov_absorption_prob()`, `markov_absorption_time()` | `src/stochastic_process.c` |
| 11 | Time-series autocorrelation analysis | `time_series_autocorr()` | `src/stochastic_process.c` |
| 12 | Spectral analysis (periodogram) | `periodogram_psd()` | `src/stochastic_process.c` |

## L7: Applications (Complete — 5 applications)

| # | Application | Implementation | Location |
|---|------------|---------------|----------|
| 1 | Radar/sonar tracking (Kalman) | `kalman_filter_run()` + Example 3 | `src/estimation_filtering.c` |
| 2 | GPS/IMU sensor fusion | `ukf_step()` (nonlinear filtering) | `src/estimation_filtering.c` |
| 3 | Quality control & process monitoring | MLE + KS test + bootstrap | `src/probability_estimation.c` |
| 4 | Weather risk modeling | Markov chain with MFPT analysis | Example 2 |
| 5 | Financial volatility estimation | `geometric_brownian_motion()` | `src/stochastic_process.c` |
| 6 | Bayesian inference via MCMC | `mc_metropolis_hastings()` | `src/monte_carlo.c` |

## L8: Advanced Methods (Complete — 4 topics)

| # | Topic | Implementation | Location |
|---|-------|---------------|----------|
| 1 | MCMC convergence diagnostics (Gelman-Rubin) | `mcmc_gelman_rubin()` | `src/monte_carlo.c` |
| 2 | Effective sample size (ESS) for MCMC | `mcmc_effective_sample_size()` | `src/monte_carlo.c` |
| 3 | Quasi-Monte Carlo (low-discrepancy sequences) | `qmc_halton_sequence()`, `qmc_sobol_sequence()` | `src/monte_carlo.c` |
| 4 | Unscented Kalman filter (sigma-point method) | `ukf_step()` | `src/estimation_filtering.c` |
| 5 | Particle filter for non-Gaussian/nonlinear systems | `particle_filter_sir()` | `src/estimation_filtering.c` |
| 6 | Importance sampling with self-normalised weights | `mc_importance_sampling()` | `src/monte_carlo.c` |

## L9: Research Frontiers (Partial — documented)

| # | Topic | Status |
|---|-------|--------|
| 1 | Hamiltonian Monte Carlo (HMC) | Documented, not implemented |
| 2 | Sequential Monte Carlo (SMC) beyond SIR | Documented, not implemented |
| 3 | Nonparametric Bayesian methods (Gaussian processes) | Documented, not implemented |
| 4 | Rough volatility models (fractional Brownian motion) | Documented, not implemented |
| 5 | Causal inference with do-calculus | Documented, not implemented |
