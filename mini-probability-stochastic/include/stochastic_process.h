/**
 * @file stochastic_process.h
 * @brief Stochastic process definitions: Markov chains, Poisson processes,
 *        Brownian motion, and related discrete-event systems.
 *
 * Knowledge coverage:
 *   L1 - Definitions:  Markov chain, transition matrix, stationary
 *                      distribution, Poisson process, Brownian motion
 *   L2 - Core Concepts: Markov property, stationarity, ergodicity,
 *                      memoryless property, hitting times
 *   L3 - Engineering Quantities: transition probabilities, hitting
 *                                probabilities, mean first-passage times
 *   L4 - Conservation Laws: Chapman-Kolmogorov equation, master equation,
 *                          detailed balance
 *
 * Reference texts:
 *   - Norris, "Markov Chains" (1997)
 *   - Karlin & Taylor, "A First Course in Stochastic Processes", 2nd ed.
 *   - Ross, "Stochastic Processes", 2nd ed. (1996)
 *   - Oksendal, "Stochastic Differential Equations", 6th ed. (2003)
 *
 * Course mapping:
 *   MIT 6.262 (Discrete Stochastic Processes)
 *   Stanford MS&E 321 (Stochastic Processes)
 *   Berkeley IEOR 263A (Applied Stochastic Processes)
 *   Cambridge Part III (Advanced Probability)
 */

#ifndef STOCHASTIC_PROCESS_H
#define STOCHASTIC_PROCESS_H

#include "probability_core.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L1: Markov chain definitions
 * ================================================================ */

/**
 * @brief Discrete-time Markov chain (DTMC) with finite state space.
 *
 * States are indexed 0..n_states-1.
 * Transition matrix P is n_states-by-n_states, row-stochastic.
 */
typedef struct {
    size_t       n_states;    /**< Number of states */
    prob_t      *P;           /**< n_states-by-n_states transition matrix,
                                   row-major: P[i*n_states + j] = P(i->j) */
    prob_t      *pi0;         /**< Initial distribution over states */
} markov_chain_dtmc_t;

/**
 * @brief Continuous-time Markov chain (CTMC) with finite state space.
 *
 * Defined by rate matrix Q (generator).
 * q_ij >= 0 for i != j, q_ii = -sum_{j!=i} q_ij.
 */
typedef struct {
    size_t       n_states;    /**< Number of states */
    prob_real_t *Q;           /**< n_states-by-n_states rate matrix,
                                   row-major. Q[i,j] = rate from i to j */
} markov_chain_ctmc_t;

/* ================================================================
 * L2: Markov chain operations
 * ================================================================ */

/**
 * @brief Create a DTMC from transition matrix and initial distribution.
 *
 * @param P    Row-stochastic transition matrix (n * n), row-major.
 * @param pi0  Initial distribution (n), must sum to 1.
 * @param n    Number of states.
 * @return Valid DTMC or NULL on invalid input.
 * Complexity: O(n^2).
 */
markov_chain_dtmc_t *markov_dtmc_create(const prob_t *P,
                                        const prob_t *pi0,
                                        size_t n);

/**
 * @brief Release a DTMC.
 */
void markov_dtmc_free(markov_chain_dtmc_t *mc);

/**
 * @brief Compute n-step transition matrix P^n.
 *
 * Uses binary exponentiation (repeated squaring) for O(n^3 log k).
 * Reference: Norris, Theorem 1.1.1 (Chapman-Kolmogorov).
 *
 * @param Pn  Output matrix (n-by-n, row-major, pre-allocated).
 * @param mc  Markov chain.
 * @param k   Number of steps.
 */
void markov_dtmc_nstep(const markov_chain_dtmc_t *mc,
                       prob_t *Pn, uint32_t k);

/**
 * @brief Compute the state distribution after k steps: pi_k = pi_0 * P^k.
 *
 * @param pi_k  Output distribution (n, pre-allocated).
 * Complexity: O(k * n^2) for small k; O(n^3 log k) if k large with
 *             binary exponentiation.
 */
void markov_dtmc_distribution(const markov_chain_dtmc_t *mc,
                              prob_t *pi_k, uint32_t k);

/**
 * @brief Find stationary distribution pi = pi * P by solving
 *        (P^T - I) * pi = 0 with sum(pi) = 1.
 *
 * Uses Gaussian elimination with partial pivoting.
 * Assumes irreducible, aperiodic chain (unique stationary distribution).
 * Reference: Norris, Theorem 1.7.7.
 *
 * @param pi  Output stationary distribution (n, pre-allocated).
 * @return 0 on success, -1 if singular/unable to solve.
 */
int markov_dtmc_stationary(const markov_chain_dtmc_t *mc, prob_t *pi);

/**
 * @brief Check detailed balance: pi_i * P_ij == pi_j * P_ji for all i, j.
 *
 * @return Maximum absolute deviation from detailed balance.
 *         0.0 means exactly reversible.
 */
prob_real_t markov_dtmc_detailed_balance(const markov_chain_dtmc_t *mc,
                                         const prob_t *pi);

/**
 * @brief Compute mean first-passage time from state i to state j.
 *
 * Solves the linear system for absorption probabilities.
 * Reference: Norris, Section 1.3.
 * Complexity: O(n^3) for Gaussian elimination.
 *
 * @return Expected number of steps to reach j from i.
 *         Returns INFINITY if j not reachable from i.
 */
prob_real_t markov_dtmc_mfpt(const markov_chain_dtmc_t *mc,
                             size_t i, size_t j);

/**
 * @brief Simulate a trajectory of the DTMC.
 *
 * @param steps    Number of steps to simulate.
 * @param states   Output array (steps+1, pre-allocated).
 * @param seed     RNG seed (0 for non-deterministic).
 * Complexity: O(steps).
 */
void markov_dtmc_simulate(const markov_chain_dtmc_t *mc,
                          size_t steps, size_t *states, uint64_t seed);

/* ================================================================
 * L1/L2: Continuous-time Markov chain operations
 * ================================================================ */

/**
 * @brief Create a CTMC from rate matrix Q.
 *
 * @param Q  Rate matrix (n * n, row-major). Must satisfy:
 *           Q[i][j] >= 0 for i != j, sum_j Q[i][j] = 0.
 * @param n  Number of states.
 * @return Valid CTMC or NULL.
 */
markov_chain_ctmc_t *markov_ctmc_create(const prob_real_t *Q, size_t n);

/**
 * @brief Release a CTMC.
 */
void markov_ctmc_free(markov_chain_ctmc_t *mc);

/**
 * @brief Compute transition probability matrix P(t) = exp(t * Q).
 *
 * Uses uniformisation (Jensen's method) with truncation.
 * Reference: Norris, Theorem 2.1.3.
 *
 * @param Pt  Output matrix (n-by-n, row-major, pre-allocated).
 * @param t   Time interval.
 */
void markov_ctmc_transition(const markov_chain_ctmc_t *mc,
                            prob_t *Pt, prob_real_t t);

/**
 * @brief Find stationary distribution for CTMC: pi * Q = 0.
 *
 * @return 0 on success, -1 on failure.
 */
int markov_ctmc_stationary(const markov_chain_ctmc_t *mc, prob_t *pi);

/* ================================================================
 * L4: Gillespie algorithm (stochastic simulation algorithm)
 * ================================================================ */

/**
 * @brief Gillespie direct method for exact stochastic simulation
 *        of chemical reaction networks / CTMC trajectories.
 *
 * Reference: Gillespie, "Exact Stochastic Simulation of Coupled
 * Chemical Reactions", J. Phys. Chem. (1977).
 *
 * @param Q          Rate matrix (n_states * n_states).
 * @param n_states   Number of states.
 * @param initial    Starting state index.
 * @param t_max      Maximum simulation time.
 * @param times      Output jump time array (allocated by caller, >= max_jumps).
 * @param states     Output state array (allocated by caller, >= max_jumps).
 * @param max_jumps  Maximum number of jumps to record.
 * @param seed       RNG seed.
 * @return Actual number of jumps recorded.
 */
size_t gillespie_simulate(const prob_real_t *Q, size_t n_states,
                          size_t initial,
                          prob_real_t t_max,
                          prob_real_t *times, size_t *states,
                          size_t max_jumps, uint64_t seed);

/* ================================================================
 * L1: Poisson process
 * ================================================================ */

/**
 * @brief Simulate a homogeneous Poisson process with rate lambda.
 *
 * Generates inter-arrival times from Exp(lambda).
 * Reference: Ross, Ch. 5.
 *
 * @param lambda  Arrival rate (> 0).
 * @param t_max   Maximum time.
 * @param events  Output array of event times (pre-allocated, >= max_events).
 * @param max_events  Max events to record.
 * @param seed    RNG seed.
 * @return Number of events generated in [0, t_max].
 */
size_t poisson_process_simulate(prob_real_t lambda, prob_real_t t_max,
                                prob_real_t *events, size_t max_events,
                                uint64_t seed);

/**
 * @brief Non-homogeneous Poisson process with rate function lambda(t).
 *
 * Uses thinning (Lewis & Shedler, 1979).
 * Reference: Ross, Ch. 5, Section 5.4.
 *
 * @param lambda_fn  Rate function lambda(t) >= 0.
 * @param lambda_max Maximum value of lambda(t) on [0, t_max].
 * @param t_max      Maximum time.
 * @param events     Output array.
 * @param max_events Max events.
 * @param seed       RNG seed.
 * @return Number of events generated.
 */
size_t poisson_nhpp_simulate(prob_real_t (*lambda_fn)(prob_real_t),
                             prob_real_t lambda_max, prob_real_t t_max,
                             prob_real_t *events, size_t max_events,
                             uint64_t seed);

/* ================================================================
 * L1: Brownian motion / Wiener process
 * ================================================================ */

/**
 * @brief Simulate a standard Brownian motion (Wiener process) path.
 *
 * Properties: W(0)=0, independent Gaussian increments,
 *             W(t) ~ N(0, t), continuous paths.
 * Reference: Oksendal, Ch. 2.
 *
 * @param n_steps  Number of time steps.
 * @param t_max    Total time.
 * @param W        Output array of length n_steps+1 (W[0]=0).
 * @param seed     RNG seed.
 */
void brownian_motion_simulate(size_t n_steps, prob_real_t t_max,
                              prob_real_t *W, uint64_t seed);

/**
 * @brief Simulate geometric Brownian motion:
 *        dS = mu * S * dt + sigma * S * dW.
 *
 * Exact solution: S(t) = S0 * exp((mu - sigma^2/2)*t + sigma*W(t)).
 * Used in Black-Scholes option pricing.
 *
 * @param S0       Initial value.
 * @param mu       Drift.
 * @param sigma    Volatility (> 0).
 * @param n_steps  Time steps.
 * @param t_max    Total time.
 * @param S        Output array of length n_steps+1.
 * @param seed     RNG seed.
 */
void geometric_brownian_motion(prob_real_t S0, prob_real_t mu,
                               prob_real_t sigma,
                               size_t n_steps, prob_real_t t_max,
                               prob_real_t *S, uint64_t seed);

/**
 * @brief Simulate Ornstein-Uhlenbeck process:
 *        dX = theta * (mu - X) * dt + sigma * dW.
 *
 * Mean-reverting process. Exact transition density known.
 * Reference: Oksendal, Ch. 5.
 *
 * @param X0       Initial value.
 * @param theta    Mean reversion speed (> 0).
 * @param mu       Long-term mean.
 * @param sigma    Volatility.
 * @param n_steps  Time steps.
 * @param t_max    Total time.
 * @param X        Output array of length n_steps+1.
 * @param seed     RNG seed.
 */
void ornstein_uhlenbeck_simulate(prob_real_t X0, prob_real_t theta,
                                 prob_real_t mu, prob_real_t sigma,
                                 size_t n_steps, prob_real_t t_max,
                                 prob_real_t *X, uint64_t seed);

/* ================================================================
 * L3: Process statistics
 * ================================================================ */

/**
 * @brief Compute sample autocorrelation function for a time series.
 *
 * rho(k) = sum_{t=1}^{N-k} (x_t - mu)(x_{t+k} - mu) /
 *          sum_{t=1}^N (x_t - mu)^2
 *
 * Reference: Box, Jenkins, Reinsel, "Time Series Analysis", 4th ed.
 *
 * @param x      Time series data.
 * @param n      Length of series.
 * @param max_lag Maximum lag to compute.
 * @param acf    Output autocorrelation array (length max_lag+1, pre-allocated).
 *               acf[0] = 1.0 always.
 */
void time_series_autocorr(const prob_real_t *x, size_t n, size_t max_lag,
                          prob_real_t *acf);

/**
 * @brief Compute power spectral density via periodogram.
 *
 * Uses the squared magnitude of the DFT (simple O(n^2) implementation).
 * Reference: Oppenheim & Schafer, "Discrete-Time Signal Processing".
 *
 * @param x     Time series.
 * @param n     Length.
 * @param psd   Output PSD array (n/2+1 frequencies, pre-allocated).
 */
void periodogram_psd(const prob_real_t *x, size_t n, prob_real_t *psd);

/* ================================================================
 * L6: Hitting probabilities and absorbing Markov chains
 * ================================================================ */

/**
 * @brief Compute absorption probabilities for a DTMC with absorbing states.
 *
 * @param P            Transition matrix (n * n).
 * @param n            Total number of states.
 * @param absorbing    Boolean array: absorbing[i] = 1 if state i is absorbing.
 * @param abs_prob     Output n-by-n_absorbing matrix (row-major).
 *                     abs_prob[i, j] = P(absorbed in state j | start in i).
 * @return 0 on success, -1 on failure.
 */
int markov_absorption_prob(const prob_t *P, size_t n,
                           const int *absorbing,
                           prob_t *abs_prob);

/**
 * @brief Compute expected time until absorption from each transient state.
 *
 * @param P          Transition matrix (n * n).
 * @param n          Total states.
 * @param absorbing  Boolean array.
 * @param exp_time   Output array (n), exp_time[i] = E[time to absorption | X_0=i].
 *                   Exp_time[i] = 0 for absorbing states.
 * @return 0 on success.
 */
int markov_absorption_time(const prob_t *P, size_t n,
                           const int *absorbing,
                           prob_real_t *exp_time);

#ifdef __cplusplus
}
#endif

#endif /* STOCHASTIC_PROCESS_H */
