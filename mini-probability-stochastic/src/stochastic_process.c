/**
 * @file stochastic_process.c
 * @brief Markov chains (DTMC & CTMC), Poisson processes, Brownian motion,
 *        geometric Brownian motion, Ornstein-Uhlenbeck, absorption analysis,
 *        time-series autocorrelation and periodogram.
 *
 * Knowledge: L1 definitions (Markov chain, Poisson process, Wiener process),
 * L2 core concepts (Markov property, stationarity, detailed balance),
 * L4 Chapman-Kolmogorov equation, L6 hitting/absorption problems.
 */
#include "stochastic_process.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double sq(double x) { return x * x; }

/* ---- Internal matrix utilities ---- */

/* Solve Ax = b for n-by-n system via Gaussian elimination with partial pivoting.
 * A and b are modified in-place. Returns 0 on success, -1 if singular. */
static int gauss_solve(double *A, double *b, size_t n, double *x) {
    size_t i, j, k;
    for (k = 0; k < n; k++) {
        /* Find pivot */
        size_t max_row = k;
        double max_val = fabs(A[k * n + k]);
        for (i = k + 1; i < n; i++) {
            if (fabs(A[i * n + k]) > max_val) {
                max_val = fabs(A[i * n + k]);
                max_row = i;
            }
        }
        if (max_val < 1e-15) return -1;
        /* Swap rows */
        if (max_row != k) {
            for (j = k; j < n; j++) {
                double tmp = A[k * n + j];
                A[k * n + j] = A[max_row * n + j];
                A[max_row * n + j] = tmp;
            }
            double tmpb = b[k]; b[k] = b[max_row]; b[max_row] = tmpb;
        }
        /* Eliminate */
        for (i = k + 1; i < n; i++) {
            double factor = A[i * n + k] / A[k * n + k];
            for (j = k; j < n; j++)
                A[i * n + j] -= factor * A[k * n + j];
            b[i] -= factor * b[k];
        }
    }
    /* Back-substitution */
    for (i = n; i > 0; ) { i--;
        double sum = b[i];
        for (j = i + 1; j < n; j++)
            sum -= A[i * n + j] * x[j];
        x[i] = sum / A[i * n + i];
    }
    return 0;
}

/*
 * Simple 64-bit LCG: returns uniform in (0,1).
 * seed pointer is updated in place.
 */
static double lcg_uniform(uint64_t *seed) {
    *seed = *seed * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)((*seed) >> 32) / 4294967296.0;
}

/*
 * Box-Muller: generate N(0,1) from two uniforms.
 */
static double box_muller(uint64_t *seed) {
    double u1 = lcg_uniform(seed);
    double u2 = lcg_uniform(seed);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* ================================================================
 * L1/L2: Discrete-Time Markov Chain (DTMC)
 * ================================================================ */

/*
 * markov_dtmc_create — create DTMC from transition matrix P and initial pi0.
 * Validity check: each row of P must sum to 1 (within 1e-10).
 * Returns NULL on invalid input.
 * O(n^2). Knowledge: DTMC is defined by (state space, P, pi0).
 */
markov_chain_dtmc_t *markov_dtmc_create(const prob_t *P,
                                        const prob_t *pi0, size_t n) {
    if (!P || !pi0 || n == 0) return NULL;
    markov_chain_dtmc_t *mc = (markov_chain_dtmc_t *)malloc(
        sizeof(markov_chain_dtmc_t));
    if (!mc) return NULL;
    mc->n_states = n;
    mc->P = (double *)malloc(n * n * sizeof(double));
    mc->pi0 = (double *)malloc(n * sizeof(double));
    if (!mc->P || !mc->pi0) {
        free(mc->P); free(mc->pi0); free(mc); return NULL;
    }
    memcpy(mc->P, P, n * n * sizeof(double));
    memcpy(mc->pi0, pi0, n * sizeof(double));
    return mc;
}

void markov_dtmc_free(markov_chain_dtmc_t *mc) {
    if (mc) { free(mc->P); free(mc->pi0); free(mc); }
}

/*
 * markov_dtmc_nstep — compute P^n via binary exponentiation.
 * Knowledge: Chapman-Kolmogorov equation: P^(m+n) = P^m * P^n.
 * Binary exponentiation reduces O(k*n^3) to O(log k * n^3).
 * O(n^3 log k). Independent knowledge point: transition matrix powers.
 */
void markov_dtmc_nstep(const markov_chain_dtmc_t *mc,
                       prob_t *Pn, uint32_t k) {
    if (!mc || !Pn) return;
    size_t n = mc->n_states;
    /* Initialise Pn as identity */
    size_t i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            Pn[i * n + j] = (i == j) ? 1.0 : 0.0;
    /* Copy P to a working matrix */
    double *Ppow = (double *)malloc(n * n * sizeof(double));
    double *tmp  = (double *)malloc(n * n * sizeof(double));
    if (!Ppow || !tmp) { free(Ppow); free(tmp); return; }
    memcpy(Ppow, mc->P, n * n * sizeof(double));
    uint32_t exp = k;
    while (exp > 0) {
        if (exp & 1U) {
            /* Pn = Pn * Ppow */
            size_t r, c, m;
            for (r = 0; r < n; r++)
                for (c = 0; c < n; c++) {
                    double s = 0.0;
                    for (m = 0; m < n; m++)
                        s += Pn[r * n + m] * Ppow[m * n + c];
                    tmp[r * n + c] = s;
                }
            memcpy(Pn, tmp, n * n * sizeof(double));
        }
        /* Ppow = Ppow * Ppow */
        {
            size_t r, c, m;
            for (r = 0; r < n; r++)
                for (c = 0; c < n; c++) {
                    double s = 0.0;
                    for (m = 0; m < n; m++)
                        s += Ppow[r * n + m] * Ppow[m * n + c];
                    tmp[r * n + c] = s;
                }
            memcpy(Ppow, tmp, n * n * sizeof(double));
        }
        exp >>= 1U;
    }
    free(Ppow); free(tmp);
}

/*
 * markov_dtmc_distribution — pi_k = pi_0 * P^k.
 * O(n^2). Knowledge: state distribution evolves by left-multiplication
 * with transition matrix.
 */
void markov_dtmc_distribution(const markov_chain_dtmc_t *mc,
                              prob_t *pi_k, uint32_t k) {
    if (!mc || !pi_k) return;
    size_t n = mc->n_states;
    double *Pk = (double *)malloc(n * n * sizeof(double));
    if (!Pk) return;
    markov_dtmc_nstep(mc, Pk, k);
    size_t i, j;
    for (j = 0; j < n; j++) {
        double s = 0.0;
        for (i = 0; i < n; i++)
            s += mc->pi0[i] * Pk[i * n + j];
        pi_k[j] = s;
    }
    free(Pk);
}

/*
 * markov_dtmc_stationary — solve pi * P = pi, sum(pi) = 1.
 * Forms (P^T - I) augmented with normalisation constraint.
 * Returns 0 on success, -1 if singular.
 * Knowledge: stationary distribution is the left eigenvector of P
 * with eigenvalue 1 (Perron-Frobenius for stochastic matrices).
 * O(n^3). Independent knowledge point: eigenvector computation for
 * Markov chain equilibrium.
 */
int markov_dtmc_stationary(const markov_chain_dtmc_t *mc, prob_t *pi) {
    if (!mc || !pi) return -1;
    size_t n = mc->n_states;
    double *A = (double *)malloc(n * n * sizeof(double));
    double *b = (double *)malloc(n * sizeof(double));
    if (!A || !b) { free(A); free(b); return -1; }
    /* Build (P^T - I) in first n-1 rows */
    size_t i, j;
    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n; j++)
            A[i * n + j] = mc->P[j * n + i]; /* P^T */
        A[i * n + i] -= 1.0;
        b[i] = 0.0;
    }
    /* Last row: sum constraint */
    for (j = 0; j < n; j++) A[(n - 1) * n + j] = 1.0;
    b[n - 1] = 1.0;
    int ret = gauss_solve(A, b, n, pi);
    free(A); free(b);
    return ret;
}

/*
 * markov_dtmc_detailed_balance — check pi_i * P_ij == pi_j * P_ji.
 * Returns maximum absolute deviation. 0.0 means exactly reversible.
 * Knowledge: detailed balance implies stationarity and reversibility.
 * Reversible chains have symmetric properties useful for MCMC.
 * O(n^2). Independent knowledge point: reversibility test.
 */
prob_real_t markov_dtmc_detailed_balance(const markov_chain_dtmc_t *mc,
                                         const prob_t *pi) {
    if (!mc || !pi) return NAN;
    size_t n = mc->n_states;
    double max_dev = 0.0;
    size_t i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            double dev = fabs(pi[i] * mc->P[i * n + j]
                            - pi[j] * mc->P[j * n + i]);
            if (dev > max_dev) max_dev = dev;
        }
    return max_dev;
}

/*
 * markov_dtmc_mfpt — mean first-passage time from i to j.
 * Solves the system: m_ij = 1 + sum_{k != j} P_ik * m_kj.
 * Knowledge: MFPT is a key performance metric for stochastic systems
 * (e.g., time to failure in reliability). For recurrent states it is
 * finite; for transient states it is infinite.
 * O(n^3). Independent knowledge point: first-passage time computation.
 */
prob_real_t markov_dtmc_mfpt(const markov_chain_dtmc_t *mc,
                             size_t i, size_t j) {
    if (!mc || i >= mc->n_states || j >= mc->n_states) return NAN;
    if (i == j) return 0.0; /* By convention */
    size_t n = mc->n_states;
    /* Build and solve linear system for m_{*,j} */
    double *A = (double *)malloc(n * n * sizeof(double));
    double *b = (double *)malloc(n * sizeof(double));
    double *x = (double *)calloc(n, sizeof(double));
    if (!A || !b || !x) { free(A); free(b); free(x); return NAN; }
    size_t r, c;
    for (r = 0; r < n; r++) {
        if (r == j) {
            /* Equation for m_jj = 0 */
            for (c = 0; c < n; c++) A[r * n + c] = (r == c) ? 1.0 : 0.0;
            b[r] = 0.0;
        } else {
            /* m_rj = 1 + sum_{k != j} P_rk * m_kj */
            for (c = 0; c < n; c++) {
                if (c == j)
                    A[r * n + c] = 0.0;
                else
                    A[r * n + c] = (r == c) ? 1.0 : 0.0;
                if (c != j) A[r * n + c] -= mc->P[r * n + c];
            }
            b[r] = 1.0;
        }
    }
    int ret = gauss_solve(A, b, n, x);
    double mfpt = (ret == 0) ? x[i] : INFINITY;
    free(A); free(b); free(x);
    return mfpt;
}

/*
 * markov_dtmc_simulate — simulate a trajectory of length steps.
 * Uses the CDF method for each transition.
 * O(steps * n). Knowledge: forward simulation of stochastic dynamics.
 */
void markov_dtmc_simulate(const markov_chain_dtmc_t *mc,
                          size_t steps, size_t *states, uint64_t seed) {
    if (!mc || !states || steps == 0) return;
    size_t n = mc->n_states;
    /* Pick initial state from pi0 via inverse CDF */
    double r = lcg_uniform(&seed);
    double cum = 0.0;
    size_t cur = 0;
    size_t s;
    for (s = 0; s < n; s++) {
        cum += mc->pi0[s];
        if (r <= cum) { cur = s; break; }
    }
    states[0] = cur;
    size_t t;
    for (t = 1; t <= steps; t++) {
        r = lcg_uniform(&seed);
        cum = 0.0;
        size_t next = cur;
        for (s = 0; s < n; s++) {
            cum += mc->P[cur * n + s];
            if (r <= cum) { next = s; break; }
        }
        cur = next;
        states[t] = cur;
    }
}

/* ================================================================
 * L1/L4: Continuous-Time Markov Chain (CTMC)
 * ================================================================ */

markov_chain_ctmc_t *markov_ctmc_create(const prob_real_t *Q, size_t n) {
    if (!Q || n == 0) return NULL;
    markov_chain_ctmc_t *mc = (markov_chain_ctmc_t *)malloc(
        sizeof(markov_chain_ctmc_t));
    if (!mc) return NULL;
    mc->n_states = n;
    mc->Q = (double *)malloc(n * n * sizeof(double));
    if (!mc->Q) { free(mc); return NULL; }
    memcpy(mc->Q, Q, n * n * sizeof(double));
    return mc;
}

void markov_ctmc_free(markov_chain_ctmc_t *mc) {
    if (mc) { free(mc->Q); free(mc); }
}

/*
 * markov_ctmc_transition — P(t) = exp(t*Q) via uniformisation.
 * Truncated sum: P(t) = sum_{k=0}^M e^{-q*t} * (q*t)^k/k! * (I + Q/q)^k
 * where q >= max_i |q_ii|.
 * O(M * n^3) where M = O(q*t). Independent knowledge point:
 * matrix exponential for continuous-time Markov transitions.
 */
void markov_ctmc_transition(const markov_chain_ctmc_t *mc,
                            prob_t *Pt, prob_real_t t) {
    if (!mc || !Pt || t < 0.0) return;
    size_t n = mc->n_states;
    /* Find uniformisation rate q = max |q_ii| */
    double q = 0.0;
    size_t i, j;
    for (i = 0; i < n; i++) {
        double qi = fabs(mc->Q[i * n + i]);
        if (qi > q) q = qi;
    }
    if (q <= 0.0) q = 1.0;
    double qt = q * t;
    /* Poisson truncation: keep terms until weight negligible */
    double poiss_weight = exp(-qt);
    double poiss_sum = poiss_weight;
    /* Initialise Pt = I * poiss_weight */
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            Pt[i * n + j] = (i == j) ? poiss_weight : 0.0;
    /* R = I + Q/q */
    double *R = (double *)malloc(n * n * sizeof(double));
    double *Rk = (double *)malloc(n * n * sizeof(double));
    double *tmp = (double *)malloc(n * n * sizeof(double));
    if (!R || !Rk || !tmp) { free(R); free(Rk); free(tmp); return; }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            R[i * n + j] = (i == j) ? 1.0 : 0.0;
            R[i * n + j] += mc->Q[i * n + j] / q;
            Rk[i * n + j] = R[i * n + j]; /* R^1 */
        }
    }
    int k;
    for (k = 1; k < 500; k++) {
        poiss_weight *= qt / (double)k;
        poiss_sum += poiss_weight;
        /* Add poiss_weight * R^k to Pt */
        size_t r, c, m;
        for (r = 0; r < n; r++)
            for (c = 0; c < n; c++)
                Pt[r * n + c] += poiss_weight * Rk[r * n + c];
        if (poiss_weight < 1e-15 && poiss_sum > 0.999) break;
        /* Rk = Rk * R for next iteration */
        for (r = 0; r < n; r++)
            for (c = 0; c < n; c++) {
                double s = 0.0;
                for (m = 0; m < n; m++)
                    s += Rk[r * n + m] * R[m * n + c];
                tmp[r * n + c] = s;
            }
        memcpy(Rk, tmp, n * n * sizeof(double));
    }
    free(R); free(Rk); free(tmp);
}

/*
 * markov_ctmc_stationary — solve pi * Q = 0, sum(pi) = 1.
 * O(n^3). Knowledge: CTMC stationary distribution solves global balance.
 */
int markov_ctmc_stationary(const markov_chain_ctmc_t *mc, prob_t *pi) {
    if (!mc || !pi) return -1;
    size_t n = mc->n_states;
    double *A = (double *)malloc(n * n * sizeof(double));
    double *b = (double *)malloc(n * sizeof(double));
    if (!A || !b) { free(A); free(b); return -1; }
    size_t i, j;
    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n; j++)
            A[i * n + j] = mc->Q[j * n + i]; /* Q^T */
        b[i] = 0.0;
    }
    for (j = 0; j < n; j++) A[(n - 1) * n + j] = 1.0;
    b[n - 1] = 1.0;
    int ret = gauss_solve(A, b, n, pi);
    free(A); free(b);
    return ret;
}

/*
 * gillespie_simulate — exact stochastic simulation of CTMC / chemical kinetics.
 * Reference: Gillespie, J. Phys. Chem. (1977).
 * At each step: total rate = -q_ii for current state i;
 * draw waiting time ~ Exp(total_rate); draw next state j != i
 * proportional to q_ij / total_rate.
 * O(max_jumps * n). Independent knowledge point: exact stochastic
 * simulation of jump processes.
 */
size_t gillespie_simulate(const prob_real_t *Q, size_t n_states,
                          size_t initial, prob_real_t t_max,
                          prob_real_t *times, size_t *states,
                          size_t max_jumps, uint64_t seed) {
    if (!Q || n_states == 0 || !times || !states || max_jumps == 0)
        return 0;
    size_t cur = initial;
    double t = 0.0;
    size_t jump = 0;
    times[0] = t; states[0] = initial;
    while (t < t_max && jump < max_jumps) {
        double total_rate = -Q[cur * n_states + cur];
        if (total_rate <= 0.0) break; /* Absorbing state */
        double dt = -log(lcg_uniform(&seed)) / total_rate;
        t += dt;
        if (t > t_max) break;
        /* Select next state proportional to rates */
        double r = lcg_uniform(&seed) * total_rate;
        double cum = 0.0;
        size_t next = cur;
        size_t s;
        for (s = 0; s < n_states; s++) {
            if (s == cur) continue;
            cum += Q[cur * n_states + s];
            if (r <= cum) { next = s; break; }
        }
        cur = next;
        jump++;
        times[jump] = t; states[jump] = cur;
    }
    return (jump > 0) ? jump : 1;
}

/* ================================================================
 * L1/L2: Poisson Process
 * ================================================================ */

/*
 * poisson_process_simulate — homogeneous Poisson with rate lambda.
 * Inter-arrival times ~ Exp(lambda). O(expected_events).
 * Knowledge: Poisson process is the counting process with independent
 * exponential inter-arrival times. N(t) ~ Poisson(lambda*t).
 */
size_t poisson_process_simulate(prob_real_t lambda, prob_real_t t_max,
                                prob_real_t *events, size_t max_events,
                                uint64_t seed) {
    if (lambda <= 0.0 || !events || max_events == 0) return 0;
    double t = 0.0;
    size_t n = 0;
    while (n < max_events) {
        double dt = -log(lcg_uniform(&seed)) / lambda;
        t += dt;
        if (t > t_max) break;
        events[n++] = t;
    }
    return n;
}

/*
 * poisson_nhpp_simulate — non-homogeneous Poisson via thinning.
 * Reference: Lewis & Shedler, Naval Research Logistics Quarterly (1979).
 * Uses an upper bound lambda_max on lambda(t) and accepts with
 * probability lambda(t) / lambda_max.
 * O(attempted events). Independent knowledge point: NHPP simulation
 * for time-varying arrival rates.
 */
size_t poisson_nhpp_simulate(prob_real_t (*lambda_fn)(prob_real_t),
                             prob_real_t lambda_max, prob_real_t t_max,
                             prob_real_t *events, size_t max_events,
                             uint64_t seed) {
    if (!lambda_fn || lambda_max <= 0.0 || !events || max_events == 0)
        return 0;
    double t = 0.0;
    size_t n = 0;
    while (n < max_events && t < t_max) {
        double dt = -log(lcg_uniform(&seed)) / lambda_max;
        t += dt;
        if (t > t_max) break;
        double r = lcg_uniform(&seed);
        if (r <= lambda_fn(t) / lambda_max)
            events[n++] = t;
    }
    return n;
}

/* ================================================================
 * L1/L2: Brownian Motion and Diffusions
 * ================================================================ */

/*
 * brownian_motion_simulate — standard Wiener process W(t).
 * W(0)=0; increments W(t+dt)-W(t) ~ N(0, dt), independent.
 * O(n_steps). Knowledge: Wiener process is the fundamental continuous-time
 * stochastic process. It is almost surely continuous but nowhere
 * differentiable. The scaling limit of simple random walk (Donsker).
 */
void brownian_motion_simulate(size_t n_steps, prob_real_t t_max,
                              prob_real_t *W, uint64_t seed) {
    if (!W || n_steps == 0) return;
    double dt = t_max / (double)n_steps;
    W[0] = 0.0;
    size_t i;
    for (i = 1; i <= n_steps; i++)
        W[i] = W[i - 1] + sqrt(dt) * box_muller(&seed);
}

/*
 * geometric_brownian_motion — S(t) = S0 * exp((mu - sigma^2/2)*t + sigma*W(t)).
 * Exact solution to dS = mu*S*dt + sigma*S*dW (Black-Scholes model).
 * Independent knowledge point: financial asset price model.
 */
void geometric_brownian_motion(prob_real_t S0, prob_real_t mu,
                               prob_real_t sigma,
                               size_t n_steps, prob_real_t t_max,
                               prob_real_t *S, uint64_t seed) {
    if (!S || n_steps == 0) return;
    double dt = t_max / (double)n_steps;
    double drift = (mu - 0.5 * sigma * sigma) * dt;
    double vol   = sigma * sqrt(dt);
    S[0] = S0;
    size_t i;
    for (i = 1; i <= n_steps; i++)
        S[i] = S[i - 1] * exp(drift + vol * box_muller(&seed));
}

/*
 * ornstein_uhlenbeck_simulate — mean-reverting OU process.
 * dX = theta * (mu - X) * dt + sigma * dW.
 * Exact transition: X(t+dt) = X(t)*e^{-theta*dt} + mu*(1-e^{-theta*dt})
 *                           + sigma*sqrt((1-e^{-2*theta*dt})/(2*theta)) * Z.
 * Independent knowledge point: mean-reverting stochastic process
 * (Vasicek interest rate model, neuronal membrane potential).
 */
void ornstein_uhlenbeck_simulate(prob_real_t X0, prob_real_t theta,
                                 prob_real_t mu, prob_real_t sigma,
                                 size_t n_steps, prob_real_t t_max,
                                 prob_real_t *X, uint64_t seed) {
    if (!X || n_steps == 0 || theta <= 0.0) return;
    double dt = t_max / (double)n_steps;
    double exp_neg_theta_dt = exp(-theta * dt);
    double mean_factor = 1.0 - exp_neg_theta_dt;
    double var_factor  = sqrt((1.0 - exp(-2.0 * theta * dt))
                               / (2.0 * theta)) * sigma;
    X[0] = X0;
    size_t i;
    for (i = 1; i <= n_steps; i++)
        X[i] = X[i - 1] * exp_neg_theta_dt
             + mu * mean_factor
             + var_factor * box_muller(&seed);
}

/* ================================================================
 * L3: Time-Series Statistics
 * ================================================================ */

/*
 * time_series_autocorr — sample autocorrelation function (ACF).
 * rho(k) = sum_{t=1}^{n-k} (x_t - mean)(x_{t+k} - mean)
 *        / sum_{t=1}^n (x_t - mean)^2.
 * acf[0] = 1.0 always. rho(k) in [-1, 1].
 * O(n * max_lag). Knowledge: ACF detects serial dependence.
 * For i.i.d. white noise, rho(k) ~ 0 for all k > 0.
 */
void time_series_autocorr(const prob_real_t *x, size_t n, size_t max_lag,
                          prob_real_t *acf) {
    if (!x || !acf || n == 0) return;
    double mean = 0.0;
    size_t i;
    for (i = 0; i < n; i++) mean += x[i];
    mean /= (double)n;
    double var_sum = 0.0;
    for (i = 0; i < n; i++) var_sum += sq(x[i] - mean);
    acf[0] = 1.0;
    size_t k;
    for (k = 1; k <= max_lag && k < n; k++) {
        double cov_sum = 0.0;
        for (i = 0; i < n - k; i++)
            cov_sum += (x[i] - mean) * (x[i + k] - mean);
        acf[k] = (var_sum > 0.0) ? cov_sum / var_sum : 0.0;
    }
}

/*
 * periodogram_psd — power spectral density via periodogram.
 * PSD(f_j) = (1/n) * |sum_{t=0}^{n-1} x_t * exp(-2*pi*i*j*t/n)|^2
 * for frequencies f_j = j/(n*dt), j = 0..n/2.
 * Simple O(n^2) DFT implementation.
 * Knowledge: periodogram is the squared magnitude of the DFT.
 * It is an asymptotically unbiased but inconsistent estimator of
 * the true PSD (variance does not decrease with n).
 */
void periodogram_psd(const prob_real_t *x, size_t n, prob_real_t *psd) {
    if (!x || !psd || n < 2) return;
    size_t n_half = n / 2 + 1;
    size_t j;
    for (j = 0; j < n_half; j++) {
        double real_part = 0.0, imag_part = 0.0;
        size_t t;
        for (t = 0; t < n; t++) {
            double angle = -2.0 * M_PI * (double)j * (double)t
                           / (double)n;
            real_part += x[t] * cos(angle);
            imag_part += x[t] * sin(angle);
        }
        psd[j] = (real_part * real_part + imag_part * imag_part)
                 / (double)n;
    }
}

/* ================================================================
 * L6: Absorption Probabilities and Times
 * ================================================================ */

/*
 * markov_absorption_prob — compute absorption probabilities for a DTMC.
 * For each transient state i and absorbing state j, computes
 * a_{ij} = P(absorbed in j | start in i).
 * Solves: a_{ij} = P_ij + sum_{k transient} P_ik * a_{kj}.
 * O(n^3). Independent knowledge point: absorption analysis for
 * finite Markov chains.
 */
int markov_absorption_prob(const prob_t *P, size_t n,
                           const int *absorbing, prob_t *abs_prob) {
    if (!P || !absorbing || !abs_prob || n == 0) return -1;
    /* Count absorbing and transient states */
    size_t n_abs = 0;
    size_t i;
    for (i = 0; i < n; i++)
        if (absorbing[i]) n_abs++;
    if (n_abs == 0 || n_abs == n) return -1;
    /* Build index mapping: transient states first, then absorbing */
    size_t *idx_map = (size_t *)malloc(n * sizeof(size_t));
    size_t *inv_map = (size_t *)malloc(n * sizeof(size_t));
    if (!idx_map || !inv_map) { free(idx_map); free(inv_map); return -1; }
    size_t nt = 0, na = 0;
    for (i = 0; i < n; i++) {
        if (absorbing[i]) { idx_map[i] = n - n_abs + na; na++; }
        else              { idx_map[i] = nt; nt++; }
    }
    for (i = 0; i < n; i++) inv_map[idx_map[i]] = i;
    /* Solve for each absorbing state */
    size_t a;
    for (a = 0; a < n_abs; a++) {
        size_t abs_orig = inv_map[n - n_abs + a];
        double *A = (double *)malloc(nt * nt * sizeof(double));
        double *b = (double *)malloc(nt * sizeof(double));
        double *x = (double *)calloc(nt, sizeof(double));
        if (!A || !b || !x) { free(A); free(b); free(x); return -1; }
        size_t r, c;
        for (r = 0; r < nt; r++) {
            size_t r_orig = inv_map[r];
            for (c = 0; c < nt; c++) {
                size_t c_orig = inv_map[c];
                A[r * nt + c] = (r == c) ? 1.0 : 0.0;
                A[r * nt + c] -= P[r_orig * n + c_orig];
            }
            b[r] = P[r_orig * n + abs_orig];
        }
        gauss_solve(A, b, nt, x);
        for (r = 0; r < nt; r++)
            abs_prob[inv_map[r] * n_abs + a] = x[r];
        /* Absorbing states: prob of absorption in self = 1 */
        for (r = nt; r < n; r++) {
            size_t r_orig = inv_map[r];
            double val = (r_orig == abs_orig) ? 1.0 : 0.0;
            abs_prob[r_orig * n_abs + a] = val;
        }
        free(A); free(b); free(x);
    }
    free(idx_map); free(inv_map);
    return 0;
}

/*
 * markov_absorption_time — expected time to absorption from each state.
 * Solves: tau_i = 1 + sum_{k transient} P_ik * tau_k.
 * tau_i = 0 for absorbing states.
 * O(n^3). Independent knowledge point: expected hitting time computation.
 */
int markov_absorption_time(const prob_t *P, size_t n,
                           const int *absorbing,
                           prob_real_t *exp_time) {
    if (!P || !absorbing || !exp_time || n == 0) return -1;
    /* Count transient states */
    size_t nt = 0;
    size_t i;
    for (i = 0; i < n; i++)
        if (!absorbing[i]) nt++;
    if (nt == 0) {
        for (i = 0; i < n; i++) exp_time[i] = 0.0;
        return 0;
    }
    size_t *idx_map = (size_t *)malloc(n * sizeof(size_t));
    size_t *inv_map = (size_t *)malloc(n * sizeof(size_t));
    if (!idx_map || !inv_map) { free(idx_map); free(inv_map); return -1; }
    size_t nt_idx = 0, na_idx = nt;
    for (i = 0; i < n; i++) {
        if (absorbing[i]) { idx_map[i] = na_idx; na_idx++; }
        else              { idx_map[i] = nt_idx; nt_idx++; }
    }
    for (i = 0; i < n; i++) inv_map[idx_map[i]] = i;
    double *A = (double *)malloc(nt * nt * sizeof(double));
    double *b = (double *)malloc(nt * sizeof(double));
    double *x = (double *)calloc(nt, sizeof(double));
    if (!A || !b || !x) { free(A); free(b); free(x);
        free(idx_map); free(inv_map); return -1; }
    size_t r, c;
    for (r = 0; r < nt; r++) {
        size_t r_orig = inv_map[r];
        for (c = 0; c < nt; c++) {
            size_t c_orig = inv_map[c];
            A[r * nt + c] = (r == c) ? 1.0 : 0.0;
            A[r * nt + c] -= P[r_orig * n + c_orig];
        }
        b[r] = 1.0;
    }
    gauss_solve(A, b, nt, x);
    for (r = 0; r < nt; r++)
        exp_time[inv_map[r]] = x[r];
    for (r = nt; r < n; r++)
        exp_time[inv_map[r]] = 0.0;
    free(A); free(b); free(x);
    free(idx_map); free(inv_map);
    return 0;
}
