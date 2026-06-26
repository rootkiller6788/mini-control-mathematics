/**
 * @file estimation_filtering.h
 * @brief State estimation and filtering: Kalman filter (linear, extended,
 *        unscented), particle filter (SIR), and maximum likelihood estimation.
 *
 * Knowledge coverage:
 *   L1 - Definitions:  state-space model, observation model, filter density
 *   L2 - Core Concepts: Bayes filter, prediction-update cycle,
 *                      recursive least squares
 *   L5 - Engineering Methods: Kalman filter, EKF, UKF, particle filter
 *   L6 - Engineering Problems: target tracking, sensor fusion, navigation
 *   L7 - Applications: GPS/IMU integration, radar tracking, financial
 *                      volatility estimation
 *
 * Reference texts:
 *   - Kalman, "A New Approach to Linear Filtering and Prediction
 *     Problems", J. Basic Eng. (1960)
 *   - Anderson & Moore, "Optimal Filtering" (1979)
 *   - Julier & Uhlmann, "Unscented Filtering and Nonlinear Estimation"
 *     Proc. IEEE (2004)
 *   - Arulampalam et al., "A Tutorial on Particle Filters", IEEE TSP (2002)
 *
 * Course mapping:
 *   MIT 6.432 (Stochastic Processes, Detection, and Estimation)
 *   Stanford EE 363 (Linear Dynamical Systems)
 *   Berkeley EE 221A (Linear Systems Theory)
 *   Cambridge 4F12 (Computer Vision ? Kalman/particle filters)
 */

#ifndef ESTIMATION_FILTERING_H
#define ESTIMATION_FILTERING_H

#include "probability_core.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L1: State-space model definitions
 * ================================================================ */

/**
 * @brief Linear Gaussian state-space model:
 *        x_k = F * x_{k-1} + B * u_k + w_k,   w_k ~ N(0, Q)
 *        y_k = H * x_k + v_k,                  v_k ~ N(0, R)
 */
typedef struct {
    size_t       state_dim;   /**< Dimension of state vector x */
    size_t       input_dim;   /**< Dimension of control input u (0 if none) */
    size_t       obs_dim;     /**< Dimension of observation y */

    prob_real_t *F;  /**< State transition matrix (state_dim * state_dim) */
    prob_real_t *B;  /**< Control matrix (state_dim * input_dim), NULL if no ctrl */
    prob_real_t *H;  /**< Observation matrix (obs_dim * state_dim) */
    prob_real_t *Q;  /**< Process noise covariance (state_dim * state_dim) */
    prob_real_t *R;  /**< Observation noise covariance (obs_dim * obs_dim) */
} kalman_model_t;

/**
 * @brief Nonlinear state-space model (for EKF/UKF):
 *        x_k = f(x_{k-1}, u_k) + w_k
 *        y_k = h(x_k) + v_k
 */
typedef struct {
    size_t state_dim;
    size_t obs_dim;

    /** State transition function f(x, u) -> x_next.
     *  x: [state_dim], u: [input_dim], output written to x_next[state_dim]. */
    void (*f)(const prob_real_t *x, const prob_real_t *u, prob_real_t *x_next);

    /** Observation function h(x) -> y.
     *  x: [state_dim], output written to y[obs_dim]. */
    void (*h)(const prob_real_t *x, prob_real_t *y);

    prob_real_t *Q;  /**< Process noise covariance */
    prob_real_t *R;  /**< Observation noise covariance */
} nlinear_model_t;

/* ================================================================
 * L5: Linear Kalman filter
 * ================================================================ */

/**
 * @brief Kalman filter single predict step.
 *
 * x_pred = F * x + B * u
 * P_pred = F * P * F^T + Q
 *
 * Complexity: O(d^3) for matrix multiplication (d = state_dim).
 */
void kalman_predict(const kalman_model_t *model,
                    const prob_real_t *x, const prob_real_t *P,
                    const prob_real_t *u, /* can be NULL */
                    prob_real_t *x_pred, prob_real_t *P_pred);

/**
 * @brief Kalman filter single update (correction) step.
 *
 * K = P_pred * H^T * (H * P_pred * H^T + R)^{-1}   (Kalman gain)
 * x = x_pred + K * (y - H * x_pred)                 (state update)
 * P = (I - K * H) * P_pred                          (covariance update)
 *
 * Complexity: O(d^3 + d^2*m + d*m^2 + m^3) where d=state_dim, m=obs_dim.
 */
void kalman_update(const kalman_model_t *model,
                   const prob_real_t *x_pred, const prob_real_t *P_pred,
                   const prob_real_t *y,
                   prob_real_t *x, prob_real_t *P,
                   prob_real_t *K /* output Kalman gain, can be NULL */);

/**
 * @brief Run Kalman filter over a sequence of observations.
 *
 * Full predict-update cycle for each time step.
 * Complexity: O(T * d^3) where T is number of time steps.
 *
 * @param observations  T * obs_dim, row-major.
 * @param controls      T * input_dim, row-major (can be NULL).
 * @param T             Number of time steps.
 * @param x_init        Initial state.
 * @param P_init        Initial covariance.
 * @param x_estimates   Output: T * state_dim, row-major.
 * @param P_estimates   Output: T * state_dim * state_dim (flattened).
 */
void kalman_filter_run(const kalman_model_t *model,
                       const prob_real_t *observations,
                       const prob_real_t *controls,
                       size_t T,
                       const prob_real_t *x_init,
                       const prob_real_t *P_init,
                       prob_real_t *x_estimates,
                       prob_real_t *P_estimates);

/**
 * @brief Kalman smoother (Rauch-Tung-Striebel backward pass).
 *
 * Computes smoothed estimates x_{k|T}, P_{k|T} for k = 0..T-1.
 * Requires forward Kalman filter output.
 *
 * Reference: Rauch, Tung, Striebel, "Maximum Likelihood Estimates
 *            of Linear Dynamic Systems", AIAA J. (1965).
 * Complexity: O(T * d^3).
 */
void kalman_smoother(const kalman_model_t *model,
                     const prob_real_t *x_filtered,
                     const prob_real_t *P_filtered,
                     const prob_real_t *x_predicted,
                     const prob_real_t *P_predicted,
                     size_t T,
                     prob_real_t *x_smoothed,
                     prob_real_t *P_smoothed);

/* ================================================================
 * L5: Extended Kalman filter (EKF)
 * ================================================================ */

/**
 * @brief EKF predict step using Jacobian linearisation.
 *
 * Computes Jacobians of f and h via finite differences.
 * Complexity: O(d^3 + d^2 * input_dim).
 */
void ekf_predict(const nlinear_model_t *model,
                 const prob_real_t *x, const prob_real_t *P,
                 const prob_real_t *u,
                 prob_real_t *x_pred, prob_real_t *P_pred,
                 prob_real_t eps /* finite diff step */);

/**
 * @brief EKF update step.
 */
void ekf_update(const nlinear_model_t *model,
                const prob_real_t *x_pred, const prob_real_t *P_pred,
                const prob_real_t *y,
                prob_real_t *x, prob_real_t *P,
                prob_real_t eps);

/* ================================================================
 * L5: Unscented Kalman filter (UKF)
 * ================================================================ */

/**
 * @brief UKF single predict-update cycle.
 *
 * Uses 2n+1 sigma points for deterministic sampling.
 * Captures mean and covariance to 2nd order (3rd for Gaussians).
 *
 * Reference: Julier & Uhlmann (2004).
 * Complexity: O(d^3) per step.
 *
 * @param alpha  Sigma-point spread parameter (typically 1e-3 to 1).
 * @param beta   Prior knowledge parameter (2 for Gaussian optimal).
 * @param kappa  Secondary scaling (typically 0 or 3-n).
 */
void ukf_step(const nlinear_model_t *model,
              const prob_real_t *x, const prob_real_t *P,
              const prob_real_t *u, const prob_real_t *y,
              prob_real_t *x_out, prob_real_t *P_out,
              prob_real_t alpha, prob_real_t beta, prob_real_t kappa);

/* ================================================================
 * L5: Particle filter (Sequential Importance Resampling)
 * ================================================================ */

/**
 * @brief Bootstrap particle filter (SIR) for a single time step.
 *
 * Reference: Gordon, Salmond, Smith, "Novel Approach to
 *            Nonlinear/Non-Gaussian Bayesian State Estimation",
 *            IEE Proc-F (1993).
 *
 * @param n_particles  Number of particles.
 * @param state_dim    State dimension.
 * @param obs_dim      Observation dimension.
 *
 * @param particles      [n_particles * state_dim], in/out: particle states.
 * @param weights        [n_particles], in/out: importance weights.
 * @param f              State transition (with process noise).
 * @param obs_likelihood p(y | x): observation likelihood function.
 * @param u              Control input (can be NULL).
 * @param y              Current observation.
 * @param seed           RNG seed pointer.
 */
void particle_filter_sir(size_t n_particles, size_t state_dim,
                         size_t obs_dim,
                         prob_real_t *particles,
                         prob_real_t *weights,
                         void (*f)(const prob_real_t*, const prob_real_t*,
                                   prob_real_t*, uint64_t*),
                         prob_real_t (*obs_likelihood)(const prob_real_t*,
                                                       const prob_real_t*,
                                                       size_t),
                         const prob_real_t *u, const prob_real_t *y,
                         uint64_t *seed);

/**
 * @brief Estimate state from weighted particles (weighted mean).
 *
 * @param particles  [n_particles * state_dim].
 * @param weights    [n_particles], normalised.
 * @param state_dim  State dimension.
 * @param x_est      Output: [state_dim].
 */
void particle_filter_estimate(const prob_real_t *particles,
                              const prob_real_t *weights,
                              size_t n_particles, size_t state_dim,
                              prob_real_t *x_est);

/* ================================================================
 * L6: Likelihood computation and model selection
 * ================================================================ */

/**
 * @brief Compute log-likelihood of observations under a Kalman model.
 *
 * Useful for parameter estimation via maximum likelihood.
 * Reference: Shumway & Stoffer, "Time Series Analysis and Its
 *            Applications", 4th ed. (2017), Ch. 6.
 * Complexity: O(T * d^3).
 */
prob_real_t kalman_log_likelihood(const kalman_model_t *model,
                                  const prob_real_t *observations,
                                  size_t T,
                                  const prob_real_t *x_init,
                                  const prob_real_t *P_init);

/**
 * @brief Akaike Information Criterion for model selection.
 *
 * AIC = 2*k - 2*ln(L), where k = number of parameters,
 * L = maximised likelihood.
 * Lower AIC indicates better model.
 *
 * Reference: Akaike, "A New Look at the Statistical Model
 *            Identification", IEEE TAC (1974).
 */
prob_real_t info_criterion_aic(prob_real_t log_likelihood, size_t n_params);

/**
 * @brief Bayesian Information Criterion.
 *
 * BIC = k * ln(n) - 2 * ln(L), where n = number of observations.
 * Penalises complexity more heavily than AIC.
 *
 * Reference: Schwarz, "Estimating the Dimension of a Model",
 *            Annals of Statistics (1978).
 */
prob_real_t info_criterion_bic(prob_real_t log_likelihood, size_t n_params,
                               size_t n_observations);

/* ================================================================
 * Utility: matrix operations on prob_real_t
 * ================================================================ */

/**
 * @brief Matrix-matrix multiply: C = A * B (all square n-by-n).
 */
void mat_mul_n(const prob_real_t *A, const prob_real_t *B, prob_real_t *C,
               size_t n);

/**
 * @brief Matrix-vector multiply: y = A * x.
 */
void mat_vec_mul(const prob_real_t *A, const prob_real_t *x, prob_real_t *y,
                 size_t rows, size_t cols);

/**
 * @brief Matrix transpose: B = A^T.
 */
void mat_transpose(const prob_real_t *A, prob_real_t *B,
                   size_t rows, size_t cols);

/**
 * @brief Solve linear system Ax = b using Cholesky decomposition.
 *        A must be symmetric positive definite.
 *
 * @return 0 on success, -1 if non-positive-definite.
 */
int mat_cholesky_solve(const prob_real_t *A, const prob_real_t *b,
                       prob_real_t *x, size_t n);

/**
 * @brief Solve linear system using Gaussian elimination with partial pivoting.
 *
 * @return 0 on success, -1 if singular.
 */
int mat_gaussian_solve(prob_real_t *A, prob_real_t *b, size_t n,
                       prob_real_t *x);

/**
 * @brief Compute matrix inverse B = A^{-1} via Gaussian elimination.
 *
 * @return 0 on success, -1 if singular.
 */
int mat_inverse(const prob_real_t *A, prob_real_t *B, size_t n);

#ifdef __cplusplus
}
#endif

#endif /* ESTIMATION_FILTERING_H */
