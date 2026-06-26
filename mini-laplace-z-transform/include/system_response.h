/**
 * @file system_response.h
 * @brief Time-Domain System Response Analysis
 *
 * Knowledge Coverage: L5 (Algorithms), L6 (Canonical Problems)
 *
 * Core response types and specifications:
 *
 *   1. Step Response — system response to u(t) input
 *      Key specs: rise time (t_r), peak time (t_p), settling time (t_s),
 *                 overshoot (M_p), steady-state error (e_ss)
 *
 *   2. Impulse Response — system response to δ(t) input
 *      Complete characterization of LTI system. Convolution kernel.
 *
 *   3. Ramp Response — system response to t·u(t) input
 *      Measures tracking capability for time-varying references.
 *
 *   4. Transient Specification for 2nd-order systems (closed-form):
 *      - t_r ≈ 1.8/ωₙ (10%-90% rise time)
 *      - t_p = π/(ωₙ√(1-ζ²))
 *      - M_p = 100·exp(-πζ/√(1-ζ²)) %
 *      - t_s ≈ 4/(ζωₙ) (2% criterion)
 *      - For ζ < 1 (underdamped)
 *
 * References:
 *   - Ogata, "Modern Control Engineering" (5th ed.), Ch. 5
 *   - Franklin, Powell, Emami-Naeini, "Feedback Control", Ch. 3
 */

#ifndef SYSTEM_RESPONSE_H
#define SYSTEM_RESPONSE_H

#include <stddef.h>
#include "laplace_core.h"
#include "transfer_function.h"

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------
 * L1 — Response Data Types
 *----------------------------------------------------------------------*/

/**
 * @enum ResponseType
 * @brief  Types of test inputs for system characterization
 */
typedef enum {
    RESPONSE_STEP    = 0,    /**< Input = u(t) unit step               */
    RESPONSE_IMPULSE = 1,    /**< Input = δ(t) Dirac delta             */
    RESPONSE_RAMP    = 2,    /**< Input = t·u(t) unit ramp             */
    RESPONSE_SINUSOID = 3   /**< Input = sin(ωt)·u(t) sinusoidal      */
} ResponseType;

/**
 * @struct TransientSpecs
 * @brief  Standard transient response specifications
 *
 * L1: These are the universal metrics for time-domain control performance.
 */
typedef struct {
    double   rise_time;              /**< t_r: 10%-90% rise time (s)      */
    double   peak_time;              /**< t_p: time to first peak (s)     */
    double   settling_time;          /**< t_s: time to stay within 2% (s) */
    double   overshoot_percent;      /**< M_p: % overshoot above steady-state */
    double   steady_state_value;     /**< y_ss: final steady-state value  */
    double   steady_state_error;     /**< e_ss: reference - y_ss          */
    double   delay_time;             /**< t_d: 50% rise time (s)          */
    int      is_valid;               /**< 1 if computation converged      */
} TransientSpecs;

/**
 * @struct TimeResponse
 * @brief  Sampled time response of a system
 */
typedef struct {
    int      n_samples;              /**< Number of sample points         */
    double   t_start;                /**< Start time (s)                  */
    double   t_step;                 /**< Time step (s)                   */
    double  *time;                   /**< Time axis values (caller free)  */
    double  *output;                 /**< y(t) values (caller free)       */
    double  *input;                  /**< u(t) values (caller free)       */
} TimeResponse;

/*----------------------------------------------------------------------
 * L5 — Response Computation
 *----------------------------------------------------------------------*/

/**
 * @brief  Compute time response of a transfer function to a standard input
 * @param  tf           Transfer function G(s)
 * @param  resp_type    Input type (step, impulse, ramp, sinusoid)
 * @param  param        Extra parameter (frequency for sinusoid, ignored otherwise)
 * @param  T_end        Simulation end time (s)
 * @param  n_samples    Number of time samples
 * @param  response     Output time response (caller-allocated, arrays will be malloc'd)
 * @return 0 on success, -1 on numerical failure
 *
 * L5: Uses numerical inverse Laplace (Dubner-Abate-Crump) at each time point.
 *     For step response: Y(s) = G(s)/s
 *     For impulse: Y(s) = G(s)
 *     For ramp: Y(s) = G(s)/s²
 *     For sinusoid: Y(s) = G(s)·ω/(s²+ω²)
 *     Complexity: O(N·M) where N = n_samples, M = Fourier terms in inverse LT.
 */
int compute_time_response(const TransferFunction *tf,
                           ResponseType resp_type,
                           double param,
                           double T_end, int n_samples,
                           TimeResponse *response);

/**
 * @brief  Compute transient specifications from step response data
 * @param  response   Previously computed step response
 * @param  specs      Output transient specifications
 * @return 0 on success, -1 if data insufficient
 *
 * L5: Extracts t_r, t_p, t_s, M_p, e_ss from sampled step response.
 *     Uses interpolation for accurate crossing detection.
 *     Complexity: O(N) scan through response data.
 */
int extract_transient_specs(const TimeResponse *response,
                             TransientSpecs *specs);

/**
 * @brief  Compute closed-form 2nd-order specs directly from ωₙ and ζ
 * @param  omega_n      Natural frequency (rad/s)
 * @param  zeta         Damping ratio
 * @param  specs        Output specifications
 * @return 0 on success, -1 if ζ or ωₙ invalid
 *
 * L5: Analytical formulas for second-order underdamped systems:
 *     - t_r = (π - acos(ζ))/(ωₙ√(1-ζ²))   [0%-100%] or ≈ 1.8/ωₙ [10%-90%]
 *     - t_p = π/(ωₙ√(1-ζ²))
 *     - M_p = exp(-πζ/√(1-ζ²)) · 100%
 *     - t_s = 4/(ζωₙ) for 2% criterion
 *     Complexity: O(1).
 */
int second_order_specs_closed_form(double omega_n, double zeta,
                                    TransientSpecs *specs);

/*----------------------------------------------------------------------
 * L6 — Canonical Problem Solvers
 *----------------------------------------------------------------------*/

/**
 * @brief  Solve second-order step response analytically at time t
 * @param  omega_n      Natural frequency ωₙ
 * @param  zeta         Damping ratio ζ
 * @param  t            Time point t
 * @return y(t); NaN if invalid
 *
 * L6: Analytical solution for y(t) of G(s) = ωₙ²/(s² + 2ζωₙs + ωₙ²)
 *     driven by unit step:
 *
 *     Underdamped (0 < ζ < 1):
 *       y(t) = 1 - (e^{-ζωₙt}/√(1-ζ²))·sin(ω_d·t + φ)
 *       where ω_d = ωₙ√(1-ζ²), φ = acos(ζ)
 *
 *     Critically damped (ζ = 1):
 *       y(t) = 1 - e^{-ωₙt}·(1 + ωₙt)
 *
 *     Overdamped (ζ > 1):
 *       y(t) = 1 + (1/(s₂-s₁))·(s₂·e^{s₁t} - s₁·e^{s₂t})
 *       where s₁,₂ = ωₙ·(-ζ ± √(ζ²-1))
 *
 *     Undamped (ζ = 0):
 *       y(t) = 1 - cos(ωₙt)
 */
double second_order_step_response_at(double omega_n, double zeta, double t);

/**
 * @brief  Compute impulse response h(t) of a rational transfer function
 *         using residue method (partial fraction decomposition)
 * @param  tf           Transfer function (strictly proper)
 * @param  T_end        End time
 * @param  n_samples    Number of samples
 * @param  response     Output impulse response
 * @return 0 on success
 *
 * L6: h(t) = Σ A_k·e^{p_k·t} for distinct poles, with residue A_k.
 *     Uses partial fraction decomposition + analytic evaluation.
 *     More accurate than numerical inverse LT for rational functions.
 */
int impulse_response_analytic(const TransferFunction *tf,
                               double T_end, int n_samples,
                               TimeResponse *response);

/**
 * @brief  Compute steady-state error for unity feedback system
 * @param  G            Plant transfer function G(s)
 * @param  input_type   RESPONSE_STEP or RESPONSE_RAMP or RESPONSE_PARABOLA
 * @param  ess          Output: steady-state error
 * @return 0 on success
 *
 * L6: For unity feedback with H(s)=1:
 *     - Step (K_p): e_ss = 1/(1+K_p), K_p = lim_{s→0} G(s)
 *     - Ramp (K_v): e_ss = 1/K_v,   K_v = lim_{s→0} s·G(s)
 *     - Parabola (K_a): e_ss = 1/K_a, K_a = lim_{s→0} s²·G(s)
 *     Complexity: O(1).
 */
int steady_state_error(const TransferFunction *G,
                        ResponseType input_type,
                        double *ess);

/**
 * @brief  Free memory allocated by compute_time_response or impulse_response_analytic
 * @param  response   Time response to free
 */
void time_response_free(TimeResponse *response);

#ifdef __cplusplus
}
#endif

#endif /* SYSTEM_RESPONSE_H */
