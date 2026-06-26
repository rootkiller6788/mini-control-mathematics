/**
 * @file transfer_function.h
 * @brief Transfer Function Analysis — Poles, Zeros, Gain, Frequency Response
 *
 * Knowledge Coverage: L1 (Definitions), L2 (Core Concepts), L5 (Methods)
 *
 * A transfer function G(s) = N(s)/D(s) characterizes an LTI system's
 * input-output relationship in the complex frequency domain.
 *
 * Key analyses:
 *   - Pole-zero map: stability assessment
 *   - DC gain: G(0) = steady-state gain for step input
 *   - Type number: number of integrators (poles at s=0)
 *   - Relative degree: excess poles over zeros
 *   - Natural frequency ωₙ and damping ratio ζ for 2nd-order systems
 *   - Frequency response: G(jω) = |G(jω)|·e^{j·∠G(jω)}
 *   - Bandwidth, resonant peak, resonant frequency
 *
 * Reference:
 *   - Franklin, Powell, Emami-Naeini, "Feedback Control of Dynamic Systems"
 *   - Dorf & Bishop, "Modern Control Systems"
 */

#ifndef TRANSFER_FUNCTION_H
#define TRANSFER_FUNCTION_H

#include <stddef.h>
#include <complex.h>
#include "laplace_core.h"
#include "z_transform.h"

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------
 * L1 — Transfer Function Type Definitions
 *----------------------------------------------------------------------*/

/**
 * @struct TransferFunction
 * @brief  Unified transfer function representation (works in s or z domain)
 *
 * L1: This is the central data structure of classical control theory.
 *     G(s) = K · Π(z - zᵢ) / Π(s - pⱼ)
 *
 * Domain is indicated by sampling_period:
 *   - T = 0.0 → continuous-time (Laplace/s-domain)
 *   - T > 0.0 → discrete-time   (Z-transform/z-domain)
 */
typedef struct {
    int      num_order;      /**< Degree of numerator N(·)           */
    int      den_order;      /**< Degree of denominator D(·)         */
    double   num_coeff[LAPLACE_MAX_ORDER + 1];  /**< N coefficients  */
    double   den_coeff[LAPLACE_MAX_ORDER + 1];  /**< D coefficients  */
    double   gain;           /**< Overall gain K                     */
    double   sampling_period; /**< 0 = continuous, >0 = discrete     */
} TransferFunction;

/**
 * @enum SystemType
 * @brief  Classification by number of free integrators (poles at origin)
 *
 * L1/L2: Type-N systems can track polynomial references of order N-1
 *        with zero steady-state error.
 */
typedef enum {
    SYSTEM_TYPE_0  = 0,   /**< No free integrator → finite step error    */
    SYSTEM_TYPE_1  = 1,   /**< One integrator → zero step error          */
    SYSTEM_TYPE_2  = 2,   /**< Two integrators → zero ramp error         */
    SYSTEM_TYPE_3  = 3    /**< Three integrators → zero parabolic error  */
} SystemType;

/**
 * @struct FrequencyResponse
 * @brief  Bode plot data at a single frequency ω
 *
 * L2/L5: Frequency-domain characterization of LTI systems.
 *        Magnitude in dB, phase in degrees.
 */
typedef struct {
    double   omega;         /**< Frequency ω in rad/s                */
    double   magnitude_db;  /**< 20·log₁₀|G(jω)| in dB              */
    double   phase_deg;     /**< ∠G(jω) in degrees                   */
    double   real_part;     /**< Re{G(jω)}                           */
    double   imag_part;     /**< Im{G(jω)}                           */
} FrequencyResponse;

/**
 * @struct SecondOrderParams
 * @brief  Standard 2nd-order system parameters
 *
 * L2: For G(s) = ωₙ²/(s² + 2ζωₙs + ωₙ²):
 *   - ωₙ: undamped natural frequency
 *   - ζ:  damping ratio (0=undamped, <1=underdamped, =1=critically, >1=overdamped)
 */
typedef struct {
    double   omega_n;       /**< Natural frequency ωₙ (rad/s)        */
    double   zeta;          /**< Damping ratio ζ                      */
    int      is_standard_form; /**< 1 if matches standard 2nd-order  */
} SecondOrderParams;

/*----------------------------------------------------------------------
 * L2 — Transfer Function Analysis
 *----------------------------------------------------------------------*/

/**
 * @brief  Evaluate transfer function G(s) or G(z) at a complex point
 * @param  tf   Transfer function
 * @param  s_z  Complex evaluation point (s or z)
 * @return G(s/z); NaN at poles
 */
double complex tf_evaluate(const TransferFunction *tf, double complex s_z);

/**
 * @brief  Compute DC gain: G(0) for continuous or G(1) for discrete
 * @param  tf   Transfer function
 * @return DC gain; NaN if denominator at 0 (or 1) is zero
 *
 * L2: DC gain determines steady-state output for a unit step input
 *     (assuming stable system).
 */
double tf_dc_gain(const TransferFunction *tf);

/**
 * @brief  Determine system type (number of poles at origin for CT, at z=1 for DT)
 * @param  tf   Transfer function
 * @return SystemType
 */
SystemType tf_system_type(const TransferFunction *tf);

/**
 * @brief  Extract 2nd-order parameters ωₙ and ζ
 * @param  tf     Transfer function (must have den_order >= 2)
 * @param  params Output parameters
 * @return 0 on success, -1 if not second-order
 *
 * L5: For a standard second-order denominator s² + 2ζωₙs + ωₙ²:
 *        ζ = a₁/(2·√a₀), ωₙ = √a₀
 *     where D(s) = a₀ + a₁s + a₂s² (normalized so a₂ = 1).
 */
int tf_extract_second_order(const TransferFunction *tf,
                             SecondOrderParams *params);

/**
 * @brief  Compute frequency response G(jω) at a given ω
 * @param  tf   Continuous transfer function
 * @param  omega Frequency in rad/s
 * @param  resp Output frequency response (caller-allocated)
 * @return 0 on success
 *
 * L5: G(jω) = N(jω)/D(jω). Magnitude = |G(jω)|, phase = atan2(Im, Re).
 *     Complexity: O(n+m).
 */
int tf_frequency_response(const TransferFunction *tf,
                           double omega,
                           FrequencyResponse *resp);

/**
 * @brief  Compute frequency response over a range (Bode data array)
 * @param  tf       Transfer function
 * @param  omega_start, omega_end, n_points  Frequency sweep range
 * @param  responses Output array (caller-allocated, size n_points)
 * @return 0 on success
 *
 * L5: Bode plot data generation.
 *     Logarithmic spacing recommended for ω.
 *     Complexity: O(N·(n+m)) for N points.
 */
int tf_frequency_sweep(const TransferFunction *tf,
                        double omega_start,
                        double omega_end,
                        int n_points,
                        FrequencyResponse *responses);

/**
 * @brief  Find resonant peak frequency ω_r and peak magnitude M_r
 * @param  tf        Transfer function
 * @param  omega_r   Output: resonant frequency (rad/s)
 * @param  M_r_db    Output: peak magnitude (dB)
 * @return 0 on success, -1 if no resonance
 *
 * L5: For 2nd-order underdamped (ζ<0.707): ω_r = ωₙ·√(1-2ζ²)
 *      For higher-order systems: numerical search on |G(jω)|.
 */
int tf_resonant_peak(const TransferFunction *tf,
                      double *omega_r,
                      double *M_r_db);

/**
 * @brief  Compute bandwidth ω_bw (frequency where magnitude drops 3dB from DC)
 * @param  tf        Transfer function
 * @param  omega_bw  Output: bandwidth frequency (rad/s)
 * @return 0 on success, -1 if cannot determine
 *
 * L5: Bandwidth is a key specification for control design.
 *     |G(jω_bw)| = |G(0)|/√2 → -3dB from DC.
 */
int tf_bandwidth(const TransferFunction *tf,
                  double *omega_bw);

/**
 * @brief  Poles and zeros extraction for TransferFunction
 * @param  tf     Transfer function
 * @param  pz     Output pole-zero structure
 * @return 0 on success
 */
int tf_pole_zero_map(const TransferFunction *tf,
                      LaplacePoleZero *pz);

/**
 * @brief  Check if transfer function is minimum-phase
 *         (all zeros in LHP for CT, inside unit circle for DT)
 * @param  tf   Transfer function
 * @return 1 if minimum-phase, 0 otherwise
 *
 * L8: Minimum-phase systems have no RHP zeros and allow certain
 *     control design simplifications.
 */
int tf_is_minimum_phase(const TransferFunction *tf);

/**
 * @brief  Compute relative degree (excess poles over zeros)
 * @param  tf   Transfer function
 * @return den_order - num_order
 *
 * L2: Relative degree determines high-frequency roll-off rate:
 *     -40dB/decade per relative degree.
 */
int tf_relative_degree(const TransferFunction *tf);

#ifdef __cplusplus
}
#endif

#endif /* TRANSFER_FUNCTION_H */
