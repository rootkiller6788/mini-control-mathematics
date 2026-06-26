/**
 * @file digital_filter.h
 * @brief Digital Filter Design — Z-domain IIR and FIR Filters
 *
 * Knowledge Coverage: L5 (Algorithms), L6 (Canonical Problems), L7 (Applications)
 *
 * Digital filter design through Z-domain specifications:
 *
 *   IIR (Infinite Impulse Response):
 *     - Butterworth: maximally flat magnitude in passband
 *     - Chebyshev Type I: equiripple in passband, monotonic stopband
 *     - Chebyshev Type II: monotonic passband, equiripple stopband
 *     - Elliptic (Cauer): equiripple in both passband and stopband
 *     - Design via bilinear transform from analog prototype
 *
 *   FIR (Finite Impulse Response):
 *     - Window method (Hann, Hamming, Blackman, Kaiser)
 *     - Frequency sampling method
 *     - Parks-McClellan (Remez exchange) algorithm
 *
 * Filter structures:
 *   - Direct Form I, Direct Form II (canonical)
 *   - Cascade of second-order sections (biquads)
 *   - Parallel form
 *
 * References:
 *   - Parks & Burrus, "Digital Filter Design" (1987)
 *   - Proakis & Manolakis, "Digital Signal Processing" (4th ed.)
 *   - Rabiner & Gold, "Theory and Application of Digital Signal Processing"
 */

#ifndef DIGITAL_FILTER_H
#define DIGITAL_FILTER_H

#include <stddef.h>
#include "z_transform.h"
#include "transfer_function.h"

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------
 * L1 — Filter Specification Types
 *----------------------------------------------------------------------*/

/**
 * @enum FilterType
 * @brief  Standard filter types
 */
typedef enum {
    FILTER_LOWPASS   = 0,    /**< Pass low frequencies, attenuate high     */
    FILTER_HIGHPASS  = 1,    /**< Pass high frequencies, attenuate low     */
    FILTER_BANDPASS  = 2,    /**< Pass a band, attenuate outside           */
    FILTER_BANDSTOP  = 3     /**< Attenuate a band, pass outside (notch)   */
} FilterType;

/**
 * @enum IIRDesignMethod
 * @brief  Analog prototype → digital IIR design methods
 */
typedef enum {
    IIR_BUTTERWORTH    = 0,  /**< Maximally flat magnitude                */
    IIR_CHEBYSHEV_I    = 1,  /**< Equiripple passband                     */
    IIR_CHEBYSHEV_II   = 2,  /**< Equiripple stopband                     */
    IIR_ELLIPTIC       = 3   /**< Equiripple passband + stopband          */
} IIRDesignMethod;

/**
 * @enum FIRWindowType
 * @brief  Window functions for FIR filter design
 */
typedef enum {
    FIR_WIN_RECTANGULAR = 0, /**< w[n] = 1 (sharpest transition, high sidelobes) */
    FIR_WIN_HANN        = 1, /**< w[n] = 0.5(1 - cos(2πn/N))                   */
    FIR_WIN_HAMMING     = 2, /**< w[n] = 0.54 - 0.46·cos(2πn/N)               */
    FIR_WIN_BLACKMAN    = 3, /**< w[n] = 0.42 - 0.5·cos(2πn/N) + 0.08·cos(4πn/N) */
    FIR_WIN_KAISER      = 4  /**< Kaiser window with adjustable β              */
} FIRWindowType;

/**
 * @struct FilterSpec
 * @brief  Filter design specification (frequency-domain)
 */
typedef struct {
    FilterType  type;               /**< LP/HP/BP/BS                    */
    double      sampling_freq;      /**< Fs in Hz                       */
    double      passband_edge;      /**< Passband edge freq (Hz)        */
    double      stopband_edge;      /**< Stopband edge freq (Hz)        */
    double      passband_ripple_db; /**< Max passband ripple (dB)       */
    double      stopband_atten_db;  /**< Min stopband attenuation (dB)  */
    /* For bandpass/bandstop, these are the second edges */
    double      passband_edge2;     /**< Second passband edge (BP/BS)   */
    double      stopband_edge2;     /**< Second stopband edge (BP/BS)   */
} FilterSpec;

/**
 * @struct FilterCoeff
 * @brief  Realized digital filter coefficients
 */
typedef struct {
    double  *b;             /**< Feed-forward coefficients b₀..b_M     */
    double  *a;             /**< Feed-back coefficients a₀..a_N (a₀=1) */
    int     M;              /**< Numerator order (number of zeros)      */
    int     N;              /**< Denominator order (number of poles)   */
    double  sampling_freq;  /**< Sampling frequency (Hz)               */
    int     is_fir;         /**< 1 if FIR (all a_i = 0 except a₀=1)   */
} FilterCoeff;

/**
 * @struct BiquadSection
 * @brief  Second-order section (biquad) for cascade filter implementation
 *
 * H(z) = (b₀ + b₁·z⁻¹ + b₂·z⁻²) / (1 + a₁·z⁻¹ + a₂·z⁻²)
 */
typedef struct {
    double   b0, b1, b2;       /**< Numerator coefficients             */
    double   a1, a2;           /**< Denominator coefficients (a₀ = 1)  */
    double   gain;              /**< Section gain                       */
} BiquadSection;

/*----------------------------------------------------------------------
 * L5 — IIR Filter Design
 *----------------------------------------------------------------------*/

/**
 * @brief  Design IIR digital filter from frequency specifications
 * @param  spec   Filter specification
 * @param  method Design method (Butterworth, Chebyshev I/II, Elliptic)
 * @param  coeff  Output filter coefficients (caller frees b and a)
 * @return 0 on success, -1 if specification cannot be met
 *
 * L5: Complete IIR design flow:
 *     1. Prewarp critical frequencies (for bilinear)
 *     2. Compute analog filter order from attenuation requirements
 *     3. Compute analog prototype poles/zeros
 *     4. Apply bilinear transform to obtain digital coefficients
 *     Complexity: O(N²) where N = filter order.
 */
int iir_filter_design(const FilterSpec *spec,
                       IIRDesignMethod method,
                       FilterCoeff *coeff);

/**
 * @brief  Compute Butterworth filter order required to meet spec
 * @param  spec   Filter specification
 * @return Required order N ≥ 1; -1 on invalid spec
 *
 * L5: N ≥ log₁₀((10^{0.1·A_s} - 1)/(10^{0.1·A_p} - 1)) / (2·log₁₀(ω_s/ω_p))
 *     where A_s = stopband attenuation dB, A_p = passband ripple dB.
 *     Complexity: O(1).
 */
int butterworth_order(const FilterSpec *spec);

/**
 * @brief  Compute Chebyshev Type I filter order
 * @param  spec   Filter specification
 * @return Required order N
 *
 * L5: N ≥ acosh(√((10^{0.1·A_s} - 1)/(10^{0.1·A_p} - 1))) / acosh(ω_s/ω_p)
 */
int chebyshev1_order(const FilterSpec *spec);

/**
 * @brief  Decompose IIR filter into cascade of biquad sections
 * @param  coeff       Direct-form filter coefficients
 * @param  sections    Output biquad array (caller allocated, size ceil(N/2))
 * @param  n_sections  Output: number of biquad sections
 * @return 0 on success
 *
 * L5: Biquad decomposition improves numerical stability for high-order IIR.
 *     Pairs complex-conjugate poles into real second-order sections.
 */
int iir_to_biquads(const FilterCoeff *coeff,
                    BiquadSection *sections,
                    int *n_sections);

/*----------------------------------------------------------------------
 * L5 — FIR Filter Design
 *----------------------------------------------------------------------*/

/**
 * @brief  Design FIR filter using window method (impulse response truncation)
 * @param  spec       Filter specification
 * @param  win_type   Window function type
 * @param  beta       Kaiser window β parameter (only for KAISER, ignored otherwise)
 * @param  coeff      Output filter coefficients (b array; a = {1})
 * @return 0 on success
 *
 * L5: Ideal impulse response h_d[n] = F⁻¹{H_d(e^{jω})}, then
 *     h[n] = h_d[n]·w[n] where w[n] is the chosen window.
 *     Filter order N estimated from transition width and stopband attenuation.
 *     Complexity: O(N).
 */
int fir_window_design(const FilterSpec *spec,
                       FIRWindowType win_type, double beta,
                       FilterCoeff *coeff);

/**
 * @brief  Estimate required FIR filter order for given spec
 * @param  spec   Filter specification
 * @return Estimated order N
 *
 * L5: Kaiser's formula: N ≈ (A_s - 7.95) / (14.36 · Δf/Fs)
 *     where Δf = transition width, A_s = stopband attenuation dB.
 */
int fir_estimate_order(const FilterSpec *spec);

/**
 * @brief  Compute window coefficient for given window type
 * @param  n          Sample index (0..N)
 * @param  N          Filter order
 * @param  win_type   Window type
 * @param  beta       Kaiser β parameter
 * @return Window coefficient w[n]
 *
 * L5: Window function evaluation.
 */
double fir_window_coeff(int n, int N, FIRWindowType win_type, double beta);

/*----------------------------------------------------------------------
 * L5 — Filter Analysis
 *----------------------------------------------------------------------*/

/**
 * @brief  Compute magnitude and phase response of a digital filter
 * @param  coeff   Filter coefficients
 * @param  omega   Normalized frequency (rad/sample, 0 to π)
 * @param  mag_db  Output: magnitude in dB
 * @param  phase   Output: phase in radians
 * @return 0 on success
 */
int filter_freq_response(const FilterCoeff *coeff, double omega,
                          double *mag_db, double *phase);

/**
 * @brief  Apply digital filter to input sequence
 * @param  coeff   Filter coefficients
 * @param  input   Input sequence x[0..length-1]
 * @param  length  Sequence length
 * @param  output  Output sequence y[0..length-1] (caller allocated)
 * @return 0 on success
 *
 * L6: Direct Form II transposed implementation.
 *     y[n] = Σ b_k·x[n-k] - Σ a_k·y[n-k]
 *     Complexity: O(N·(M+Nmax)) per sample.
 */
int filter_apply(const FilterCoeff *coeff,
                  const double *input, int length,
                  double *output);

/**
 * @brief  Apply cascade of biquad sections (more numerically stable)
 * @param  sections    Array of biquad sections
 * @param  n_sections  Number of sections
 * @param  input       Input sequence
 * @param  length      Sequence length
 * @param  output      Output sequence
 * @return 0 on success
 */
int biquad_filter_apply(const BiquadSection *sections, int n_sections,
                         const double *input, int length,
                         double *output);

/**
 * @brief  Free coefficient arrays allocated by design functions
 * @param  coeff   Filter coefficients to free
 */
void filter_coeff_free(FilterCoeff *coeff);

#ifdef __cplusplus
}
#endif

#endif /* DIGITAL_FILTER_H */
