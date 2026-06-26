/**
 * @file bilinear_transform.h
 * @brief Bilinear (Tustin) Transform — Continuous ↔ Discrete Mapping
 *
 * Knowledge Coverage: L5 (Algorithms), L3 (Mathematical Structures)
 *
 * The bilinear transform maps the continuous s-domain to the discrete z-domain:
 *
 *           2    1 - z⁻¹
 *   s = ───── · ─────────     (Tustin approximation)
 *          T    1 + z⁻¹
 *
 * Inverse mapping (z → s):
 *
 *        2 + sT
 *   z = ─────────
 *        2 - sT
 *
 * Key properties:
 *   - Maps the entire LHP (stable region) into the unit circle → stability preserved
 *   - Maps jω-axis to unit circle → frequency response preserved (with warping)
 *   - One-to-one mapping → no aliasing
 *   - Frequency warping: ω_a = (2/T)·tan(ω_d·T/2)
 *
 * Discretization methods compared:
 *   1. Bilinear (Tustin) — stability preserving, frequency warping
 *   2. Impulse Invariance — preserves impulse response, aliasing possible
 *   3. Zero-Order Hold (ZOH) — exact for step inputs
 *   4. Forward Euler — simple, may not preserve stability
 *   5. Backward Euler — preserves stability, less accurate
 *
 * References:
 *   - Tustin (1947), "A method of analysing the behaviour of linear systems"
 *   - Oppenheim & Schafer, "Discrete-Time Signal Processing"
 *   - Åström & Wittenmark, "Computer-Controlled Systems"
 */

#ifndef BILINEAR_TRANSFORM_H
#define BILINEAR_TRANSFORM_H

#include <stddef.h>
#include "laplace_core.h"
#include "z_transform.h"
#include "transfer_function.h"

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------
 * L3 — Mapping Parameters and Types
 *----------------------------------------------------------------------*/

/**
 * @enum DiscretizationMethod
 * @brief  Methods for converting continuous to discrete transfer functions
 */
typedef enum {
    DISC_METHOD_BILINEAR        = 0, /**< Tustin (bilinear transform)     */
    DISC_METHOD_IMPULSE_INV     = 1, /**< Impulse invariance               */
    DISC_METHOD_ZOH             = 2, /**< Zero-Order Hold equivalent       */
    DISC_METHOD_FORWARD_EULER   = 3, /**< Forward Euler (s → (z-1)/T)     */
    DISC_METHOD_BACKWARD_EULER  = 4, /**< Backward Euler (s → (1-z⁻¹)/T)  */
    DISC_METHOD_MATCHED_Z       = 5  /**< Matched pole-zero mapping        */
} DiscretizationMethod;

/*----------------------------------------------------------------------
 * L5 — Bilinear (Tustin) Transform
 *----------------------------------------------------------------------*/

/**
 * @brief  Convert continuous transfer function G(s) to discrete G(z)
 *         using the bilinear (Tustin) transform.
 * @param  ct_tf   Continuous-time transfer function G(s)
 * @param  T       Sampling period (seconds)
 * @param  prewarp_freq  Prewarp frequency ω_p (rad/s); 0 = no prewarping
 * @param  dt_tf   Output: discrete-time transfer function G(z)
 * @return 0 on success, -1 on overflow
 *
 * L5: The Tustin substitution s = (2/T)·(z-1)/(z+1) into G(s) yields G(z).
 *     Each s^k term expands via binomial theorem, requiring careful
 *     polynomial manipulation.
 *     Complexity: O(n²) where n = max(numerator_order, denominator_order).
 *
 *     Theorem: The bilinear transform preserves stability — LHP maps to
 *     unit circle interior.
 */
int bilinear_discretize(const TransferFunction *ct_tf,
                         double T, double prewarp_freq,
                         TransferFunction *dt_tf);

/**
 * @brief  Compute frequency prewarping: ω_p = (2/T)·tan(ω_d·T/2)
 * @param  omega_desired  Desired analog frequency ω_d (rad/s)
 * @param  T              Sampling period
 * @return Prewarped frequency ω_p
 *
 * L5: Without prewarping, the bilinear transform causes nonlinear
 *     frequency distortion. Prewarping corrects one critical frequency.
 *     Used when exact matching at a specific ω is required (e.g., notch filter).
 */
double bilinear_prewarp_frequency(double omega_desired, double T);

/**
 * @brief  Warp a digital frequency back to its equivalent analog frequency
 * @param  omega_digital  Digital frequency ω_d (rad/sample → rad/s using T)
 * @param  T              Sampling period
 * @return Equivalent analog frequency ω_a = (2/T)·tan(ω_d·T/2)
 *
 * L5: Inverse of the prewarp relationship. Used to interpret digital
 *     filter characteristics in analog terms.
 */
double bilinear_analog_frequency(double omega_digital, double T);

/*----------------------------------------------------------------------
 * L5 — Alternative Discretization Methods
 *----------------------------------------------------------------------*/

/**
 * @brief  Discretize using impulse invariance (h[n] = T·h_c(nT))
 * @param  ct_tf   Continuous G(s)
 * @param  T       Sampling period
 * @param  dt_tf   Output discrete G(z)
 * @return 0 on success
 *
 * L5: Impulse invariance preserves the impulse response exactly at
 *     sampling instants. Requires partial fraction decomposition.
 *     G(z) = Σ T·A_k/(1 - e^{p_k·T}·z⁻¹) for poles p_k with residues A_k.
 *     Limitation: aliasing if G(s) is not bandlimited.
 */
int impulse_invariance_discretize(const TransferFunction *ct_tf,
                                   double T, TransferFunction *dt_tf);

/**
 * @brief  Discretize using Zero-Order Hold (ZOH) equivalent
 * @param  ct_tf   Continuous G(s)
 * @param  T       Sampling period
 * @param  dt_tf   Output: G(z) = (1-z⁻¹)·Z{G(s)/s}
 * @return 0 on success
 *
 * L5: ZOH models the DAC output: holds each sample constant for T seconds.
 *     G_zoh(z) = (z-1)/z · Z{G(s)/s}
 *     Most commonly used in digital control system design.
 */
int zoh_discretize(const TransferFunction *ct_tf,
                    double T, TransferFunction *dt_tf);

/**
 * @brief  Forward Euler discretization: s ≈ (z - 1) / T
 * @param  ct_tf   Continuous G(s)
 * @param  T       Sampling period
 * @param  dt_tf   Output discrete G(z)
 * @return 0 on success
 *
 * L5: Forward Euler is the simplest mapping. Substitutes s → (z-1)/T.
 *     Warning: may map stable CT poles to unstable DT poles for large T.
 *     Stability requires T < 2/|p_max| for each pole.
 */
int forward_euler_discretize(const TransferFunction *ct_tf,
                              double T, TransferFunction *dt_tf);

/**
 * @brief  Backward Euler discretization: s ≈ (1 - z⁻¹) / T
 * @param  ct_tf   Continuous G(s)
 * @param  T       Sampling period
 * @param  dt_tf   Output discrete G(z)
 * @return 0 on success
 *
 * L5: Backward Euler substitutes s → (z-1)/(Tz).
 *     Preserves stability (maps LHP into a circle of radius 1/2 centered at 0.5).
 *     More accurate than forward Euler for stiff systems.
 */
int backward_euler_discretize(const TransferFunction *ct_tf,
                               double T, TransferFunction *dt_tf);

/**
 * @brief  Matched pole-zero mapping
 * @param  ct_tf   Continuous G(s)
 * @param  T       Sampling period
 * @param  dt_tf   Output discrete G(z)
 * @return 0 on success
 *
 * L5: Maps each pole/zero via z = e^{sT}. Zeros at s=∞ map to z=-1.
 *     Gain adjusted so DC gains match.
 */
int matched_z_discretize(const TransferFunction *ct_tf,
                          double T, TransferFunction *dt_tf);

/*----------------------------------------------------------------------
 * L5 — Continuous ↔ Discrete Conversion (Unified API)
 *----------------------------------------------------------------------*/

/**
 * @brief  Discretize using specified method (unified interface)
 * @param  ct_tf   Continuous transfer function
 * @param  T       Sampling period
 * @param  method  Discretization method
 * @param  prewarp_freq  Prewarp frequency (only for bilinear)
 * @param  dt_tf   Output discrete transfer function
 * @return 0 on success
 */
int discretize(const TransferFunction *ct_tf,
                double T, DiscretizationMethod method,
                double prewarp_freq, TransferFunction *dt_tf);

/**
 * @brief  Convert continuous transfer function back from discrete
 *         using inverse bilinear: z → (2+sT)/(2-sT)
 * @param  dt_tf   Discrete transfer function G(z)
 * @param  T       Original sampling period
 * @param  ct_tf   Output continuous G(s)
 * @return 0 on success
 *
 * L5: Inverse Tustin transform. Used to analyze discrete designs in
 *     the continuous domain.
 */
int inverse_bilinear(const TransferFunction *dt_tf,
                      double T, TransferFunction *ct_tf);

#ifdef __cplusplus
}
#endif

#endif /* BILINEAR_TRANSFORM_H */
