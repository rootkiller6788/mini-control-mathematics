/**
 * @file z_transform.h
 * @brief Z-Transform — Discrete-Time Signal Analysis
 *
 * Knowledge Coverage: L1 (Definitions), L2 (Core Concepts), L4 (Fundamental Laws)
 *
 * The Z-transform maps a discrete-time sequence x[n] to the complex z-plane:
 *   X(z) = Z{x[n]} = Σ_{n=0}^∞ x[n] · z^(-n)
 *
 * Key properties:
 *   - Linearity:       Z{a·x[n] + b·y[n]} = a·X(z) + b·Y(z)
 *   - Time shift:      Z{x[n - k]} = z^(-k)·X(z)   (causal)
 *   - Time advance:    Z{x[n + k]} = z^k·X(z) - Σ_{m=0}^{k-1} x[m]·z^{k-m}
 *   - Scaling:         Z{a^n·x[n]} = X(z/a)
 *   - Convolution:     Z{x[n] * y[n]} = X(z)·Y(z)
 *   - Initial Value:   x[0] = lim_{z→∞} X(z)
 *   - Final Value:     x[∞] = lim_{z→1} (z-1)·X(z)  (if poles in |z|<1)
 *
 * Reference:
 *   - Oppenheim & Schafer, "Discrete-Time Signal Processing" (3rd ed.)
 *   - Ogata, "Discrete-Time Control Systems" (2nd ed.)
 */

#ifndef Z_TRANSFORM_H
#define Z_TRANSFORM_H

#include <stddef.h>
#include <complex.h>
#include "laplace_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------
 * L1 — Core Type Definitions
 *----------------------------------------------------------------------*/

#define Z_MAX_ORDER  64

/**
 * @struct ZPoly
 * @brief  Polynomial in z^(-1): P(z) = b₀ + b₁·z^(-1) + ... + bₙ·z^(-n)
 *
 * L1: Discrete transfer functions use negative powers of z.
 */
typedef struct {
    int    order;                      /**< Degree n                     */
    double coeff[Z_MAX_ORDER + 1];     /**< Coefficients b₀..bₙ          */
} ZPoly;

/**
 * @struct ZRational
 * @brief  Rational function in z: G(z) = N(z)/D(z) = Σbₖz^(-k) / Σaₖz^(-k)
 *
 * L1: Every discrete-time LTI system has a rational Z-transform.
 */
typedef struct {
    ZPoly numerator;        /**< N(z) numerator           */
    ZPoly denominator;      /**< D(z) denominator         */
    double sampling_period; /**< T in seconds              */
} ZRational;

/**
 * @enum ZROCType
 * @brief Region of Convergence type for Z-transform
 *
 * L1: ROC shape determines uniqueness and causality of the sequence.
 */
typedef enum {
    Z_ROC_EXTERIOR   = 0,  /**< |z| > r (causal sequence)            */
    Z_ROC_INTERIOR   = 1,  /**< |z| < r (anti-causal)                */
    Z_ROC_ANNULUS    = 2,  /**< r₁ < |z| < r₂ (two-sided)            */
    Z_ROC_ENTIRE     = 3,  /**< Entire z-plane except possibly z=0,∞  */
    Z_ROC_EMPTY      = 4   /**< No ROC (no Z-transform exists)        */
} ZROCType;

/**
 * @struct ZPoleZero
 * @brief  Pole-zero description in the z-plane
 */
typedef struct {
    int            num_poles;
    int            num_zeros;
    double complex poles[Z_MAX_ORDER];
    double complex zeros[Z_MAX_ORDER];
    double         gain;
} ZPoleZero;

/**
 * @struct DiscreteSignal
 * @brief  Finite-length discrete-time sequence x[n]
 */
typedef struct {
    int      length;         /**< Number of samples        */
    int      start_index;    /**< Index of x[0]            */
    double  *values;         /**< Array x[0..length-1]     */
} DiscreteSignal;

/*----------------------------------------------------------------------
 * L2 — Core Property Functions
 *----------------------------------------------------------------------*/

/**
 * @brief  Evaluate Z-poly at complex z: P(z) = Σ pₖ·z^(-k)
 * @param  poly  Polynomial in z^(-1)
 * @param  z     Complex evaluation point
 * @return P(z)
 *
 * L2: Polynomial evaluation using Horner's method in z^(-1).
 *     Complexity: O(n).
 */
double complex z_poly_eval(const ZPoly *poly, double complex z);

/**
 * @brief  Evaluate Z-rational at z: G(z) = N(z)/D(z)
 * @param  rat   Rational function
 * @param  z     Complex point
 * @return G(z); NaN if D(z)=0 (pole)
 */
double complex z_rational_eval(const ZRational *rat, double complex z);

/**
 * @brief  Verify linearity numerically
 */
int z_verify_linearity(double a, double b,
                        double complex Xz, double complex Yz,
                        double complex combined,
                        double tolerance);

/**
 * @brief  Time-shift property: multiply by z^(-k) for delay of k samples
 * @return z^(-k)
 *
 * L2: Pure delay of k samples ↔ z^(-k) factor.
 *     Used extensively in digital control loops and FIR filter design.
 */
double complex z_delay_factor(double complex z, int k);

/**
 * @brief  Time-advance property: multiply by z^k subtracting initial conditions
 * @param  X        X(z) before advance
 * @param  k        Advance amount
 * @param  initials Initial values x[0]..x[k-1]
 * @param  result   Output: Z{x[n+k]}
 * @return 0 on success
 *
 * L2: Time advancing requires initial condition terms.
 *     Z{x[n+k]} = z^k·X(z) - Σ_{m=0}^{k-1} x[m]·z^{k-m}
 */
int z_time_advance(const ZRational *X, int k,
                    const double *initials,
                    ZRational *result);

/**
 * @brief  Frequency scaling: Z{a^n·x[n]} = X(z/a)
 * @param  X        Original transform
 * @param  a        Scaling factor (real, a ≠ 0)
 * @param  result   Output: X(z/a)
 * @return 0 on success
 *
 * L2: Exponential modulation in discrete time.
 */
int z_frequency_scaling(const ZRational *X, double a,
                         ZRational *result);

/**
 * @brief  Convolution property: Z{x[n] * y[n]} = X(z)·Y(z)
 * @param  X,Y      Z-transforms
 * @param  result   X(z)·Y(z)
 * @return 0 on success
 *
 * L4: Discrete convolution ↔ multiplication in z-domain.
 *     This is the basis of FIR and IIR filter analysis.
 */
int z_convolution(const ZRational *X, const ZRational *Y,
                   ZRational *result);

/**
 * @brief  Initial value theorem: x[0] = lim_{z→∞} X(z)
 * @param  X     Z-transform X(z)
 * @return x[0]; NaN if limit diverges
 *
 * L4: Requires X(z) to be a proper rational function.
 */
double z_initial_value(const ZRational *X);

/**
 * @brief  Final value theorem: x[∞] = lim_{z→1} (z-1)·X(z)
 * @param  X     Z-transform X(z)
 * @return x[∞]; NaN if any pole on or outside unit circle
 *
 * L4: Only valid when all poles are inside the unit circle.
 *     Complexity: O(n) — checks stability first, then computes limit.
 */
double z_final_value(const ZRational *X);

/*----------------------------------------------------------------------
 * L5 — Inverse Z-Transform Methods
 *----------------------------------------------------------------------*/

/**
 * @brief  Inverse Z-transform by power series expansion (long division)
 * @param  X        Rational Z-transform X(z)
 * @param  n_terms  Number of sequence terms to compute
 * @param  x        Output array x[0..n_terms-1] (caller allocated)
 * @return 0 on success
 *
 * L5: Long division of numerator by denominator yields x[n] directly.
 *     Complexity: O(N·d) where d = denominator order.
 *     Suitable for first few terms only (numerical error accumulates).
 */
int z_inverse_power_series(const ZRational *X,
                            int n_terms,
                            double *x);

/**
 * @brief  Inverse Z-transform by partial fraction expansion + table lookup
 * @param  X        Rational Z-transform
 * @param  n_terms  Output sequence length
 * @param  x        Output array x[0..n_terms-1]
 * @return 0 on success, -1 if denominator cannot be factored
 *
 * L5: Decompose into first-order terms, then use known inverse pairs.
 *     More stable than power series for moderate to large n.
 *     Complexity: O(N·p) where p = number of poles.
 */
int z_inverse_partial_fraction(const ZRational *X,
                                int n_terms,
                                double *x);

/**
 * @brief  Inverse Z-transform using residue method (Cauchy integral)
 * @param  X        Rational Z-transform
 * @param  n        Time index n
 * @return x[n]     via sum of residues inside ROC
 *
 * L4/L5: x[n] = (1/2πj)∮ X(z)·z^(n-1) dz = Σ Res{X(z)·z^(n-1)}
 *        Evaluated at poles inside the ROC.
 *        Complexity: O(M²) for M poles.
 */
double z_inverse_residue(const ZRational *X, int n);

/*----------------------------------------------------------------------
 * L3 — Polynomial Operations for Z-domain
 *----------------------------------------------------------------------*/

int z_poly_multiply(const ZPoly *P, const ZPoly *Q, ZPoly *R);
int z_poly_add(const ZPoly *P, const ZPoly *Q, ZPoly *R);

/**
 * @brief  Map Z-domain poles back to s-domain via s = ln(z)/T
 * @param  zp       Z-plane pole-zero description
 * @param  T        Sampling period
 * @param  sp       Output: equivalent s-plane poles/zeros
 * @return 0 on success
 *
 * L3/L5: inverse mapping z → s, used for continuous-time interpretation
 *        of discrete-time designs.
 */
int z_poles_to_s_plane(const ZPoleZero *zp, double T,
                        LaplacePoleZero *sp);

#ifdef __cplusplus
}
#endif

#endif /* Z_TRANSFORM_H */
