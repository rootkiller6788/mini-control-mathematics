/**
 * @file transform_pairs.h
 * @brief Common Laplace and Z-Transform Pairs — Lookup and Generation
 *
 * Knowledge Coverage: L1 (Definitions), L3 (Mathematical Structures)
 *
 * Core transform pairs for control engineering:
 *
 *   Time Domain f(t) / x[n]         Laplace F(s)              Z-Transform X(z)
 *   ─────────────────────────────────────────────────────────────────────────
 *   δ(t) / δ[n]                      1                         1
 *   u(t) / u[n] (unit step)          1/s                       1/(1 - z⁻¹)
 *   t·u(t) / n·u[n]                  1/s²                      z⁻¹/(1 - z⁻¹)²
 *   t²·u(t) / n²·u[n]                2/s³                      z⁻¹(1+z⁻¹)/(1-z⁻¹)³
 *   e^(-at)·u(t) / aⁿ·u[n]           1/(s + a)                 1/(1 - a·z⁻¹)
 *   t·e^(-at)·u(t) / n·aⁿ·u[n]      1/(s + a)²                a·z⁻¹/(1 - a·z⁻¹)²
 *   sin(ωt)·u(t) / sin(ωn)·u[n]     ω/(s² + ω²)              z⁻¹sin(ω)/(1-2z⁻¹cos(ω)+z⁻²)
 *   cos(ωt)·u(t) / cos(ωn)·u[n]     s/(s² + ω²)              (1-z⁻¹cos(ω))/(1-2z⁻¹cos(ω)+z⁻²)
 *   e^(-at)sin(ωt) / aⁿsin(ωn)      ω/((s+a)² + ω²)          a·z⁻¹sin(ω)/(1-2a·z⁻¹cos(ω)+a²z⁻²)
 *   e^(-at)cos(ωt) / aⁿcos(ωn)      (s+a)/((s+a)² + ω²)     (1-a·z⁻¹cos(ω))/(1-2a·z⁻¹cos(ω)+a²z⁻²)
 *
 * References:
 *   - Churchill & Brown, "Complex Variables and Applications"
 *   - Jury, "Theory and Application of the Z-Transform Method"
 */

#ifndef TRANSFORM_PAIRS_H
#define TRANSFORM_PAIRS_H

#include <stddef.h>
#include "laplace_core.h"
#include "z_transform.h"

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------
 * L1 — Transform Pair Enumeration
 *----------------------------------------------------------------------*/

/**
 * @enum LaplacePairID
 * @brief  Enumerates common Laplace transform pairs
 */
typedef enum {
    LAPLACE_PAIR_IMPULSE = 0,       /**< δ(t)              ↔ 1              */
    LAPLACE_PAIR_STEP,              /**< u(t)              ↔ 1/s            */
    LAPLACE_PAIR_RAMP,              /**< t·u(t)            ↔ 1/s²           */
    LAPLACE_PAIR_PARABOLA,          /**< t²·u(t)           ↔ 2/s³           */
    LAPLACE_PAIR_EXP_DECAY,         /**< e^(-at)·u(t)      ↔ 1/(s+a)        */
    LAPLACE_PAIR_EXP_RAMP,          /**< t·e^(-at)·u(t)    ↔ 1/(s+a)²       */
    LAPLACE_PAIR_SINE,              /**< sin(ωt)·u(t)      ↔ ω/(s²+ω²)      */
    LAPLACE_PAIR_COSINE,            /**< cos(ωt)·u(t)      ↔ s/(s²+ω²)      */
    LAPLACE_PAIR_DAMPED_SINE,       /**< e^(-at)sin(ωt)    ↔ ω/((s+a)²+ω²)  */
    LAPLACE_PAIR_DAMPED_COSINE,     /**< e^(-at)cos(ωt)    ↔ (s+a)/((s+a)²+ω²) */
    LAPLACE_PAIR_SINH,              /**< sinh(at)·u(t)     ↔ a/(s²-a²)      */
    LAPLACE_PAIR_COSH,              /**< cosh(at)·u(t)     ↔ s/(s²-a²)      */
    LAPLACE_PAIR_COUNT
} LaplacePairID;

/**
 * @enum ZPairID
 * @brief  Enumerates common Z-transform pairs
 */
typedef enum {
    Z_PAIR_IMPULSE = 0,             /**< δ[n]              ↔ 1              */
    Z_PAIR_STEP,                    /**< u[n]              ↔ 1/(1-z⁻¹)      */
    Z_PAIR_RAMP,                    /**< n·u[n]            ↔ z⁻¹/(1-z⁻¹)²   */
    Z_PAIR_PARABOLA,                /**< n²·u[n]           ↔ z⁻¹(1+z⁻¹)/(1-z⁻¹)³ */
    Z_PAIR_EXP_DECAY,               /**< aⁿ·u[n]           ↔ 1/(1-a·z⁻¹)    */
    Z_PAIR_EXP_RAMP,                /**< n·aⁿ·u[n]         ↔ a·z⁻¹/(1-a·z⁻¹)² */
    Z_PAIR_SINE,                    /**< sin(ωn)·u[n]      ↔ z⁻¹sin(ω)/(..)  */
    Z_PAIR_COSINE,                  /**< cos(ωn)·u[n]      ↔ (1-z⁻¹cos(ω))/(..) */
    Z_PAIR_DAMPED_SINE,             /**< aⁿsin(ωn)·u[n]    ↔ a·z⁻¹sin(ω)/(..) */
    Z_PAIR_DAMPED_COSINE,           /**< aⁿcos(ωn)·u[n]    ↔ (1-a·z⁻¹cos(ω))/(..) */
    Z_PAIR_COUNT
} ZPairID;

/**
 * @struct LaplacePairDef
 * @brief  Single Laplace transform pair entry
 */
typedef struct {
    LaplacePairID   id;
    const char     *time_name;       /**< "δ(t)", "u(t)", ...          */
    const char     *s_domain_name;   /**< "1", "1/s", ...              */
    LaplaceRational transform;       /**< F(s) rational representation */
} LaplacePairDef;

/**
 * @struct ZPairDef
 * @brief  Single Z-transform pair entry
 */
typedef struct {
    ZPairID       id;
    const char   *time_name;       /**< "δ[n]", "u[n]", ...          */
    const char   *z_domain_name;   /**< "1", "1/(1-z⁻¹)", ...        */
    ZRational     transform;       /**< X(z) rational representation */
} ZPairDef;

/*----------------------------------------------------------------------
 * L3 — Pair Lookup and Custom Pair Creation
 *----------------------------------------------------------------------*/

/**
 * @brief  Look up a standard Laplace transform pair by ID
 * @param  id    Pair identifier
 * @param  param Parameter (a for exponential, ω for sinusoidal, etc.)
 * @param  pair  Output pair with pre-filled F(s)
 * @return 0 on success, -1 if param invalid
 */
int laplace_pair_get(LaplacePairID id, double param, LaplacePairDef *pair);

/**
 * @brief  Look up a standard Z-transform pair by ID
 * @param  id    Pair identifier
 * @param  param Parameter (a or ω, depending on pair)
 * @param  T     Sampling period for digital frequency mapping
 * @param  pair  Output pair
 * @return 0 on success
 */
int z_pair_get(ZPairID id, double param, double T, ZPairDef *pair);

/**
 * @brief  Create a rational representation for an arbitrary Laplace pair
 *         F(s) = (bₘsᵐ + ... + b₀) / (sⁿ + aₙ₋₁sⁿ⁻¹ + ... + a₀)
 * @param  num_order   Numerator degree
 * @param  num_coeff   Numerator coefficients b₀..bₘ
 * @param  den_order   Denominator degree
 * @param  den_coeff   Denominator coefficients a₀..aₙ₋₁ (monic)
 * @param  result      Output rational function
 * @return 0 on success
 */
int laplace_pair_custom(int num_order, const double *num_coeff,
                         int den_order, const double *den_coeff,
                         LaplaceRational *result);

/**
 * @brief  List all available Laplace pair IDs and their descriptions
 * @param  buffer Output text buffer (caller allocated, min 2048 bytes)
 * @return Number of pairs listed
 */
int laplace_pair_list(char *buffer);

/**
 * @brief  List all available Z-transform pair IDs
 * @param  buffer Output text buffer
 * @return Number of pairs listed
 */
int z_pair_list(char *buffer);

/**
 * @brief  Generate a time signal for a chosen Laplace pair
 *         (numerical evaluation of f(t) for given t in [0, T_max])
 * @param  id        Pair ID
 * @param  param     Pair parameter
 * @param  n_points  Number of time samples
 * @param  T_max     Maximum time
 * @param  sig       Output time signal (caller-allocated values array)
 * @return 0 on success
 */
int laplace_pair_time_signal(LaplacePairID id, double param,
                              int n_points, double T_max,
                              TimeSignal *sig);

/**
 * @brief  Generate a discrete sequence for a chosen Z-transform pair
 * @param  id      Pair ID
 * @param  param   Pair parameter
 * @param  length  Number of samples
 * @param  seq     Output sequence (caller-allocated values array)
 * @return 0 on success
 */
int z_pair_sequence(ZPairID id, double param, int length,
                     DiscreteSignal *seq);

#ifdef __cplusplus
}
#endif

#endif /* TRANSFORM_PAIRS_H */
