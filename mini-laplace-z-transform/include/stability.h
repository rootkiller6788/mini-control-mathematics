/**
 * @file stability.h
 * @brief Stability Analysis — Routh-Hurwitz, Jury Test, Nyquist Criterion
 *
 * Knowledge Coverage: L4 (Fundamental Laws), L5 (Algorithms)
 *
 * Core stability theorems implemented:
 *
 *   1. Routh-Hurwitz Criterion (continuous-time):
 *      A polynomial a₀sⁿ + a₁sⁿ⁻¹ + ... + aₙ (a₀ > 0) has all roots in the
 *      open Left Half-Plane (LHP) iff all elements in the first column of the
 *      Routh array are positive.
 *
 *   2. Jury Stability Test (discrete-time):
 *      A polynomial A(z) = a₀zⁿ + a₁zⁿ⁻¹ + ... + aₙ (a₀ > 0) has all roots
 *      inside the unit circle iff all Jury table conditions are satisfied.
 *      Conditions involve n+1 constraints for an n-th order system.
 *
 *   3. Nyquist Stability Criterion:
 *      For a closed-loop system with loop gain L(s), the number of unstable
 *      closed-loop poles is N = Z - P, where Z = encirclements of -1 point
 *      by L(jω) and P = unstable open-loop poles.
 *
 *   4. Gain Margin & Phase Margin (L5):
 *      Quantify relative stability from open-loop frequency response.
 *
 * References:
 *   - Routh (1877), "A Treatise on the Stability of a Given State of Motion"
 *   - Jury (1964), "Theory and Application of the Z-Transform Method"
 *   - Nyquist (1932), "Regeneration Theory"
 *   - Ogata, "Modern Control Engineering" (5th ed.)
 */

#ifndef STABILITY_H
#define STABILITY_H

#include <stddef.h>
#include "laplace_core.h"
#include "transfer_function.h"

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------
 * L1 — Stability Result Types
 *----------------------------------------------------------------------*/

/**
 * @struct RouthResult
 * @brief  Routh array and stability verdict
 */
typedef struct {
    int     is_stable;                    /**< 1 if all first-col > 0   */
    int     n_rhp_poles;                  /**< Number of RHP poles      */
    int     n_jw_poles;                   /**< Number of jω-axis poles  */
    double  first_column[LAPLACE_MAX_ORDER + 1]; /**< First col entries */
} RouthResult;

/**
 * @struct JuryResult
 * @brief  Jury stability test verdict
 */
typedef struct {
    int     is_stable;                    /**< 1 if all conditions met  */
    int     n_unstable_roots;             /**< Roots outside unit circle*/
    int     conditions_passed;            /**< How many of n+1 pass     */
    int     conditions_total;             /**< n+1 total conditions     */
} JuryResult;

/**
 * @struct NyquistResult
 * @brief  Nyquist stability analysis result
 */
typedef struct {
    int     is_stable;                    /**< 1 = closed-loop stable   */
    int     encirclements;                /**< N (CW encirclements of -1)*/
    double  gain_margin_db;               /**< GM in dB                 */
    double  phase_margin_deg;             /**< PM in degrees            */
    double  gain_crossover_freq;          /**< ω where ∠L = -180°      */
    double  phase_crossover_freq;         /**< ω where |L| = 1 (0 dB)  */
} NyquistResult;

/**
 * @struct RootLocusPoint
 * @brief  Single point on the root locus
 */
typedef struct {
    double complex s;                     /**< Closed-loop pole location */
    double         K;                     /**< Gain at this point        */
} RootLocusPoint;

/*----------------------------------------------------------------------
 * L4 — Routh-Hurwitz Stability Criterion
 *----------------------------------------------------------------------*/

/**
 * @brief  Construct Routh array and assess stability
 * @param  den    Denominator polynomial (characteristic equation)
 * @param  result Output Routh analysis
 * @return 0 on success, -1 if polynomial order exceeds MAX
 *
 * L4: The Routh-Hurwitz criterion is a necessary and sufficient condition
 *     for continuous-time LTI system stability.
 *     Algorithm: Generate Routh table row by row.
 *     Complexity: O(n²) for n-th order polynomial.
 *
 *     Theorem (Routh, 1877): All roots of a₀sⁿ + ... + aₙ are in LHP
 *     iff all first-column entries of the Routh array have the same sign.
 */
int routh_stability(const LaplacePolynomial *den, RouthResult *result);

/**
 * @brief  Count sign changes in Routh first column → RHP poles
 * @param  result  Previously computed Routh result
 * @return Number of RHP poles detected
 *
 * L4: Each sign change in the first column corresponds to one RHP pole.
 *     This is the Routh-Hurwitz theorem's secondary claim.
 */
int routh_count_rhp_poles(const RouthResult *result);

/**
 * @brief  Check Routh special cases: zero in first column, entire row of zeros
 * @param  row      Routh array row index where issue occurs
 * @param  row_data Row entries
 * @param  n_cols   Number of columns in this row
 * @param  epsilon  Small perturbation to resolve zero first-column entry
 * @return 0 = normal, 1 = zero-first-col (use epsilon), 2 = all-zero row
 *
 * L5: Handling Routh array special cases:
 *     Case 1: Zero in first column → replace with small ε, continue
 *     Case 2: Entire row zero → use auxiliary polynomial from prior row
 */
int routh_special_case(int row, const double *row_data, int n_cols,
                        double epsilon);

/*----------------------------------------------------------------------
 * L4 — Jury Stability Test (Discrete-Time)
 *----------------------------------------------------------------------*/

/**
 * @brief  Apply Jury stability test to discrete characteristic polynomial
 * @param  poly    Characteristic polynomial D(z) = d₀zⁿ + ... + dₙ, d₀ > 0
 * @param  result  Output Jury verdict
 * @return 0 on success
 *
 * L4: Jury test is the discrete-time counterpart of Routh-Hurwitz.
 *     Theorem (Jury, 1964): All roots of D(z) are inside |z| < 1 iff
 *     n+1 constraints are satisfied (for n-th order).
 *     Complexity: O(n²).
 *
 *     Conditions:
 *     (1) D(1) > 0, (2) (-1)^n·D(-1) > 0,
 *     (3) |dₙ| < d₀, plus inner matrix constraints.
 */
int jury_stability(const ZPoly *poly, JuryResult *result);

/*----------------------------------------------------------------------
 * L4/L5 — Nyquist Criterion and Stability Margins
 *----------------------------------------------------------------------*/

/**
 * @brief  Analyze closed-loop stability using Nyquist criterion
 * @param  loop_tf  Open-loop transfer function L(s)
 * @param  num_freq Number of frequency points for contour evaluation
 * @param  result   Output Nyquist analysis
 * @return 0 on success
 *
 * L4: Nyquist criterion: Z = N + P
 *     Z = #unstable closed-loop poles, P = #unstable open-loop poles,
 *     N = #CW encirclements of -1 point by L(jω).
 *     Complexity: O(N·(n+m)) for N frequency points.
 */
int nyquist_stability(const TransferFunction *loop_tf,
                       int num_freq, NyquistResult *result);

/**
 * @brief  Compute gain margin from open-loop frequency response
 * @param  loop_tf  Open-loop L(s)
 * @param  gm_db    Output: gain margin in dB
 * @param  pm_deg   Output: phase margin in degrees
 * @return 0 on success
 *
 * L5: GM = -20·log₁₀|L(jω_π)| where ∠L(jω_π) = -180°
 *     PM = 180° + ∠L(jω_c) where |L(jω_c)| = 1
 *     Complexity: O(N·(n+m)) via frequency sweep.
 */
int stability_margins(const TransferFunction *loop_tf,
                       double *gm_db, double *pm_deg);

/*----------------------------------------------------------------------
 * L5 — Root Locus Computation
 *----------------------------------------------------------------------*/

/**
 * @brief  Compute root locus points for K ∈ [K_min, K_max]
 * @param  tf        Open-loop transfer function G(s)
 * @param  K_min     Minimum gain
 * @param  K_max     Maximum gain
 * @param  n_points  Number of gain values
 * @param  locus     Output array of pole arrays (caller-allocated)
 * @return Number of branches (number of closed-loop poles)
 *
 * L5: Root locus shows how closed-loop poles move as gain K varies.
 *     1 + K·G(s) = 0 → poles are roots of D(s) + K·N(s) = 0
 *     Complexity: O(n_points · M²) for M poles per iteration.
 */
int root_locus_compute(const TransferFunction *tf,
                        double K_min, double K_max, int n_points,
                        RootLocusPoint *locus);

/**
 * @brief  Find gain K that places a dominant pole at desired location
 * @param  tf             Transfer function
 * @param  desired_pole   Desired closed-loop pole location in s-plane
 * @param  K              Output: required gain
 * @return 0 on success, -1 if cannot place
 *
 * L5: Magnitude condition: |K·G(s)| = 1 at closed-loop pole
 *     → K = |1/G(s)| (for negative feedback)
 */
int root_locus_find_gain(const TransferFunction *tf,
                          double complex desired_pole,
                          double *K);

#ifdef __cplusplus
}
#endif

#endif /* STABILITY_H */
