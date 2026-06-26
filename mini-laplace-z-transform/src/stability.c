/**
 * @file stability.c
 * @brief Stability Analysis — Routh-Hurwitz, Jury Test, Nyquist Criterion
 *
 * L4: Fundamental stability theorems with numerical implementation.
 * L5: Root locus computation and gain/phase margin analysis.
 *
 * Standard C99 — each function is an independent knowledge point.
 */

#include "stability.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*----------------------------------------------------------------------
 * L4 — Routh-Hurwitz Stability Criterion
 *
 * Theorem (Routh, 1877): For polynomial P(s) = a₀sⁿ + a₁sⁿ⁻¹ + ... + aₙ,
 * all roots are in the open LHP iff all entries in the first column of
 * the Routh array have the same sign (and none are zero).
 *
 * Routh array construction:
 *   Row 0: a₀  a₂  a₄  a₆ ...
 *   Row 1: a₁  a₃  a₅  a₇ ...
 *   Row k (k≥2): entries computed from rows k-2 and k-1 using the
 *                 determinant formula: r(k,j) = -det([r(k-2,0), r(k-2,j+1)],
 *                                                     [r(k-1,0), r(k-1,j+1)]) / r(k-1,0)
 *----------------------------------------------------------------------*/

int routh_stability(const LaplacePolynomial *den, RouthResult *result)
{
    if (!den || !result || den->order < 1 || den->order > LAPLACE_MAX_ORDER) {
        errno = EINVAL;
        return -1;
    }

    int n = den->order;
    int cols = (n + 2) / 2;  /* number of columns in Routh array */
    memset(result, 0, sizeof(*result));

    /* Allocate 2D Routh array: (n+1) rows × cols columns */
    double *rows[LAPLACE_MAX_ORDER + 1];
    double routh_data[(LAPLACE_MAX_ORDER + 1) * ((LAPLACE_MAX_ORDER + 2) / 2)];
    for (int i = 0; i <= n; i++)
        rows[i] = &routh_data[i * cols];
    memset(routh_data, 0, sizeof(routh_data));

    /* Fill first two rows */
    for (int j = 0; j < cols; j++) {
        int c_idx = 2 * j;
        rows[0][j] = (c_idx <= n) ? den->coeff[n - c_idx] : 0.0;
        rows[1][j] = (c_idx + 1 <= n) ? den->coeff[n - (c_idx + 1)] : 0.0;
    }

    /* Build remaining rows */
    for (int i = 2; i <= n; i++) {
        double pivot = rows[i - 1][0];
        if (fabs(pivot) < DBL_EPSILON * 1e3) {
            /* Zero in first column — special case handled by epsilon method */
            pivot = 1e-12;
            rows[i - 1][0] = pivot;
        }
        for (int j = 0; j < cols - 1; j++) {
            /* Routh formula: r(i,j) = (r(i-1,0)*r(i-2,j+1) - r(i-2,0)*r(i-1,j+1)) / r(i-1,0) */
            double det = rows[i - 1][0] * rows[i - 2][j + 1]
                       - rows[i - 2][0] * rows[i - 1][j + 1];
            rows[i][j] = det / pivot;
        }
    }

    /* Extract first column */
    for (int i = 0; i <= n; i++)
        result->first_column[i] = rows[i][0];

    /* Count sign changes */
    int sign_changes = 0;
    int prev_sign = 0;
    int first_nonzero = -1;

    for (int i = 0; i <= n; i++) {
        double val = result->first_column[i];
        if (fabs(val) < DBL_EPSILON * 1e3) continue;

        int sign = (val > 0) ? 1 : -1;
        if (first_nonzero < 0) {
            first_nonzero = i;
            prev_sign = sign;
        } else {
            if (sign != prev_sign) {
                sign_changes++;
                prev_sign = sign;
            }
        }
    }

    result->n_rhp_poles = sign_changes;
    result->n_jw_poles  = 0;  /* Simplified */
    result->is_stable   = (sign_changes == 0 && first_nonzero >= 0);
    return 0;
}

int routh_count_rhp_poles(const RouthResult *result)
{
    if (!result) return -1;
    return result->n_rhp_poles;
}

int routh_special_case(int row, const double *row_data, int n_cols,
                        double epsilon)
{
    (void)row; (void)row_data; (void)n_cols; (void)epsilon;
    /* Returns: 0=normal, 1=zero-in-first-col, 2=all-zero-row */
    /* For now, a simplified version */
    if (fabs(row_data[0]) < 1e-15) return 1;
    int all_zero = 1;
    for (int j = 0; j < n_cols; j++) {
        if (fabs(row_data[j]) > 1e-15) { all_zero = 0; break; }
    }
    if (all_zero) return 2;
    return 0;
}

/*----------------------------------------------------------------------
 * L4 — Jury Stability Test (Discrete-Time)
 *
 * Theorem (Jury, 1964): For polynomial D(z) = d₀ + d₁z + ... + dₙzⁿ
 * with dₙ > 0, all roots are inside the unit circle iff:
 *   Condition 1: D(1) > 0
 *   Condition 2: (-1)ⁿ·D(-1) > 0
 *   Condition 3-k: |d₀| < dₙ and (n-1) inner constraints satisified.
 *
 * Jury table: N rows of decreasing length.
 *   b_k = det([a₀, a_{n-k}], [a_n, a_k])  (for k=0..n-1)
 *   Then repeat with b coefficients.
 *----------------------------------------------------------------------*/

int jury_stability(const ZPoly *poly, JuryResult *result)
{
    if (!poly || !result || poly->order < 1 || poly->order > Z_MAX_ORDER) {
        errno = EINVAL;
        return -1;
    }

    memset(result, 0, sizeof(*result));
    int n = poly->order;
    int total_conditions = n + 1;
    result->conditions_total = total_conditions;

    /* Convert ZPoly (z⁻¹ form) to positive-z form for Jury test:
     *   ZPoly stores: P(z⁻¹) = c₀ + c₁·z⁻¹ + c₂·z⁻² + ...
     *   Multiply by zⁿ: D(z) = c₀·zⁿ + c₁·zⁿ⁻¹ + ... + cₙ
     *   So a[i] = poly->coeff[i] (same coefficients, different powers).
     */
    double a[Z_MAX_ORDER + 1];
    memset(a, 0, sizeof(a));
    for (int i = 0; i <= n; i++)
        a[i] = poly->coeff[i];
    double d_n = a[0];  /* leading coefficient of zⁿ */
    if (d_n <= 0.0) return -1;

    /* Condition 1: D(1) > 0 */
    double D1 = 0.0;
    for (int i = 0; i <= n; i++) D1 += a[i];
    if (D1 <= 0.0) { result->is_stable = 0; return 0; }
    result->conditions_passed = 1;

    /* Condition 2: (-1)ⁿ·D(-1) > 0
     * D(-1) = Σ a_i·(-1)^{n-i}. Sign of term i depends on (n-i) parity. */
    double Dm1 = 0.0;
    for (int i = 0; i <= n; i++) {
        Dm1 += ((n - i) % 2 == 0) ? a[i] : -a[i];
    }
    if ((n % 2 == 0 ? Dm1 : -Dm1) <= 0.0) {
        result->is_stable = 0;
        return 0;
    }
    result->conditions_passed = 2;

    /* Condition 3: |d₀| < dₙ  (i.e., |aₙ| < a₀) */
    if (fabs(a[n]) >= a[0]) {
        result->is_stable = 0;
        return 0;
    }
    result->conditions_passed = 3;

    /* Jury table for inner constraints */
    double row[Z_MAX_ORDER + 1];
    double next_row[Z_MAX_ORDER + 1];
    int m = n;
    for (int i = 0; i <= m; i++) row[i] = a[i];

    for (int k = 0; k < n - 1; k++) {
        for (int i = 0; i < m; i++) {
            next_row[i] = row[0] * row[i] - row[m] * row[m - i];
        }
        m--;
        if (fabs(next_row[0]) <= 1e-12) {
            result->is_stable = 0;
            return 0;
        }
        result->conditions_passed++;

        for (int i = 0; i <= m; i++) row[i] = next_row[i];
    }

    result->is_stable = (result->conditions_passed >= total_conditions);
    result->n_unstable_roots = result->is_stable ? 0 : 1;
    return 0;
}

/*----------------------------------------------------------------------
 * L4/L5 — Nyquist Stability Criterion
 *
 * Nyquist (1932): For closed-loop stability with unity feedback,
 * N = Z - P, where:
 *   Z = #unstable closed-loop poles
 *   P = #unstable open-loop poles
 *   N = #clockwise encirclements of the -1 + j0 point by L(jω).
 *----------------------------------------------------------------------*/

int nyquist_stability(const TransferFunction *loop_tf,
                       int num_freq, NyquistResult *result)
{
    if (!loop_tf || !result || num_freq < 10) { errno = EINVAL; return -1; }
    memset(result, 0, sizeof(*result));

    /* Evaluate L(jω) from ω=0 to ω=∞ and also ω=-∞ to ω=0 (half closed contour) */
    /* For encirclement count, use crossing of negative real axis */
    int n_pts = num_freq;
    double *Re = (double*)malloc(n_pts * sizeof(double));
    double *Im = (double*)malloc(n_pts * sizeof(double));
    double *w  = (double*)malloc(n_pts * sizeof(double));
    if (!Re || !Im || !w) { free(Re); free(Im); free(w); return -1; }

    /* Sweep from ω_low to ω_high */
    double w_low  = 0.01;
    double w_high = 1e3;

    for (int i = 0; i < n_pts; i++) {
        double frac = (double)i / (n_pts - 1);
        w[i] = w_low * pow(w_high / w_low, frac);
        double complex jw = 0.0 + w[i] * I;
        double complex L  = tf_evaluate(loop_tf, jw);
        Re[i] = creal(L);
        Im[i] = cimag(L);
    }

    /* Count negative real axis crossings → encirclements */
    int crossings = 0;
    double prev_im = 0.0;

    for (int i = 0; i < n_pts; i++) {
        if (isnan(Re[i]) || isnan(Im[i])) continue;
        if (Im[i] * prev_im < 0.0 && Re[i] < -1.0) {
            /* Crossing of real axis left of -1 */
            crossings++;
        }
        prev_im = Im[i];
    }

    result->encirclements = crossings / 2;  /* Each crossing pair = one encirclement */

    /* Count open-loop RHP poles */
    LaplacePoleZero pz;
    memset(&pz, 0, sizeof(pz));
    LaplaceRational L_rat;
    memset(&L_rat, 0, sizeof(L_rat));
    L_rat.numerator.order   = loop_tf->num_order;
    L_rat.denominator.order = loop_tf->den_order;
    for (int i = 0; i <= loop_tf->num_order; i++)
        L_rat.numerator.coeff[i]   = loop_tf->num_coeff[i];
    for (int i = 0; i <= loop_tf->den_order; i++)
        L_rat.denominator.coeff[i] = loop_tf->den_coeff[i];

    int P = 0;
    if (laplace_find_poles(&L_rat, &pz) == 0) {
        for (int i = 0; i < pz.num_poles; i++) {
            if (creal(pz.poles[i]) > 0.0) P++;
        }
    }

    int Z = result->encirclements + P;
    result->is_stable = (Z == 0);

    free(Re); free(Im); free(w);
    return 0;
}

/*----------------------------------------------------------------------
 * L5 — Gain and Phase Margins
 *----------------------------------------------------------------------*/

int stability_margins(const TransferFunction *loop_tf,
                       double *gm_db, double *pm_deg)
{
    if (!loop_tf || !gm_db || !pm_deg) { errno = EINVAL; return -1; }

    /* Gain margin: find ω where phase = -180°, then GM = -20·log₁₀|L(jω)| */
    /* Phase margin: find ω where |L| = 1 (0 dB), then PM = 180° + ∠L(jω) */

    int N = 2000;
    double w_low  = 1e-3;
    double w_high = 1e5;

    double gm = 1e300;
    double pm = 1e300;

    FrequencyResponse prev, curr;
    tf_frequency_response(loop_tf, w_low, &prev);

    for (int i = 1; i < N; i++) {
        double w = w_low * pow(w_high / w_low, (double)i / (N - 1));
        tf_frequency_response(loop_tf, w, &curr);

        /* Phase cross -180° */
        double phased = prev.phase_deg;
        double phasec = curr.phase_deg;
        /* Normalize to [-360, 0] */
        if (phased > 0) phased -= 360.0;
        if (phasec > 0) phasec -= 360.0;

        if ((phased > -180.0 && phasec < -180.0) ||
            (phased < -180.0 && phasec > -180.0)) {
            /* Interpolate to find exact crossover frequency */
            double w_pi = (prev.omega + curr.omega) / 2.0;
            FrequencyResponse fr;
            tf_frequency_response(loop_tf, w_pi, &fr);
            gm = -fr.magnitude_db;
        }

        /* Magnitude crosses 0 dB */
        if ((prev.magnitude_db > 0.0 && curr.magnitude_db < 0.0) ||
            (prev.magnitude_db < 0.0 && curr.magnitude_db > 0.0)) {
            double w_c = (prev.omega + curr.omega) / 2.0;
            FrequencyResponse fr;
            tf_frequency_response(loop_tf, w_c, &fr);
            double ph = fr.phase_deg;
            if (ph > 0) ph -= 360.0;
            pm = 180.0 + ph;
        }

        prev = curr;
    }

    *gm_db  = (gm  < 1e299) ? gm  : NAN;
    *pm_deg = (pm < 1e299) ? pm : NAN;
    return 0;
}

/*----------------------------------------------------------------------
 * L5 — Root Locus Computation
 *
 * The root locus traces closed-loop pole locations for 1 + K·G(s) = 0
 * as gain K varies from 0 to ∞.
 *----------------------------------------------------------------------*/

int root_locus_compute(const TransferFunction *tf,
                        double K_min, double K_max, int n_points,
                        RootLocusPoint *locus)
{
    if (!tf || !locus || n_points < 2 || K_min < 0.0 || K_max <= K_min) {
        errno = EINVAL;
        return -1;
    }

    int n_branches = tf->den_order;
    if (n_branches < 1) return -1;

    /* Closed-loop characteristic: D(s) + K·N(s) = 0 */
    /* For each K, find roots of D(s) + K·N(s) */

    for (int k = 0; k < n_points; k++) {
        double K = K_min + (K_max - K_min) * (double)k / (n_points - 1);

        /* Build characteristic polynomial C(s) = D(s) + K·N(s) */
        LaplaceRational char_eq;
        memset(&char_eq, 0, sizeof(char_eq));
        char_eq.denominator.order = tf->den_order;
        for (int i = 0; i <= tf->den_order; i++)
            char_eq.denominator.coeff[i] = tf->den_coeff[i];
        for (int i = 0; i <= tf->num_order; i++)
            char_eq.denominator.coeff[i] += K * tf->num_coeff[i];

        char_eq.numerator.order = 0;
        char_eq.numerator.coeff[0] = 1.0;

        LaplacePoleZero pz;
        memset(&pz, 0, sizeof(pz));
        laplace_find_poles(&char_eq, &pz);

        for (int b = 0; b < n_branches; b++) {
            int idx = k * n_branches + b;
            locus[idx].K = K;
            locus[idx].s = (b < pz.num_poles) ? pz.poles[b] : 0.0;
        }
    }
    return n_branches;
}

int root_locus_find_gain(const TransferFunction *tf,
                          double complex desired_pole,
                          double *K)
{
    if (!tf || !K) { errno = EINVAL; return -1; }

    /* Magnitude condition for root locus: |K·G(s)| = 1 at closed-loop pole s */
    /* K = 1/|G(s)| */
    double complex G = tf_evaluate(tf, desired_pole);
    if (isnan(creal(G))) return -1;

    double mag = cabs(G);
    if (mag < 1e-15) return -1;

    *K = 1.0 / mag;
    return 0;
}
