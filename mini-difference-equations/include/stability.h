/**
 * stability.h — Stability Analysis of Discrete-Time Systems
 *
 * Domain: BIBO stability, asymptotic stability, Jury/Schur-Cohn criteria,
 *         Lyapunov method for discrete-time systems, discretization methods
 *         and their stability preservation properties.
 *
 * Knowledge Level Coverage:
 *   L1: Stability definitions (BIBO, asymptotic, marginal)
 *   L2: Stability concepts (internal vs external, Lyapunov stability)
 *   L3: Gain/phase margins (discrete), dominant pole analysis
 *   L4: Jury criterion (necessary & sufficient), Schur-Cohn test
 *   L5: Jury table computation, eigenvalue check, bilinear transform,
 *       discretization methods (Forward Euler, Backward Euler, Tustin, ZOH)
 *
 * Reference: Jury, E.I. "Inners and Stability of Dynamic Systems" (1974)
 * Reference: Åström & Wittenmark "Computer-Controlled Systems" (1997)
 * Reference: Ogata, K. "Discrete-Time Control Systems" (1995)
 *
 * Course Mapping:
 *   MIT 6.241 — Dynamic Systems (discrete stability)
 *   Caltech CDS 110 — Linear Systems (Jury test)
 *   ETH 151-0563 — Robust Control (discrete Lyapunov)
 *   Cambridge 3F2 — Systems & Control (discrete stability criteria)
 */

#ifndef STABILITY_H
#define STABILITY_H

#include <stddef.h>
#include <complex.h>
#include "ztransform.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ── L1: Stability Type Definitions ──────────────────────────────────── */

/**
 * Stability classifications for discrete LTI systems.
 */
typedef enum {
    STABLE_ASYMPTOTIC,   /* All |λ_i| < 1 → decays to zero */
    STABLE_MARGINAL,     /* All |λ_i| ≤ 1, some = 1 → bounded, no decay */
    UNSTABLE,            /* Any |λ_i| > 1 → grows unbounded */
    STABLE_BIBO,         /* Bounded-input bounded-output stable */
    STABLE_LYAPUNOV,     /* Lyapunov stable (i.s.L.) */
    UNKNOWN_STABILITY
} StabilityType;

/**
 * Jury stability table for testing if all roots of:
 *   P(z) = a_0 z^n + a_1 z^{n-1} + ... + a_n  (a_0 > 0)
 * lie inside the unit circle |z| < 1.
 *
 * Table structure (n+1 columns, 2n-3 rows for n>2):
 *   Row 0: [a_0  a_1  ...  a_n]
 *   Row 1: [a_n  a_{n-1}  ...  a_0]  (reversed)
 *   Row 2k: b_k = det([row(2k-2)[0], row(2k-2)[n+1-k]; ...])
 *   Row 2k+1: reversed of row 2k
 */
typedef struct {
    int       order;        /* Polynomial order n */
    double   *coeffs;       /* Original coefficients [a_0..a_n] (n+1) */
    double  **table;        /* Jury table rows (2n-3 rows, varying columns) */
    int       num_rows;     /* Number of table rows = 2n-3 */
    int       is_stable;    /* 1 if Jury test passes, 0 otherwise */
} JuryTable;

/**
 * Discretization method enumeration.
 */
typedef enum {
    DISC_FORWARD_EULER,    /* s ← (z-1)/T */
    DISC_BACKWARD_EULER,   /* s ← (z-1)/(zT) */
    DISC_TUSTIN,           /* s ← (2/T)·(z-1)/(z+1) — bilinear */
    DISC_ZOH,              /* Zero-order hold: z = e^{sT} */
    DISC_MATCHED_POLE_ZERO /* Matched pole-zero mapping */
} DiscretizationMethod;

/**
 * Continuous-time transfer function (for discretization):
 *   G(s) = (b_0 s^m + ... + b_m) / (a_0 s^n + ... + a_n)
 */
typedef struct {
    int     num_order;   /* Numerator order m */
    int     den_order;   /* Denominator order n */
    double *num_coeffs;  /* [b_0, b_1, ..., b_m], length m+1 */
    double *den_coeffs;  /* [a_0, a_1, ..., a_n], length n+1 */
} ContTransferFunc;

/**
 * Discrete-time root locus parameters.
 */
typedef struct {
    int             num_poles;
    int             num_zeros;
    double complex *poles;
    double complex *zeros;
    double          gain_min;
    double          gain_max;
    int             num_gain_steps;
    double complex *locus_points;  /* Closed-loop poles for each gain */
    int             num_locus_points;
    double          gain_margin;   /* GM in dB */
    double          phase_margin;  /* PM in degrees */
} DiscRootLocus;

/* ── API: L4/L5 — Jury Stability Criterion ───────────────────────────── */

/**
 * Build the Jury stability table for polynomial P(z).
 *
 * Jury's Necessary Conditions (n=order):
 *   1. P(1) > 0
 *   2. (-1)^n P(-1) > 0
 *   3. |a_0| > |a_n| (for stability, a_0 > 0 is assumed)
 *
 * Jury's Sufficient Conditions (from table):
 *   |b_0| > |b_{n-1}|, |c_0| > |c_{n-2}|, ... decreasing column count.
 *
 * Complexity: O(n²)
 */
int jury_table_build(int order, const double *coeffs, JuryTable *jt);

/**
 * Run the complete Jury stability test.
 * Returns 1 if all roots are inside the unit circle, 0 otherwise.
 * Complexity: O(n²)
 */
int jury_stability_test(int order, const double *coeffs);

/**
 * Check the necessary conditions only (quick pre-screen).
 * Returns 1 if passed, 0 if violated (unstable guaranteed).
 * Complexity: O(n)
 */
int jury_necessary_conditions(int order, const double *coeffs);

/**
 * Free memory for JuryTable.
 */
void jury_table_free(JuryTable *jt);

/* ── API: L4/L5 — Eigenvalue-Based Stability ─────────────────────────── */

/**
 * Check stability by computing eigenvalues of companion matrix.
 * All |λ_i| < 1 → asymptotically stable.
 * Complexity: O(n³) for eigenvalue computation
 */
StabilityType eigenvalue_stability_check(int order, const double *a_coeffs);

/**
 * Compute spectral radius ρ(A) = max |λ_i(A)|.
 * Complexity: O(n³)
 */
double spectral_radius(const double *A, int n);

/**
 * Check BIBO stability: system is BIBO stable iff
 * Σ_{n=0}^{∞} |h[n]| < ∞ (impulse response absolutely summable).
 * Complexity: O(num_samples)
 */
int bibo_stability_check(const double *impulse_response, int num_samples,
                         double tol);

/**
 * Check asymptotic stability from impulse response decay.
 * Complexity: O(num_samples)
 */
int asymptotic_stability_check(const double *impulse_response,
                               int num_samples, double tol);

/* ── API: L5 — Bilinear (Tustin) Transform ───────────────────────────── */

/**
 * Bilinear transform: map continuous s to discrete z.
 *   s = (2/T) · (z-1)/(z+1)
 *   z = (1 + sT/2) / (1 - sT/2)
 *
 * Stability preserved: left-half s-plane ↔ inside unit circle in z-plane.
 * Complexity: O(1) for point mapping
 */
double complex bilinear_s_to_z(double complex s, double T);
double complex bilinear_z_to_s(double complex z, double T);

/**
 * Verify that bilinear transform preserves stability:
 * Re(s) < 0  ↔  |z| < 1
 * Complexity: O(1)
 */
int bilinear_stability_preservation(double complex s, double T, double tol);

/* ── API: L5 — Discretization Methods ─────────────────────────────────── */

/**
 * Forward Euler discretization: s ← (z-1)/T
 * Maps stable s-poles only if |1 + sT| < 1 → limited stability region.
 * Complexity: O(m + n)
 */
int forward_euler_discretize(const ContTransferFunc *ctf, double T,
                             ZTransferFunc *dtf);

/**
 * Backward Euler discretization: s ← (z-1)/(zT)
 * Maps entire left-half plane inside unit circle → A-stable.
 * Complexity: O(m + n)
 */
int backward_euler_discretize(const ContTransferFunc *ctf, double T,
                               ZTransferFunc *dtf);

/**
 * Tustin (bilinear) discretization: s ← (2/T)·(z-1)/(z+1)
 * A-stable, frequency warping occurs at high frequencies.
 * Complexity: O(m·n) for polynomial manipulation
 */
int tustin_discretize(const ContTransferFunc *ctf, double T,
                      ZTransferFunc *dtf);

/**
 * Zero-Order Hold (ZOH) discretization:
 *   H(z) = (1 - z^{-1}) · Z{L^{-1}{G(s)/s}_{t=kT}}
 * Exact for step-invariant discretization.
 * Complexity: O(n³) for pole computations
 */
int zoh_discretize(const ContTransferFunc *ctf, double T,
                   ZTransferFunc *dtf);

/**
 * Matched pole-zero mapping: z_k = e^{s_k T} for all poles and zeros.
 * Complexity: O(m + n)
 */
int matched_pole_zero_discretize(const ContTransferFunc *ctf, double T,
                                  ZTransferFunc *dtf);

/* ── API: L3 — Discrete Margins ──────────────────────────────────────── */

/**
 * Compute discrete gain margin from Nyquist plot of open-loop H(e^{jω}).
 * GM = 1/|H(e^{jω_pc})| where ω_pc is the phase crossover (∠H = -π).
 * Complexity: O(num_freq_points)
 */
double discrete_gain_margin(const FreqResponse *ol_freq_resp);

/**
 * Compute discrete phase margin from open-loop frequency response.
 * PM = π + ∠H(e^{jω_gc}) where ω_gc is the gain crossover (|H| = 1).
 * Complexity: O(num_freq_points)
 */
double discrete_phase_margin(const FreqResponse *ol_freq_resp);

/* ── API: L5 — Schur-Cohn Stability Test ─────────────────────────────── */

/**
 * Schur-Cohn test: an alternative to Jury test using determinants
 * of Schur-Cohn matrices. Computationally heavier than Jury but
 * conceptually simpler.
 * Complexity: O(n³)
 */
int schur_cohn_stability_test(int order, const double *coeffs);

/* ── Utility ─────────────────────────────────────────────────────────── */

void cont_transfer_free(ContTransferFunc *ctf);
ContTransferFunc* cont_transfer_alloc(int num_order, int den_order);
void disc_root_locus_free(DiscRootLocus *rl);

/**
 * Print stability verdict with reasoning.
 */
void stability_print_verdict(StabilityType type);

/**
 * Print Jury table.
 */
void jury_table_print(const JuryTable *jt);

#ifdef __cplusplus
}
#endif

#endif /* STABILITY_H */
