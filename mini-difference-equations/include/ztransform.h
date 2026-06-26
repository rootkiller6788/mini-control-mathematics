/**
 * ztransform.h — Z-Transform for Discrete-Time Systems
 *
 * Domain: The z-transform is the discrete-time counterpart of the
 * Laplace transform. It maps a discrete sequence x[n] to a complex
 * function X(z) = Σ x[n] z^{-n}, enabling algebraic analysis of
 * difference equations and discrete LTI systems.
 *
 * Knowledge Level Coverage:
 *   L1: z-transform definition, region of convergence (ROC)
 *   L2: Transfer function H(z) = B(z)/A(z) as pulse transfer function
 *   L3: Pole-zero representation, frequency response H(e^{jωT})
 *   L4: Final/Initial Value Theorems, Convolution Theorem in z-domain
 *   L5: Inverse z-transform (partial fractions, long division, contour)
 *
 * Reference: Jury, E.I. "Theory and Application of the z-Transform" (1964)
 * Reference: Ogata, K. "Discrete-Time Control Systems" (1995)
 * Reference: Oppenheim & Schafer "Discrete-Time Signal Processing" (2010)
 *
 * Course Mapping:
 *   MIT 6.003 — Signals and Systems (z-transform)
 *   Stanford EE 264 — Digital Signal Processing
 *   Berkeley EE 123 — Digital Signal Processing
 *   ETH 227-0427 — Signal and System Theory II
 */

#ifndef ZTRANSFORM_H
#define ZTRANSFORM_H

#include <stddef.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── L1: Z-Transform Definition ──────────────────────────────────────── */

/**
 * Z-transform pair: X(z) = Σ_{n=0}^{∞} x[n] z^{-n}
 *
 * The z-transform maps a discrete-time sequence to a complex function.
 * ROC (Region of Convergence): |z| > R for causal sequences.
 */
typedef struct {
    double complex *x;      /* Sequence values x[0], x[1], ..., x[N-1] */
    int             length; /* Sequence length */
    double          roc_radius; /* |z| > roc_radius for causal sequences */
} ZSequence;

/**
 * Z-domain transfer function (pulse transfer function):
 *   H(z) = B(z)/A(z)
 *        = (b_0 + b_1 z^{-1} + ... + b_m z^{-m}) /
 *          (1   + a_1 z^{-1} + ... + a_n z^{-n})
 *
 * In positive powers of z:
 *   H(z) = (b_0 z^n + b_1 z^{n-1} + ... + b_m z^{n-m}) /
 *          (z^n + a_1 z^{n-1} + ... + a_n)
 */
typedef struct {
    int             num_order;    /* Numerator order m */
    int             den_order;    /* Denominator order n */
    double         *num_coeffs;   /* Numerator b[0..m] (length m+1) */
    double         *den_coeffs;   /* Denominator a[1..n] (length n, a[0]=1) */
    double complex *zeros;        /* Zeros (roots of B(z)) */
    double complex *poles;        /* Poles (roots of A(z)) */
    int             num_zeros;
    int             num_poles;
    double          gain;         /* DC gain = H(1) if stable */
} ZTransferFunc;

/**
 * Partial fraction expansion term:
 *   H(z) = Σ r_k / (1 - p_k z^{-1}) + polynomial part
 *
 * For inverse z-transform:
 *   If p_k is real:     r_k · p_k^n · u[n]
 *   If p_k is complex:  2|r_k||p_k|^n cos(n∠p_k + ∠r_k) · u[n]
 */
typedef struct {
    double complex *residues;   /* r_k for each pole */
    double complex *poles;      /* p_k pole locations */
    int             num_terms;  /* Number of partial fraction terms */
    int             has_direct; /* 1 if direct feedthrough term exists */
    double          direct_term; /* Direct term (for proper/improper TF) */
} PartialFractionExpansion;

/**
 * Inverse z-transform method enumeration.
 */
typedef enum {
    IZT_PARTIAL_FRACTION,  /* Partial fraction expansion + table lookup */
    IZT_LONG_DIVISION,     /* Power series / long division */
    IZT_CONTOUR_INTEGRAL,  /* Residue theorem / contour integration */
    IZT_DIFFERENCE_EQ      /* Solve difference equation recursively */
} InverseZMethod;

/**
 * Frequency response: H(e^{jω}) at discrete frequencies.
 */
typedef struct {
    double complex *response;   /* H(e^{jω}) values */
    double         *magnitude;   /* |H(e^{jω})| in linear scale */
    double         *phase;       /* ∠H(e^{jω}) in radians */
    double         *frequencies; /* ω values in rad/sample */
    int             num_points;  /* Number of frequency points */
} FreqResponse;

/* ── API: L5 — Direct Z-Transform ────────────────────────────────────── */

/**
 * Compute z-transform X(z) at a specific complex point z.
 * X(z) = Σ_{n=0}^{N-1} x[n] z^{-n}
 * Complexity: O(N)
 */
double complex ztransform_at(const ZSequence *seq, double complex z);

/**
 * Compute z-transform over a vector of z values.
 * Complexity: O(N * num_z)
 */
void ztransform_vector(const ZSequence *seq,
                       const double complex *z_values, int num_z,
                       double complex *Xz);

/**
 * Evaluate transfer function H(z) at a specific complex point z.
 * Complexity: O(max(num_order, den_order))
 */
double complex ztransfer_eval(const ZTransferFunc *tf, double complex z);

/* ── API: L5 — Inverse Z-Transform ───────────────────────────────────── */

/**
 * Inverse z-transform via partial fraction expansion.
 * Decomposes X(z) → Σ r_k/(1 - p_k z^{-1}) then inverts via table.
 * Complexity: O(n³) for pole-finding, O(n²) for residue computation
 */
int inverse_ztransform_partial_fraction(const ZTransferFunc *tf,
                                        int num_samples,
                                        double *y);

/**
 * Inverse z-transform via long division (power series method).
 * Computes y[0], y[1], ..., y[N-1] by polynomial division B(z)/A(z).
 * Complexity: O(N * den_order)
 */
int inverse_ztransform_long_division(const ZTransferFunc *tf,
                                     int num_samples,
                                     double *y);

/**
 * Perform partial fraction expansion of H(z) = B(z)/A(z).
 * Handles real distinct poles, repeated poles, and complex conjugate pairs.
 * Complexity: O(n³)
 */
int partial_fraction_expand(const ZTransferFunc *tf,
                            PartialFractionExpansion *pfe);

/**
 * Free memory allocated for PartialFractionExpansion.
 */
void pfe_free(PartialFractionExpansion *pfe);

/* ── API: L3 — Frequency Response ────────────────────────────────────── */

/**
 * Compute frequency response H(e^{jω}) at N equally spaced points
 * from ω = 0 to π (Nyquist frequency).
 * Complexity: O(N * max(num_order, den_order))
 */
int ztransfer_freq_response(const ZTransferFunc *tf,
                            int num_points,
                            FreqResponse *fr);

/**
 * Find DC gain = H(1) and high-frequency gain = H(-1) if they exist.
 * Complexity: O(max(num_order, den_order))
 */
void ztransfer_dc_gain(const ZTransferFunc *tf,
                       double *dc_gain, double *hf_gain);

/**
 * Compute magnitude in dB from linear magnitude.
 * Complexity: O(num_points)
 */
void freq_response_to_db(FreqResponse *fr);

/**
 * Free memory for FreqResponse.
 */
void freq_response_free(FreqResponse *fr);

/* ── API: L4 — Z-Transform Theorems ──────────────────────────────────── */

/**
 * Verify time-shift property:
 *   Z{x[n-k]} = z^{-k} X(z)  (for k > 0, causal shift)
 * Complexity: O(N * num_z_test_points)
 */
int ztransform_verify_time_shift(const ZSequence *seq, int k,
                                 int num_test_points, double tol);

/**
 * Verify linearity property:
 *   Z{αx₁[n] + βx₂[n]} = αZ{x₁[n]} + βZ{x₂[n]}
 * Complexity: O(N * num_z)
 */
int ztransform_verify_linearity(const ZSequence *seq1,
                                const ZSequence *seq2,
                                double alpha, double beta,
                                const double complex *z_vals, int num_z,
                                double tol);

/**
 * Apply final value theorem:
 *   lim_{n→∞} x[n] = lim_{z→1} (z-1)X(z)
 * Provided the limit exists (all poles of (z-1)X(z) inside unit circle).
 * Complexity: O(den_order)
 */
double ztransform_final_value(const ZTransferFunc *tf);

/**
 * Apply initial value theorem:
 *   x[0] = lim_{z→∞} X(z)
 * Complexity: O(1)
 */
double ztransform_initial_value(const ZTransferFunc *tf);

/**
 * Verify convolution theorem:
 *   Z{x₁[n] * x₂[n]} = X₁(z) · X₂(z)
 * Complexity: O(N * num_z)
 */
int ztransform_verify_convolution(const ZSequence *seq1,
                                  const ZSequence *seq2,
                                  const double complex *z_vals, int num_z,
                                  double tol);

/* ── API: L3 — Pole/Zero Analysis ────────────────────────────────────── */

/**
 * Compute poles and zeros of transfer function by finding roots of
 * A(z) and B(z) using companion matrix eigenvalue method.
 * Complexity: O(n³)
 */
int ztransfer_pole_zero_analysis(ZTransferFunc *tf);

/**
 * Check if transfer function is minimum-phase (all zeros inside unit circle).
 * Complexity: O(num_zeros)
 */
int ztransfer_is_minimum_phase(const ZTransferFunc *tf);

/**
 * Check if transfer function is stable (all poles inside unit circle).
 * Complexity: O(num_poles)
 */
int ztransfer_is_stable(const ZTransferFunc *tf);

/**
 * Map continuous-time poles s_k to discrete-time poles z_k via:
 *   z_k = e^{s_k T}  (exact discretization with ZOH)
 * Complexity: O(num_poles)
 */
void continuous_to_discrete_poles(const double complex *s_poles,
                                  int num_poles, double T,
                                  double complex *z_poles);

/**
 * Map discrete-time poles z_k to continuous-time poles s_k via:
 *   s_k = (1/T) ln(z_k)
 * Complexity: O(num_poles)
 */
void discrete_to_continuous_poles(const double complex *z_poles,
                                  int num_poles, double T,
                                  double complex *s_poles);

/* ── Utility ─────────────────────────────────────────────────────────── */

void ztransfer_free(ZTransferFunc *tf);
void zsequence_free(ZSequence *seq);
ZTransferFunc* ztransfer_alloc(int num_order, int den_order);
ZSequence* zsequence_alloc(int length);

/**
 * Set coefficients for transfer function. Copies data internally.
 */
int ztransfer_set_coeffs(ZTransferFunc *tf,
                         const double *num, const double *den);

/**
 * Print transfer function H(z) in human-readable form.
 */
void ztransfer_print(const ZTransferFunc *tf);

#ifdef __cplusplus
}
#endif

#endif /* ZTRANSFORM_H */
