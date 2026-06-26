/**
 * @file laplace_core.h
 * @brief Laplace Transform — Core Definitions and Properties
 *
 * Knowledge Coverage: L1 (Definitions), L2 (Core Concepts), L4 (Fundamental Laws)
 *
 * The Laplace transform maps a time-domain signal f(t) to the complex s-plane:
 *   F(s) = L{f(t)} = ∫₀^∞ f(t) · e^(-st) dt
 *
 * Key properties implemented:
 *   - Linearity:  L{a·f(t) + b·g(t)} = a·F(s) + b·G(s)
 *   - Time shift: L{f(t - τ)·u(t - τ)} = e^(-sτ)·F(s)
 *   - Frequency shift: L{e^(at)·f(t)} = F(s - a)
 *   - Convolution: L{f(t) * g(t)} = F(s)·G(s)
 *   - Differentiation: L{f'(t)} = s·F(s) - f(0⁻)
 *   - Integration: L{∫₀^t f(τ)dτ} = F(s)/s
 *   - Initial Value: f(0⁺) = lim_{s→∞} s·F(s)
 *   - Final Value:   f(∞) = lim_{s→0} s·F(s)  (if poles in LHP)
 *
 * Reference textbooks:
 *   - Oppenheim & Willsky, "Signals and Systems" (2nd ed.)
 *   - Dorf & Bishop, "Modern Control Systems" (13th ed.)
 */

#ifndef LAPLACE_CORE_H
#define LAPLACE_CORE_H

#include <stddef.h>   /* size_t */
#include <complex.h>  /* C99 complex types */

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------
 * L1 — Core Type Definitions
 *----------------------------------------------------------------------*/

/** Maximum order of rational transfer functions supported */
#define LAPLACE_MAX_ORDER  64

/**
 * @struct LaplacePolynomial
 * @brief  Represents a polynomial in s: P(s) = a₀ + a₁·s + ... + aₙ·sⁿ
 *
 * L1: Polynomial representation is the fundamental data structure for
 *     rational Laplace-domain functions.
 */
typedef struct {
    int      order;                        /**< Degree n of the polynomial   */
    double   coeff[LAPLACE_MAX_ORDER + 1]; /**< Coefficients a₀..aₙ        */
} LaplacePolynomial;

/**
 * @struct LaplaceRational
 * @brief  Rational function in s-plane: F(s) = N(s) / D(s)
 *
 * L1: Every proper rational function in control theory is expressed this way.
 */
typedef struct {
    LaplacePolynomial numerator;   /**< N(s) numerator polynomial   */
    LaplacePolynomial denominator; /**< D(s) denominator polynomial */
} LaplaceRational;

/**
 * @enum LaplaceROCTypes
 * @brief Region of Convergence type classification
 *
 * L1: The ROC determines uniqueness of the inverse Laplace transform.
 *     For rational transforms, ROC is bounded by poles.
 */
typedef enum {
    LAPLACE_ROC_RIGHT_SIDED = 0,  /**< ROC: Re{s} > σ₀ (causal signals)    */
    LAPLACE_ROC_LEFT_SIDED  = 1,  /**< ROC: Re{s} < σ₀ (anti-causal)       */
    LAPLACE_ROC_STRIP       = 2,  /**< ROC: σ₁ < Re{s} < σ₂ (two-sided)     */
    LAPLACE_ROC_ENTIRE      = 3,  /**< ROC: entire s-plane (finite-duration)*/
    LAPLACE_ROC_EMPTY       = 4   /**< No ROC exists (no inverse)           */
} LaplaceROCTypes;

/**
 * @struct LaplacePoleZero
 * @brief  Pole-zero description of a transfer function
 *
 * L2: Poles and zeros characterize system behavior — stability, transient
 *     response, and frequency response.
 */
typedef struct {
    int      num_poles;                              /**< Number of poles              */
    int      num_zeros;                              /**< Number of zeros              */
    double complex poles[LAPLACE_MAX_ORDER];         /**< Pole locations in s-plane    */
    double complex zeros[LAPLACE_MAX_ORDER];         /**< Zero locations in s-plane    */
    double   gain;                                   /**< DC gain K                    */
} LaplacePoleZero;

/**
 * @struct TimeSignal
 * @brief  Sampled representation of a continuous-time signal f(t)
 *
 * L1: Time-domain signal representation for numerical computation.
 */
typedef struct {
    int      n_samples;                 /**< Number of time samples         */
    double   t_start;                   /**< Start time t₀                  */
    double   t_step;                    /**< Sampling interval Δt           */
    double  *values;                    /**< f(t₀ + k·Δt), k=0..n_samples-1 */
} TimeSignal;

/*----------------------------------------------------------------------
 * L2 — Core Property Verification Functions
 *----------------------------------------------------------------------*/

/**
 * @brief  Evaluate polynomial P(s) at complex s = σ + jω
 * @param  poly  Pointer to polynomial structure
 * @param  s     Complex frequency point
 * @return P(s)  = Σ aₖ · s^k
 *
 * L2: Polynomial evaluation is the most fundamental operation for
 *     frequency response computation and pole-zero analysis.
 *     Algorithm: Horner's method, O(n) complex multiplications.
 */
double complex laplace_poly_eval(const LaplacePolynomial *poly,
                                  double complex s);

/**
 * @brief  Evaluate rational function F(s) = N(s)/D(s) at complex s
 * @param  rat   Pointer to rational function
 * @param  s     Complex frequency point
 * @return F(s)  value; NaN if D(s) = 0 (pole location)
 *
 * L4: Implements the fundamental transfer function evaluation.
 *     Returns NaN with errno set for pole locations.
 *     Complexity: O(n+m) where n,m are numerator/denominator orders.
 */
double complex laplace_rational_eval(const LaplaceRational *rat,
                                      double complex s);

/**
 * @brief  Verify linearity property: L{a·f + b·g} = a·L{f} + b·L{g}
 * @param  a,b        Scalar multipliers
 * @param  Fs, Gs     Laplace transforms F(s), G(s) at point s
 * @param  combined   L{a·f + b·g} computed directly
 * @param  tolerance  Acceptable floating-point error bound
 * @return 1 if |combined - (a·F + b·G)| < tolerance, 0 otherwise
 *
 * L4: Linearity is the foundational property of all linear time-invariant
 *     (LTI) system analysis. This function provides numerical verification.
 */
int laplace_verify_linearity(double a, double b,
                              double complex Fs, double complex Gs,
                              double complex combined,
                              double tolerance);

/**
 * @brief  Compute time-shift factor in s-domain: e^{-sτ}
 * @param  s     Complex frequency
 * @param  tau   Time shift τ (positive for delay)
 * @return e^{-sτ}  Time-shift factor
 *
 * L2: Time shifting ↔ multiplication by exponential in s-domain.
 */
double complex laplace_time_shift_factor(double complex s, double tau);

/**
 * @brief  Compute frequency-shifted transform: F(s - a)
 * @param  F     Original rational function F(s)
 * @param  a     Frequency shift (real, can be positive or negative)
 * @param  result Output rational function F(s - a)
 * @return 0 on success, -1 if result order exceeds LAPLACE_MAX_ORDER
 *
 * L2: Frequency shifting is used for exponential modulation analysis.
 *     Complexity: O(n²) due to binomial expansion of (s - a)^k.
 */
int laplace_frequency_shift(const LaplaceRational *F,
                             double a,
                             LaplaceRational *result);

/**
 * @brief  Apply Initial Value Theorem: f(0⁺) = lim_{s→∞} s·F(s)
 * @param  F     Rational function F(s) in strictly proper form
 * @return f(0⁺) initial value; NaN if not strictly proper
 *
 * L4: The Initial Value Theorem connects s→∞ behavior to t=0⁺ behavior.
 *     Requires F(s) to be strictly proper (numerator order < denominator order).
 */
double laplace_initial_value(const LaplaceRational *F);

/**
 * @brief  Apply Final Value Theorem: f(∞) = lim_{s→0} s·F(s)
 * @param  F     Rational function F(s)
 * @return f(∞)  final value; NaN if any pole in RHP or on jω-axis (except s=0)
 *
 * L4: The Final Value Theorem provides steady-state value without inverse transform.
 *     Valid only if all poles of s·F(s) lie in the open left half-plane.
 *     Complexity: O(n) — checks pole locations then evaluates limit.
 */
double laplace_final_value(const LaplaceRational *F);

/**
 * @brief  Apply differentiation property: L{f'(t)} = s·F(s) - f(0⁻)
 * @param  F     F(s) = L{f(t)}
 * @param  f0    Initial condition f(0⁻)
 * @param  dF    Output: L{f'(t)} (numerator order may increase by 1)
 * @return 0 on success, -1 on overflow
 *
 * L2: Differentiation in time domain ↔ multiplication by s in s-domain.
 */
int laplace_differentiation(const LaplaceRational *F,
                             double f0,
                             LaplaceRational *dF);

/**
 * @brief  Apply integration property: L{∫₀^t f(τ)dτ} = F(s)/s
 * @param  F        F(s) = L{f(t)}
 * @param  integral Output: L{∫₀^t f(τ)dτ}
 * @return 0 on success, -1 on overflow
 *
 * L2: Integration in time domain ↔ division by s in s-domain.
 */
int laplace_integration(const LaplaceRational *F,
                         LaplaceRational *integral);

/**
 * @brief  Convolution theorem: L{f * g} = F(s)·G(s)
 * @param  F        F(s) = L{f(t)}
 * @param  G        G(s) = L{g(t)}
 * @param  conv     Output: F(s)·G(s) (product of rationals)
 * @return 0 on success, -1 on overflow
 *
 * L2/L4: Convolution in time ↔ multiplication in s-domain.
 *        This is arguably the most important property for LTI systems.
 *        Complexity: O(n·m) for polynomial multiplication.
 */
int laplace_convolution(const LaplaceRational *F,
                         const LaplaceRational *G,
                         LaplaceRational *conv);

/*----------------------------------------------------------------------
 * L3 — Polynomial Operations (Mathematical Structures)
 *----------------------------------------------------------------------*/

/**
 * @brief  Multiply two polynomials: R(s) = P(s) · Q(s)
 * @param  P,Q    Input polynomials
 * @param  R      Result polynomial (pre-allocated)
 * @return 0 on success, -1 if result order > LAPLACE_MAX_ORDER
 *
 * L3: Polynomial multiplication is the core operation for combining
 *     transfer functions in series.
 *     Complexity: O(np·nq) convolution of coefficient arrays.
 */
int laplace_poly_multiply(const LaplacePolynomial *P,
                           const LaplacePolynomial *Q,
                           LaplacePolynomial *R);

/**
 * @brief  Add two polynomials: R(s) = P(s) + Q(s)
 * @param  P,Q    Input polynomials
 * @param  R      Result polynomial
 * @return 0 on success
 */
int laplace_poly_add(const LaplacePolynomial *P,
                      const LaplacePolynomial *Q,
                      LaplacePolynomial *R);

/**
 * @brief  Polynomial derivative: P'(s) = dP/ds
 * @param  P      Input polynomial
 * @param  dP     Derivative polynomial
 * @return 0 on success
 *
 * L3: Polynomial derivative is used in residue computation for
 *     partial fraction expansion.
 */
int laplace_poly_derivative(const LaplacePolynomial *P,
                             LaplacePolynomial *dP);

/**
 * @brief  Compute s-plane pole locations by finding roots of denominator
 * @param  F     Rational function
 * @param  pz    Output pole-zero structure (poles filled in)
 * @return 0 on success, -1 if polynomial order > LAPLACE_MAX_ORDER
 *
 * L3/L5: Root-finding for denominator polynomial.
 *        Up to 4th order: closed-form (quadratic/cubic/quartic).
 *        Higher orders: Aberth-Ehrlich iteration.
 *        Complexity: O(n²) per iteration for n-order polynomial.
 */
int laplace_find_poles(const LaplaceRational *F,
                        LaplacePoleZero *pz);

/**
 * @brief  Compute zeros by finding roots of numerator polynomial
 * @param  F     Rational function
 * @param  pz    Pole-zero structure (zeros filled in)
 * @return 0 on success
 *
 * L3: Mirror of pole-finding for the numerator.
 */
int laplace_find_zeros(const LaplaceRational *F,
                        LaplacePoleZero *pz);

/**
 * @brief  Compute Laplace transform at a single point s numerically
 *         using composite Simpson's rule for the integral ∫₀^∞ f(t)e^{-st}dt
 * @param  sig     Pointer to sampled time signal f(t)
 * @param  s       Complex frequency s = σ + jω
 * @param  T_max   Upper integration limit (truncation)
 * @return Approximation of F(s)
 *
 * L3: Numerical Laplace transform computation for arbitrary f(t).
 *     Complexity: O(N) where N is number of samples.
 */
double complex laplace_numerical_transform(const TimeSignal *sig,
                                            double complex s,
                                            double T_max);

/*----------------------------------------------------------------------
 * L4 — Inverse Laplace (Bromwich Integral Approximation)
 *----------------------------------------------------------------------*/

/**
 * @brief  Numerical inverse Laplace transform using Fourier-series method
 *         (Dubner-Abate-Crump algorithm).
 * @param  F          Rational function F(s)
 * @param  t          Time point t > 0
 * @param  n_terms    Number of Fourier-series terms
 * @return Approximation of f(t)
 *
 * L4: Inverse Laplace via Bromwich integral approximation.
 *     f(t) ≈ (e^{σt}/T) · [F(σ)/2 + Σ_{k=1}^N Re{F(σ + jkπ/T)·e^{jkπt/T}}]
 *     Complexity: O(N) per time point.
 */
double laplace_inverse_numerical(const LaplaceRational *F,
                                  double t,
                                  int n_terms);

#ifdef __cplusplus
}
#endif

#endif /* LAPLACE_CORE_H */
