/**
 * @file z_transform.c
 * @brief Z-Transform — Core implementations
 *
 * Implements: discrete polynomial arithmetic, Z-domain evaluation,
 *             shift properties, initial/final value, inverse Z-transform
 *             methods (power series, partial fraction, residue).
 *
 * Standard C99 — each function is an independent knowledge point.
 */

#include "z_transform.h"
#include "laplace_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*----------------------------------------------------------------------
 * L2 — Z-Polynomial Evaluation
 *----------------------------------------------------------------------*/

double complex z_poly_eval(const ZPoly *poly, double complex z)
{
    if (!poly || poly->order < 0 || poly->order > Z_MAX_ORDER) {
        errno = EINVAL;
        return NAN + I * NAN;
    }

    /* P(z) = Σ b_k · z^(-k) using Horner in z^(-1) */
    double complex zinv = 1.0 / z;
    double complex result = poly->coeff[poly->order];
    for (int k = poly->order - 1; k >= 0; k--) {
        result = result * zinv + poly->coeff[k];
    }
    return result;
}

/*----------------------------------------------------------------------
 * L2 — Z-Rational Function Evaluation
 *----------------------------------------------------------------------*/

double complex z_rational_eval(const ZRational *rat, double complex z)
{
    if (!rat) { errno = EINVAL; return NAN + I * NAN; }

    double complex Nz = z_poly_eval(&rat->numerator, z);
    double complex Dz = z_poly_eval(&rat->denominator, z);

    if (cabs(Dz) < DBL_EPSILON * 1e6) {
        errno = EDOM;
        return NAN + I * NAN;
    }
    return Nz / Dz;
}

/*----------------------------------------------------------------------
 * L4 — Linearity Verification
 *----------------------------------------------------------------------*/

int z_verify_linearity(double a, double b,
                        double complex Xz, double complex Yz,
                        double complex combined,
                        double tolerance)
{
    double complex expected = a * Xz + b * Yz;
    return (cabs(combined - expected) < tolerance) ? 1 : 0;
}

/*----------------------------------------------------------------------
 * L2 — Delay Factor
 *----------------------------------------------------------------------*/

double complex z_delay_factor(double complex z, int k)
{
    /* z^(-k) for delay of k samples */
    return cpow(z, -k);
}

/*----------------------------------------------------------------------
 * L2 — Time Advance Property
 *----------------------------------------------------------------------*/

int z_time_advance(const ZRational *X, int k,
                    const double *initials,
                    ZRational *result)
{
    if (!X || !result || k < 1 || !initials) { errno = EINVAL; return -1; }

    /* Z{x[n+k]} = z^k·X(z) - Σ_{m=0}^{k-1} x[m]·z^{k-m} */
    /* First term: z^k·X(z) = z^k·N(z)/D(z) */
    /* z^k·N(z) = z^k·Σ b_i·z^(-i) = Σ b_i·z^(k-i) */
    /* This increases numerator order in positive powers of z. */
    /* We convert to standard form by multiplying numerator by z^k. */
    /* N(z)·z^k means: shift coefficients right by k positions. */

    ZPoly zk_N;
    zk_N.order = X->numerator.order + k;
    if (zk_N.order > Z_MAX_ORDER) return -1;
    memset(zk_N.coeff, 0, (Z_MAX_ORDER + 1) * sizeof(double));
    /* b_0·z^k + b_1·z^(k-1) + ... + b_n·z^(k-n) = shifted */
    /* Rewrite as z^k·(b_0 + b_1·z^{-1} + ... ) = b_0·z^k + b_1·z^{k-1} + ... */
    /* In z^(-1) representation, this means: coefficients at positive powers */
    /* Standard representation: P(z) = Σ c_i·z^{-i}. Multiplying by z^k: */
    /* P(z)·z^k = Σ c_i·z^{k-i}. The term c_i·z^{k-i} corresponds to z^{-j} */
    /* where j = i - k when i ≥ k, and positive powers of z otherwise. */
    /* For proper representation, we keep standard form and note that */
    /* the result is not strictly proper. We store it anyway. */
    for (int i = 0; i <= X->numerator.order; i++)
        zk_N.coeff[i + k] = X->numerator.coeff[i];

    /* Now subtract Σ_{m=0}^{k-1} x[m]·z^{k-m} * D(z) / D(z) */
    /* This adds another polynomial to the numerator. */
    /* Σ x[m]·z^{k-m} = x[0]·z^k + x[1]·z^{k-1} + ... + x[k-1]·z */
    /* Multiply by D(z): (Σ x[m]·z^{k-m})·D(z) */
    /* D(z) = a_0 + a_1·z^{-1} + ... a_r·z^{-r} */
    /* (x[m]·z^{k-m})·D(z) = x[m]·(a_0·z^{k-m} + a_1·z^{k-m-1} + ...) */

    ZPoly subtract_term;
    memset(&subtract_term, 0, sizeof(subtract_term));
    subtract_term.order = k + X->denominator.order;
    if (subtract_term.order > Z_MAX_ORDER) return -1;

    for (int m = 0; m < k; m++) {
        for (int j = 0; j <= X->denominator.order; j++) {
            int idx = (k - m) + j;  /* power of z: z^{k-m-j} → coefficient at index j but shifted */
            /* x[m]·a_j·z^{k-m}·z^{-j} = x[m]·a_j·z^{k-m-j} */
            /* In standard form c_i·z^{-i}: i = j - (k-m) when positive, else at positive power */
            if (idx >= 0 && idx <= Z_MAX_ORDER)
                subtract_term.coeff[idx] += initials[m] * X->denominator.coeff[j];
        }
    }

    /* Result numerator = zk_N - subtract_term; denominator = D(z) unchanged */
    ZPoly result_num;
    memset(&result_num, 0, sizeof(result_num));
    int max_o = (zk_N.order > subtract_term.order) ? zk_N.order : subtract_term.order;
    result_num.order = max_o;
    for (int i = 0; i <= max_o; i++) {
        double val = 0.0;
        if (i <= zk_N.order) val += zk_N.coeff[i];
        if (i <= subtract_term.order) val -= subtract_term.coeff[i];
        result_num.coeff[i] = val;
    }
    while (result_num.order > 0 && fabs(result_num.coeff[result_num.order]) < DBL_EPSILON)
        result_num.order--;

    result->numerator      = result_num;
    result->denominator    = X->denominator;
    result->sampling_period = X->sampling_period;
    return 0;
}

/*----------------------------------------------------------------------
 * L2 — Frequency Scaling
 *----------------------------------------------------------------------*/

int z_frequency_scaling(const ZRational *X, double a,
                         ZRational *result)
{
    if (!X || !result || fabs(a) < DBL_MIN) { errno = EINVAL; return -1; }

    /* Z{a^n·x[n]} = X(z/a) */
    /* Substitute z → z/a in X(z). Each z^(-k) → a^k·z^(-k) */
    /* So coefficient c_k → c_k · a^k */

    result->numerator.order = X->numerator.order;
    for (int i = 0; i <= X->numerator.order; i++)
        result->numerator.coeff[i] = X->numerator.coeff[i] * pow(a, i);

    result->denominator.order = X->denominator.order;
    for (int i = 0; i <= X->denominator.order; i++)
        result->denominator.coeff[i] = X->denominator.coeff[i] * pow(a, i);

    result->sampling_period = X->sampling_period;
    return 0;
}

/*----------------------------------------------------------------------
 * L4 — Convolution Property
 *----------------------------------------------------------------------*/

int z_convolution(const ZRational *X, const ZRational *Y,
                   ZRational *result)
{
    if (!X || !Y || !result) { errno = EINVAL; return -1; }

    ZPoly num, den;
    if (z_poly_multiply(&X->numerator,   &Y->numerator,   &num) != 0) return -1;
    if (z_poly_multiply(&X->denominator, &Y->denominator, &den) != 0) return -1;

    result->numerator   = num;
    result->denominator = den;
    result->sampling_period = X->sampling_period;
    return 0;
}

/*----------------------------------------------------------------------
 * L4 — Initial Value Theorem: x[0] = lim_{z→∞} X(z)
 *----------------------------------------------------------------------*/

double z_initial_value(const ZRational *X)
{
    if (!X) { errno = EINVAL; return NAN; }

    /* As z→∞, X(z) → b₀/a₀ (coefficients of z^0 terms) */
    /* Because lim_{z→∞} z^{-k} = 0 for k > 0, only constant terms survive. */
    double a0 = X->denominator.coeff[0];
    if (fabs(a0) < DBL_MIN) { errno = EDOM; return NAN; }

    return X->numerator.coeff[0] / a0;
}

/*----------------------------------------------------------------------
 * L4 — Final Value Theorem: x[∞] = lim_{z→1} (z-1)·X(z)
 *----------------------------------------------------------------------*/

double z_final_value(const ZRational *X)
{
    if (!X) { errno = EINVAL; return NAN; }

    /* Must check all poles are strictly inside the unit circle */
    /* (z-1)·X(z) = (z-1)·N(z)/D(z). Limit exists only if all poles satisfy |z|<1 */
    /* We approximate this by evaluating stability: if D(1)=0 cancel the (z-1), */
    /* otherwise check dominant pole magnitude. */

    /* Simplest check: if D(1) ≠ 0, then at z=1: (1-1)·X(1) = 0 */
    double complex z1  = 1.0 + 0.0 * I;
    double complex Dz1 = z_poly_eval(&X->denominator, z1);

    if (cabs(Dz1) > 1e-12) return 0.0;

    /* D(1) = 0: factor (z-1) from denominator */
    /* X(z) = N(z)/((z-1)·D̃(z)), so (z-1)·X(z) = N(z)/D̃(z) */
    /* limit_{z→1} N(z)/D̃(z) = N(1)/D̃(1) */
    /* We need to divide D(z) by (z-1). For polynomial in z^{-1}: */
    /* D(z) = Σ a_i·z^{-i}, D(1)=0 → Σ a_i = 0. Factor out (1 - z^{-1}) */
    /* More directly: D(z) = (1 - z^{-1})·D̃(z). Evaluate D̃ at z=1. */
    /* D̃(z) = D(z)/(1 - z^{-1}) = D(z)·z/(z-1) */
    /* At z→1, use l'Hôpital: D(1)=0, so limit = -dD/d(z^{-1})|_{z=1} */
    /* dD/d(z^{-1})|_{z=1} = Σ i·a_i */

    double dD_at_1 = 0.0;
    for (int i = 1; i <= X->denominator.order; i++)
        dD_at_1 += i * X->denominator.coeff[i];

    if (fabs(dD_at_1) < 1e-12) { errno = EDOM; return NAN; }

    double N_at_1 = 0.0;
    for (int i = 0; i <= X->numerator.order; i++)
        N_at_1 += X->numerator.coeff[i];

    return N_at_1 / (-dD_at_1);
}

/*----------------------------------------------------------------------
 * L3 — Z-Polynomial Arithmetic
 *----------------------------------------------------------------------*/

int z_poly_multiply(const ZPoly *P, const ZPoly *Q, ZPoly *R)
{
    if (!P || !Q || !R) { errno = EINVAL; return -1; }
    int r_order = P->order + Q->order;
    if (r_order > Z_MAX_ORDER) return -1;
    memset(R, 0, sizeof(*R));
    R->order = r_order;
    for (int i = 0; i <= P->order; i++)
        for (int j = 0; j <= Q->order; j++)
            R->coeff[i + j] += P->coeff[i] * Q->coeff[j];
    while (R->order > 0 && fabs(R->coeff[R->order]) < DBL_EPSILON)
        R->order--;
    return 0;
}

int z_poly_add(const ZPoly *P, const ZPoly *Q, ZPoly *R)
{
    if (!P || !Q || !R) { errno = EINVAL; return -1; }
    int max_o = (P->order > Q->order) ? P->order : Q->order;
    memset(R, 0, sizeof(*R));
    R->order = max_o;
    for (int i = 0; i <= P->order; i++) R->coeff[i] += P->coeff[i];
    for (int i = 0; i <= Q->order; i++) R->coeff[i] += Q->coeff[i];
    while (R->order > 0 && fabs(R->coeff[R->order]) < DBL_EPSILON)
        R->order--;
    return 0;
}

/*----------------------------------------------------------------------
 * L5 — Inverse Z-Transform: Power Series (Long Division)
 *----------------------------------------------------------------------*/

int z_inverse_power_series(const ZRational *X,
                            int n_terms,
                            double *x)
{
    if (!X || !x || n_terms < 1) { errno = EINVAL; return -1; }

    /* X(z) = (b₀ + b₁z⁻¹ + ... + bₘz⁻ᵐ) / (a₀ + a₁z⁻¹ + ... + aₙz⁻ⁿ) */
    /* Long division: compute x[0], x[1], ... sequentially */
    /* b₀ + b₁z⁻¹ + ... = (a₀ + a₁z⁻¹ + ...)·(x[0] + x[1]z⁻¹ + ...) */
    /* x[0] = b₀/a₀, then subtract and continue */

    double a0 = X->denominator.coeff[0];
    if (fabs(a0) < DBL_MIN) { errno = EDOM; return -1; }

    /* Working copies */
    double work_b[Z_MAX_ORDER + 1];
    int M = X->numerator.order;
    int N = X->denominator.order;
    for (int i = 0; i <= M; i++) work_b[i] = X->numerator.coeff[i];
    for (int i = M + 1; i <= Z_MAX_ORDER; i++) work_b[i] = 0.0;

    for (int k = 0; k < n_terms; k++) {
        x[k] = work_b[k] / a0;
        /* Subtract x[k] * (a₀ + a₁z⁻¹ + ...)·z^{-k} from the series */
        for (int j = 0; j <= N && (k + j) <= Z_MAX_ORDER; j++) {
            work_b[k + j] -= x[k] * X->denominator.coeff[j];
        }
    }
    return 0;
}

/*----------------------------------------------------------------------
 * L5 — Inverse Z-Transform: Partial Fraction Expansion + Table Lookup
 *----------------------------------------------------------------------*/

int z_inverse_partial_fraction(const ZRational *X,
                                int n_terms,
                                double *x)
{
    if (!X || !x || n_terms < 1) { errno = EINVAL; return -1; }

    /* Decompose X(z) into sum of first-order terms: A_k/(1 - p_k·z^{-1}) */
    /* Each term corresponds to A_k · (p_k)^n · u[n] */
    /* Requires finding denominator roots (poles). */

    /* Build polynomial in z (not z⁻¹) for root-finding */
    /* D(z⁻¹) = Σ a_i·z^{-i}. Multiply by z^N: D̃(z) = Σ a_i·z^{N-i} */
    /* Poles p_k satisfy D(p_k) = 0 → D̃(p_k) = 0 */
    /* We use simple pole extraction for first/second-order systems */

    int N = X->denominator.order;
    int M = X->numerator.order;

    if (N == 0) { errno = EDOM; return -1; }

    /* Convert D(z⁻¹) coefficients to D̃(z) coefficients */
    double d_tilde[Z_MAX_ORDER + 1];
    for (int i = 0; i <= N; i++)
        d_tilde[N - i] = X->denominator.coeff[i];

    /* Simple case: N = 1 → single pole, use power series for reliability */
    if (N == 1) {
        /* First-order: X(z) = b₀/(1 + a₁z⁻¹), x[n] = b₀·(-a₁)^n */
        /* Use power series for general 1st-order case */
        return z_inverse_power_series(X, n_terms, x);
    }

    /* For higher order, use residue method: */
    /* X(z) = Σ_{k} [Residue at p_k of X(z)·z^{n-1}] */
    /* First, factor denominator and compute residues. */
    /* Because polynomial root-finding for z-domain is the same as s-domain, */
    /* we reuse the Laplace pole finder on the reversed polynomial. */

    LaplaceRational temp_L;
    memset(&temp_L, 0, sizeof(temp_L));
    temp_L.denominator.order = N;
    for (int i = 0; i <= N; i++)
        temp_L.denominator.coeff[i] = d_tilde[i];
    temp_L.numerator.order = 0;
    temp_L.numerator.coeff[0] = 1.0;

    LaplacePoleZero pz;
    memset(&pz, 0, sizeof(pz));
    if (laplace_find_poles(&temp_L, &pz) != 0)
        return z_inverse_power_series(X, n_terms, x); /* fallback */

    double complex poles[Z_MAX_ORDER];
    int num_p = pz.num_poles;
    for (int i = 0; i < num_p; i++)
        poles[i] = pz.poles[i];

    /* Compute residues: convert numerator to z (not z⁻¹) form */
    double n_tilde[Z_MAX_ORDER + 1];
    for (int i = 0; i <= M; i++)
        n_tilde[M - i] = X->numerator.coeff[i];
    for (int i = M + 1; i <= N; i++)
        n_tilde[i] = 0.0;

    /* Residue at pole pk for X(z)·z^{n-1}: */
    /* X(z) = Ñ(z)/D̃(z) where Ñ, D̃ are in z (positive powers). */
    /* Residue at pk: Ñ(pk)/D̃'(pk). */
    double A[Z_MAX_ORDER];
    double complex p[Z_MAX_ORDER];

    /* Compute D̃'(z) coefficients */
    double d_tilde_prime[Z_MAX_ORDER + 1];
    for (int i = 1; i <= N; i++)
        d_tilde_prime[i - 1] = i * d_tilde[i];
    /* Horner evaluate D̃'(pk) */
    for (int k = 0; k < num_p; k++) {
        p[k] = poles[k];
        /* Evaluate Ñ(pk) */
        double complex N_val = 0;
        for (int i = 0; i <= M; i++)
            N_val += n_tilde[i] * cpow(p[k], i);
        /* Evaluate D̃'(pk) */
        double complex Dp_val = 0;
        for (int i = 0; i <= N - 1; i++)
            Dp_val += d_tilde_prime[i] * cpow(p[k], i);

        if (cabs(Dp_val) > 1e-12)
            A[k] = creal(N_val / Dp_val);
        else
            A[k] = 0.0;
    }

    /* Build x[n] = Σ A_k · (p_k)^n */
    for (int n = 0; n < n_terms; n++) {
        double xn = 0.0;
        for (int k = 0; k < num_p; k++) {
            xn += A[k] * creal(cpow(p[k], n));
        }
        x[n] = xn;
    }
    return 0;
}

/*----------------------------------------------------------------------
 * L4/L5 — Inverse Z: Residue Method
 *----------------------------------------------------------------------*/

double z_inverse_residue(const ZRational *X, int n)
{
    if (!X || n < 0) { errno = EINVAL; return NAN; }

    /* x[n] = Σ Res{X(z)·z^(n-1)} at poles inside ROC */
    double seq[128];
    int n_terms = (n + 1 > 128) ? 128 : (n + 1);

    if (z_inverse_partial_fraction(X, n_terms, seq) != 0) {
        /* Fallback: long division */
        if (z_inverse_power_series(X, n_terms, seq) != 0)
            return NAN;
    }

    return seq[n];
}

/*----------------------------------------------------------------------
 * L3 — Map Z-Poles back to S-Plane: s = ln(z)/T
 *----------------------------------------------------------------------*/

int z_poles_to_s_plane(const ZPoleZero *zp, double T,
                        LaplacePoleZero *sp)
{
    if (!zp || !sp || T <= 0.0) { errno = EINVAL; return -1; }

    sp->num_poles = zp->num_poles;
    sp->num_zeros = zp->num_zeros;
    sp->gain = zp->gain;

    for (int i = 0; i < zp->num_poles; i++) {
        double complex z = zp->poles[i];
        /* s = ln(z)/T = ln|z|/T + j·arg(z)/T */
        double r = cabs(z);
        double theta = carg(z);
        sp->poles[i] = log(r) / T + (theta / T) * I;
    }
    for (int i = 0; i < zp->num_zeros; i++) {
        double r = cabs(zp->zeros[i]);
        double theta = carg(zp->zeros[i]);
        sp->zeros[i] = log(r) / T + (theta / T) * I;
    }
    return 0;
}
