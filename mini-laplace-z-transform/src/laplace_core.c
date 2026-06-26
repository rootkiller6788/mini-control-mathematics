/**
 * @file laplace_core.c
 * @brief Laplace Transform — Core implementations
 *
 * Implements: polynomial arithmetic, rational evaluation, Laplace properties,
 *             Initial/Final Value theorems, numerical Laplace and inverse.
 *
 * Each function represents an independent knowledge point (no filler).
 * Standard C99 — no compiler extensions.
 */

#include "laplace_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <float.h>

/*----------------------------------------------------------------------
 * Static Helpers (no GCC extensions — plain C99)
 *----------------------------------------------------------------------*/

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Custom cube root for portability */
static double my_cbrt(double x) {
    if (x >= 0.0) return pow(x, 1.0 / 3.0);
    else return -pow(-x, 1.0 / 3.0);
}

/* Linear interpolation helper for sampled signals */
static double signal_get_sample(const TimeSignal *sig, double t)
{
    if (!sig || sig->n_samples < 1) return 0.0;
    double t_end = sig->t_start + (sig->n_samples - 1) * sig->t_step;
    if (t < sig->t_start || t > t_end) return 0.0;
    int idx = (int)((t - sig->t_start) / sig->t_step + 0.5);
    if (idx >= sig->n_samples) idx = sig->n_samples - 1;
    if (idx < 0) idx = 0;
    return sig->values[idx];
}

/* Polynomial binomial shift: substitute s → s - shift_b in polynomial P */
static int poly_binomial_shift(const LaplacePolynomial *P,
                                double shift_b,
                                LaplacePolynomial *Q)
{
    memset(Q, 0, sizeof(*Q));
    Q->order = P->order;
    if (Q->order > LAPLACE_MAX_ORDER) return -1;

    for (int k = 0; k <= P->order; k++) {
        if (P->coeff[k] == 0.0) continue;
        /* Expand (s - shift_b)^k = Σ_{j=0}^k C(k,j) · s^j · (-shift_b)^{k-j} */
        double binom = 1.0;         /* C(k,0) = 1 */
        double b_pow = 1.0;         /* (-shift_b)^k */
        for (int m = 0; m < k; m++) b_pow *= (-shift_b);

        for (int j = 0; j <= k; j++) {
            Q->coeff[j] += P->coeff[k] * binom * b_pow;
            if (j < k) {
                binom = binom * (k - j) / (j + 1.0);
                /* b_pow goes from (-shift_b)^{k-j} to (-shift_b)^{k-j-1} */
                if (fabs(shift_b) > 1e-30) b_pow /= (-shift_b);
                else b_pow = 0.0;
            }
        }
    }
    while (Q->order > 0 && fabs(Q->coeff[Q->order]) < DBL_EPSILON)
        Q->order--;
    return 0;
}

/* Quadratic formula */
static int solve_quadratic(double a, double b, double c,
                            double complex roots[2])
{
    if (fabs(a) < DBL_EPSILON) {
        if (fabs(b) < DBL_EPSILON) return 0;
        roots[0] = -c / b;
        return 1;
    }
    double disc = b * b - 4.0 * a * c;
    if (disc >= 0) {
        double sqrt_d = sqrt(disc);
        roots[0] = (-b + sqrt_d) / (2.0 * a);
        roots[1] = (-b - sqrt_d) / (2.0 * a);
    } else {
        double re = -b / (2.0 * a);
        double im = sqrt(-disc) / (2.0 * a);
        roots[0] = re + im * I;
        roots[1] = re - im * I;
    }
    return 2;
}

/* Cubic: Cardano — x³ + a·x² + b·x + c = 0 */
static int solve_cubic(double a, double b, double c,
                        double complex roots[3])
{
    double p = b - a * a / 3.0;
    double q = c - a * b / 3.0 + 2.0 * a * a * a / 27.0;
    double disc = (q * q / 4.0) + (p * p * p / 27.0);

    double complex y[3];
    if (disc > 0.0) {
        double sqrt_d = sqrt(disc);
        double u = my_cbrt(-q / 2.0 + sqrt_d);
        double v = my_cbrt(-q / 2.0 - sqrt_d);
        y[0] = u + v;
        y[1] = -(u + v) / 2.0 + (u - v) * sqrt(3.0) / 2.0 * I;
        y[2] = -(u + v) / 2.0 - (u - v) * sqrt(3.0) / 2.0 * I;
    } else if (disc < 0.0) {
        double r = sqrt(-p * p * p / 27.0);
        double phi = acos(-q / (2.0 * r));
        double two_cbrt = 2.0 * sqrt(-p / 3.0);
        y[0] = two_cbrt * cos(phi / 3.0);
        y[1] = two_cbrt * cos((phi + 2.0 * M_PI) / 3.0);
        y[2] = two_cbrt * cos((phi + 4.0 * M_PI) / 3.0);
    } else {
        double u = my_cbrt(-q / 2.0);
        y[0] = 2.0 * u;
        y[1] = -u;
        y[2] = -u;
    }

    double shift = a / 3.0;
    for (int i = 0; i < 3; i++) {
        roots[i] = y[i] - shift;
    }
    return 3;
}

/*======================================================================
 * L2 — Polynomial Evaluation (Horner's Method)
 *======================================================================*/

double complex laplace_poly_eval(const LaplacePolynomial *poly,
                                  double complex s)
{
    if (!poly || poly->order < 0 || poly->order > LAPLACE_MAX_ORDER) {
        errno = EINVAL;
        return NAN + I * NAN;
    }
    double complex result = poly->coeff[poly->order];
    for (int k = poly->order - 1; k >= 0; k--) {
        result = result * s + poly->coeff[k];
    }
    return result;
}

/*======================================================================
 * L2 — Rational Function Evaluation
 *======================================================================*/

double complex laplace_rational_eval(const LaplaceRational *rat,
                                      double complex s)
{
    if (!rat) { errno = EINVAL; return NAN + I * NAN; }
    double complex num_val = laplace_poly_eval(&rat->numerator, s);
    double complex den_val = laplace_poly_eval(&rat->denominator, s);

    if (cabs(den_val) < DBL_EPSILON * 1e6) {
        errno = EDOM;
        return NAN + I * NAN;
    }
    return num_val / den_val;
}

/*======================================================================
 * L4 — Linearity Verification
 *======================================================================*/

int laplace_verify_linearity(double a, double b,
                              double complex Fs, double complex Gs,
                              double complex combined,
                              double tolerance)
{
    double complex expected = a * Fs + b * Gs;
    return (cabs(combined - expected) < tolerance) ? 1 : 0;
}

/*======================================================================
 * L2 — Time Shift Factor
 *======================================================================*/

double complex laplace_time_shift_factor(double complex s, double tau)
{
    return cexp(-s * tau);
}

/*======================================================================
 * L2 — Frequency Shift
 *======================================================================*/

int laplace_frequency_shift(const LaplaceRational *F,
                             double a,
                             LaplaceRational *result)
{
    if (!F || !result) { errno = EINVAL; return -1; }

    /* s → s+a = s-(-a), so shift_b = -a */
    double shift_b = -a;
    LaplacePolynomial temp_num, temp_den;
    if (poly_binomial_shift(&F->numerator,   shift_b, &temp_num) != 0) return -1;
    if (poly_binomial_shift(&F->denominator, shift_b, &temp_den) != 0) return -1;

    result->numerator   = temp_num;
    result->denominator = temp_den;
    return 0;
}

/*======================================================================
 * L4 — Initial Value Theorem: f(0⁺) = lim_{s→∞} s·F(s)
 *======================================================================*/

double laplace_initial_value(const LaplaceRational *F)
{
    if (!F) { errno = EINVAL; return NAN; }
    int nN = F->numerator.order;
    int nD = F->denominator.order;
    if (nN >= nD) { errno = EDOM; return NAN; }

    double a_n = F->numerator.coeff[nN];
    double b_m = F->denominator.coeff[nD];
    if (fabs(b_m) < DBL_MIN) { errno = EDOM; return NAN; }

    /* s·F(s) ~ s·(a_n·s^nN)/(b_m·s^nD) = (a_n/b_m)·s^{1+nN-nD} */
    /* nN < nD: if nN == nD-1 → limit = a_n/b_m; else → 0 */
    return (nN == nD - 1) ? (a_n / b_m) : 0.0;
}

/*======================================================================
 * L4 — Final Value Theorem: f(∞) = lim_{s→0} s·F(s)
 *======================================================================*/

double laplace_final_value(const LaplaceRational *F)
{
    if (!F) { errno = EINVAL; return NAN; }

    LaplacePoleZero pz;
    if (laplace_find_poles(F, &pz) != 0) return NAN;

    for (int i = 0; i < pz.num_poles; i++) {
        double re = creal(pz.poles[i]);
        if (re > 1e-12) { errno = EDOM; return NAN; }         /* RHP pole */
        if (fabs(re) < 1e-12 && fabs(cimag(pz.poles[i])) > 1e-12)
            { errno = EDOM; return NAN; }                      /* jω-axis oscillatory */
    }

    double complex s0 = 0.0;
    double complex D0 = laplace_poly_eval(&F->denominator, s0);
    double complex N0 = laplace_poly_eval(&F->numerator, s0);

    if (cabs(D0) > 1e-12) return 0.0;  /* D(0)≠0 → s·F(s)|_{s=0}=0 */

    /* D(0)=0: use l'Hôpital — limit = N(0)/D'(0) */
    LaplacePolynomial dD;
    laplace_poly_derivative(&F->denominator, &dD);
    double complex Dp0 = laplace_poly_eval(&dD, s0);
    if (cabs(Dp0) < 1e-12) { errno = EDOM; return NAN; }

    return creal(N0) / creal(Dp0);
}

/*======================================================================
 * L2 — Differentiation Property
 *======================================================================*/

int laplace_differentiation(const LaplaceRational *F,
                             double f0,
                             LaplaceRational *dF)
{
    if (!F || !dF) { errno = EINVAL; return -1; }

    LaplacePolynomial sN;
    sN.order = F->numerator.order + 1;
    if (sN.order > LAPLACE_MAX_ORDER) return -1;
    sN.coeff[0] = 0.0;
    for (int i = 0; i <= F->numerator.order; i++)
        sN.coeff[i + 1] = F->numerator.coeff[i];

    LaplacePolynomial f0D;
    f0D.order = F->denominator.order;
    for (int i = 0; i <= F->denominator.order; i++)
        f0D.coeff[i] = f0 * F->denominator.coeff[i];

    int max_o = (sN.order > f0D.order) ? sN.order : f0D.order;
    LaplacePolynomial num;
    memset(&num, 0, sizeof(num));
    num.order = max_o;
    for (int i = 0; i <= sN.order; i++) num.coeff[i] += sN.coeff[i];
    for (int i = 0; i <= f0D.order; i++) num.coeff[i] -= f0D.coeff[i];

    while (num.order > 0 && fabs(num.coeff[num.order]) < DBL_EPSILON)
        num.order--;

    dF->numerator   = num;
    dF->denominator = F->denominator;
    return 0;
}

/*======================================================================
 * L2 — Integration Property
 *======================================================================*/

int laplace_integration(const LaplaceRational *F,
                         LaplaceRational *integral)
{
    if (!F || !integral) { errno = EINVAL; return -1; }
    integral->numerator = F->numerator;

    LaplacePolynomial sD;
    sD.order = F->denominator.order + 1;
    if (sD.order > LAPLACE_MAX_ORDER) return -1;
    sD.coeff[0] = 0.0;
    for (int i = 0; i <= F->denominator.order; i++)
        sD.coeff[i + 1] = F->denominator.coeff[i];
    integral->denominator = sD;
    return 0;
}

/*======================================================================
 * L2/L4 — Convolution Theorem
 *======================================================================*/

int laplace_convolution(const LaplaceRational *F,
                         const LaplaceRational *G,
                         LaplaceRational *conv)
{
    if (!F || !G || !conv) { errno = EINVAL; return -1; }
    LaplacePolynomial num, den;
    if (laplace_poly_multiply(&F->numerator,   &G->numerator,   &num) != 0) return -1;
    if (laplace_poly_multiply(&F->denominator, &G->denominator, &den) != 0) return -1;
    conv->numerator   = num;
    conv->denominator = den;
    return 0;
}

/*======================================================================
 * L3 — Polynomial Multiplication (Convolution of Coefficients)
 *======================================================================*/

int laplace_poly_multiply(const LaplacePolynomial *P,
                           const LaplacePolynomial *Q,
                           LaplacePolynomial *R)
{
    if (!P || !Q || !R) { errno = EINVAL; return -1; }
    int r_order = P->order + Q->order;
    if (r_order > LAPLACE_MAX_ORDER) return -1;

    memset(R, 0, sizeof(*R));
    R->order = r_order;

    for (int i = 0; i <= P->order; i++) {
        if (P->coeff[i] == 0.0) continue;
        for (int j = 0; j <= Q->order; j++) {
            R->coeff[i + j] += P->coeff[i] * Q->coeff[j];
        }
    }
    while (R->order > 0 && fabs(R->coeff[R->order]) < DBL_EPSILON)
        R->order--;
    return 0;
}

/*======================================================================
 * L3 — Polynomial Addition
 *======================================================================*/

int laplace_poly_add(const LaplacePolynomial *P,
                      const LaplacePolynomial *Q,
                      LaplacePolynomial *R)
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

/*======================================================================
 * L3 — Polynomial Derivative
 *======================================================================*/

int laplace_poly_derivative(const LaplacePolynomial *P,
                             LaplacePolynomial *dP)
{
    if (!P || !dP) { errno = EINVAL; return -1; }
    if (P->order <= 0) { dP->order = 0; dP->coeff[0] = 0.0; return 0; }
    dP->order = P->order - 1;
    for (int i = 0; i <= dP->order; i++)
        dP->coeff[i] = (i + 1.0) * P->coeff[i + 1];
    return 0;
}

/*======================================================================
 * L3/L5 — Pole Finding (Quadratic/Cubic/Quartic + Aberth-Ehrlich)
 *======================================================================*/

int laplace_find_poles(const LaplaceRational *F,
                        LaplacePoleZero *pz)
{
    if (!F || !pz) { errno = EINVAL; return -1; }

    const LaplacePolynomial *D = &F->denominator;
    int n = D->order;
    if (n <= 0) return -1;

    double lead = D->coeff[n];
    if (fabs(lead) < DBL_MIN) return -1;

    pz->num_poles = 0;
    pz->gain = 1.0;

    if (n == 1) {
        pz->poles[0]  = -D->coeff[0] / lead;
        pz->num_poles = 1;
    } else if (n == 2) {
        pz->num_poles = solve_quadratic(lead, D->coeff[1], D->coeff[0], pz->poles);
    } else if (n == 3) {
        double a = D->coeff[2] / lead;
        double b = D->coeff[1] / lead;
        double c = D->coeff[0] / lead;
        pz->num_poles = solve_cubic(a, b, c, pz->poles);
    } else {
        /* Higher order: Aberth-Ehrlich simultaneous iteration */
        double R = 0.0;
        for (int i = 0; i <= n; i++) {
            double abs_c = fabs(D->coeff[i] / lead);
            if (abs_c > R) R = abs_c;
        }
        R += 1.0;

        double complex z[LAPLACE_MAX_ORDER];
        for (int i = 0; i < n; i++) {
            double theta = 2.0 * M_PI * i / n;
            z[i] = R * cos(theta) + R * sin(theta) * I;
        }

        LaplacePolynomial dP;
        laplace_poly_derivative(D, &dP);

        for (int iter = 0; iter < 50; iter++) {
            int converged = 1;
            for (int i = 0; i < n; i++) {
                double complex Pz  = laplace_poly_eval(D, z[i]);
                double complex Ppz = laplace_poly_eval(&dP, z[i]);
                if (cabs(Ppz) < 1e-14) { converged = 0; continue; }

                double complex ratio = Pz / Ppz;
                double complex sum_ab = 0;
                for (int j = 0; j < n; j++) {
                    if (j == i) continue;
                    double complex diff = z[i] - z[j];
                    if (cabs(diff) > 1e-14) sum_ab += 1.0 / diff;
                }

                double complex dz = ratio / (1.0 - ratio * sum_ab);
                z[i] -= dz;
                if (cabs(dz) > 1e-10) converged = 0;
            }
            if (converged) break;
        }

        pz->num_poles = n;
        for (int i = 0; i < n; i++) pz->poles[i] = z[i];
    }
    return 0;
}

/*======================================================================
 * L3 — Zero Finding
 *======================================================================*/

int laplace_find_zeros(const LaplaceRational *F,
                        LaplacePoleZero *pz)
{
    if (!F || !pz) { errno = EINVAL; return -1; }

    LaplaceRational dummy;
    memset(&dummy, 0, sizeof(dummy));
    /* Assign denominator = numerator polynomial to reuse pole finder */
    dummy.denominator = F->numerator;
    dummy.numerator.order = 0;
    dummy.numerator.coeff[0] = 1.0;

    LaplacePoleZero zp;
    memset(&zp, 0, sizeof(zp));
    int ret = laplace_find_poles(&dummy, &zp);
    if (ret != 0) return ret;

    pz->num_zeros = zp.num_poles;
    for (int i = 0; i < zp.num_poles; i++)
        pz->zeros[i] = zp.poles[i];
    return 0;
}

/*======================================================================
 * L3 — Numerical Laplace Transform (Composite Simpson's Rule)
 *======================================================================*/

double complex laplace_numerical_transform(const TimeSignal *sig,
                                            double complex s,
                                            double T_max)
{
    if (!sig || sig->n_samples < 2) { errno = EINVAL; return NAN + I * NAN; }
    if (T_max <= 0.0) T_max = sig->t_start + sig->n_samples * sig->t_step;

    double dt_orig = sig->t_step;
    int N = (int)((T_max - sig->t_start) / dt_orig);
    if (N % 2 != 0) N++;
    if (N < 4) N = 4;

    double h = (T_max - sig->t_start) / N;
    double complex integral = 0;

    double t = sig->t_start;
    integral += signal_get_sample(sig, t) * cexp(-s * t);

    for (int i = 1; i < N; i++) {
        t = sig->t_start + i * h;
        double w = (i % 2 == 0) ? 2.0 : 4.0;
        integral += w * signal_get_sample(sig, t) * cexp(-s * t);
    }

    t = sig->t_start + N * h;
    integral += signal_get_sample(sig, t) * cexp(-s * t);

    return integral * (h / 3.0);
}

/*======================================================================
 * L4 — Numerical Inverse Laplace (Dubner-Abate-Crump method)
 *======================================================================*/

double laplace_inverse_numerical(const LaplaceRational *F,
                                  double t,
                                  int n_terms)
{
    if (!F || t < 0.0 || n_terms < 2) { errno = EINVAL; return NAN; }
    if (t == 0.0) return laplace_initial_value(F);

    double sigma  = 0.1;
    double Tper   = 2.0 * t;
    if (Tper < 1.0) Tper = 1.0;

    double sum = 0.5 * creal(laplace_rational_eval(F, sigma + 0.0 * I));

    for (int k = 1; k <= n_terms; k++) {
        double omega_k = k * M_PI / Tper;
        double complex s_k = sigma + omega_k * I;
        double complex Fk  = laplace_rational_eval(F, s_k);

        /* Skip NaN from pole crossings */
        if (isnan(creal(Fk)) || isnan(cimag(Fk))) continue;

        double term = creal(Fk) * cos(omega_k * t);

        /* Lanczos sigma factor for convergence acceleration */
        double sigma_fac = 1.0;
        if (n_terms > 1) {
            double x = M_PI * k / (n_terms + 1.0);
            sigma_fac = sin(x) / x;
        }
        sum += term * sigma_fac;
    }

    return (2.0 * exp(sigma * t) / Tper) * sum;
}
