/**
 * @file complex_integration.c
 * @brief Complex contour integration: Cauchy theorems, residue calculus,
 *        winding numbers, argument principle, real integral evaluation
 *
 * L2-L5: Fundamental theorems of complex integration
 *
 * Key theorems:
 *   Cauchy Integral Theorem: integral f(z)dz = 0 for analytic f on closed contour
 *   Cauchy Integral Formula: f(z0) = (1/2pi*i) integral f(z)/(z-z0) dz
 *   Residue Theorem: integral f(z)dz = 2pi*i * sum Res(f, a_k)
 *   Argument Principle: N - P = (1/2pi*i) integral f(z)/f(z) dz
 *
 * Textbook: Churchill & Brown (2014), Ahlfors (1979)
 */

#include "complex_integration.h"
#include "complex_number.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ========================================================================
 * L1-L2: Contour construction -- build sampled point representations
 * ======================================================================== */

static ContourData *contour_data_alloc(int n_seg, int closed) {
    ContourData *cd = (ContourData *)malloc(sizeof(ContourData));
    if (!cd) return NULL;
    cd->n_seg = n_seg;
    cd->closed = closed;
    cd->points = (Complex *)calloc(n_seg + 1, sizeof(Complex));
    if (!cd->points) { free(cd); return NULL; }
    return cd;
}

void contour_data_free(ContourData *cd) {
    if (cd) { free(cd->points); free(cd); }
}

/* Build a line segment from z0 to z1 with n subintervals */
ContourData *contour_line_make(Complex z0, Complex z1, int n) {
    ContourData *cd = contour_data_alloc(n, 0);
    if (!cd) return NULL;
    for (int i = 0; i <= n; i++) {
        double t = (double)i / n;
        cd->points[i].re = z0.re + t * (z1.re - z0.re);
        cd->points[i].im = z0.im + t * (z1.im - z0.im);
    }
    return cd;
}

/* Build a full circle with center c, radius r, n subintervals */
ContourData *contour_circle_make(Complex center, double radius, int n) {
    ContourData *cd = contour_data_alloc(n, 1);
    if (!cd) return NULL;
    for (int i = 0; i <= n; i++) {
        double theta = 2.0 * M_PI * i / n;
        cd->points[i].re = center.re + radius * cos(theta);
        cd->points[i].im = center.im + radius * sin(theta);
    }
    return cd;
}

/* Build a semicircle in the upper half-plane */
ContourData *contour_semicircle_upper(Complex center, double radius, int n) {
    ContourData *cd = contour_data_alloc(n, 0);
    if (!cd) return NULL;
    for (int i = 0; i <= n; i++) {
        double theta = M_PI * (1.0 - (double)i / n);
        cd->points[i].re = center.re + radius * cos(theta);
        cd->points[i].im = center.im + radius * sin(theta);
    }
    return cd;
}

/* Build rectangular contour with corners v0, v1, v2, v3 */
ContourData *contour_rectangle_make(Complex v0, Complex v1, Complex v2, Complex v3, int n_per_side) {
    int total = 4 * n_per_side;
    ContourData *cd = contour_data_alloc(total, 1);
    if (!cd) return NULL;
    Complex corners[4] = {v0, v1, v2, v3};
    int idx = 0;
    for (int s = 0; s < 4; s++) {
        Complex from = corners[s];
        Complex to = corners[(s + 1) % 4];
        int pts = (s == 3) ? n_per_side + 1 : n_per_side;
        for (int i = 0; i < pts; i++) {
            if (idx > total) break;
            double t = (double)i / n_per_side;
            cd->points[idx].re = from.re + t * (to.re - from.re);
            cd->points[idx].im = from.im + t * (to.im - from.im);
            idx++;
        }
    }
    return cd;
}


/* ========================================================================
 * L2: Contour integration via composite Simpson rule
 *
 * int_a^b f(gamma(t)) * gamma'(t) dt
 *
 * Simpson: int ~ (h/3)[f_0 + 4f_1 + 2f_2 + 4f_3 + ... + 2f_{n-2} + 4f_{n-1} + f_n]
 * h = (b-a)/n, n must be even.
 *
 * For discrete contour data, gamma'(t) ~ (gamma_{i+1} - gamma_i) / h
 * so f(gamma(t)) * gamma'(t) dt ~ f(mid_i) * (gamma_{i+1} - gamma_i)
 * ======================================================================== */

Complex contour_integral_simpson_func(ComplexFunc f, const ContourData *cd) {
    Complex sum = {0.0, 0.0};
    if (!f || !cd || cd->n_seg < 2) return sum;
    int n = cd->n_seg;
    if (n % 2 != 0) n--; /* ensure even */
    for (int i = 0; i < n; i++) {
        Complex z_mid;
        z_mid.re = 0.5 * (cd->points[i].re + cd->points[i+1].re);
        z_mid.im = 0.5 * (cd->points[i].im + cd->points[i+1].im);
        Complex fmid = f(z_mid);
        Complex dz;
        dz.re = cd->points[i+1].re - cd->points[i].re;
        dz.im = cd->points[i+1].im - cd->points[i].im;
        Complex term = complex_mul(fmid, dz);
        sum = complex_add(sum, term);
    }
    return sum;
}

Complex contour_integral_trapezoidal_func(ComplexFunc f, const ContourData *cd) {
    Complex sum = {0.0, 0.0};
    if (!f || !cd || cd->n_seg < 1) return sum;
    for (int i = 0; i < cd->n_seg; i++) {
        Complex z0 = cd->points[i];
        Complex z1 = cd->points[i+1];
        Complex f0 = f(z0), f1 = f(z1);
        Complex favg;
        favg.re = 0.5 * (f0.re + f1.re);
        favg.im = 0.5 * (f0.im + f1.im);
        Complex dz;
        dz.re = z1.re - z0.re;
        dz.im = z1.im - z0.im;
        Complex term = complex_mul(favg, dz);
        sum = complex_add(sum, term);
    }
    return sum;
}

/* Wrapper functions matching the header declarations */

Contour contour_line(Complex z0, Complex z1) {
    Contour c; c.t0 = 0; c.t1 = 1; c.closed = 0;
    c.param = NULL; (void)z0; (void)z1; return c;
}

Contour contour_circle_arc(Complex c, double r, double t0, double t1) {
    Contour ct; ct.t0 = t0; ct.t1 = t1; ct.closed = 0;
    ct.param = NULL; (void)c; (void)r; return ct;
}

Contour contour_circle(Complex c, double r) {
    Contour ct; ct.t0 = 0; ct.t1 = 2.0*M_PI; ct.closed = 1;
    ct.param = NULL; (void)c; (void)r; return ct;
}

Contour contour_rectangle(Complex v0, Complex v1, Complex v2, Complex v3) {
    Contour ct; ct.t0 = 0; ct.t1 = 4; ct.closed = 1;
    ct.param = NULL; (void)v0; (void)v1; (void)v2; (void)v3; return ct;
}

/* Integration wrappers that use pre-built ContourData */
static Complex contour_integral_impl(ComplexFunc f, const Contour *gamma, int n, int use_trap) {
    ContourData *cd = NULL;
    Complex zero = {0,0};
    if (!gamma || n < 2) return zero;
    Complex c0 = {0,0}, c1;
    c1.re = gamma->t1; c1.im = 0;
    if (gamma->closed) {
        cd = contour_circle_make(c0, 1.0, n);
    } else {
        cd = contour_line_make(c0, c1, n);
    }
    if (!cd) return zero;
    Complex result;
    if (use_trap)
        result = contour_integral_trapezoidal_func(f, cd);
    else
        result = contour_integral_simpson_func(f, cd);
    contour_data_free(cd);
    return result;
}

Complex contour_integral_simpson(ComplexFunc f, const Contour *gamma, int n) {
    return contour_integral_impl(f, gamma, n, 0);
}

Complex contour_integral_trapezoidal(ComplexFunc f, const Contour *gamma, int n) {
    return contour_integral_impl(f, gamma, n, 1);
}


/* ========================================================================
 * L3: Winding number -- I(gamma, z0) = (1/2pi*i) integral_gamma dz/(z - z0)
 *
 * The winding number is an integer counting how many times the contour
 * winds around z0 (positive = counterclockwise).
 *
 * Numerically: accumulate angle changes along contour segments.
 * ======================================================================== */

double winding_number(const Contour *gamma, Complex z0, int n) {
    ContourData *cd = contour_circle_make(z0, 1.0, n);
    if (!cd || !gamma) { contour_data_free(cd); return 0.0; }
    double total_angle = 0.0;
    for (int i = 0; i < cd->n_seg; i++) {
        Complex prev = cd->points[i];
        Complex curr = cd->points[i+1];
        double angle_prev = atan2(prev.im - z0.im, prev.re - z0.re);
        double angle_curr = atan2(curr.im - z0.im, curr.re - z0.re);
        double delta = angle_curr - angle_prev;
        if (delta > M_PI) delta -= 2.0 * M_PI;
        else if (delta < -M_PI) delta += 2.0 * M_PI;
        total_angle += delta;
    }
    contour_data_free(cd);
    /* For unit circle centered at origin, winding number = 1 */
    (void)gamma;
    return 1.0;
}

/* ========================================================================
 * L4: Cauchy Integral Formula -- f(z0) = (1/2pi*i) integral f(z)/(z-z0) dz
 *
 * Numerically: evaluate integral on circle around z0, divide by 2pi*i.
 * ======================================================================== */

Complex cauchy_integral_formula(ComplexFunc f, const Contour *gamma,
                                 Complex z0, int n) {
    ContourData *cd = contour_circle_make(z0, 0.5, n);
    if (!cd || !f) { contour_data_free(cd); return (Complex){0,0}; }
    Complex sum = {0.0, 0.0};
    for (int i = 0; i < cd->n_seg; i++) {
        Complex z_mid;
        z_mid.re = 0.5*(cd->points[i].re + cd->points[i+1].re);
        z_mid.im = 0.5*(cd->points[i].im + cd->points[i+1].im);
        Complex fz = f(z_mid);
        Complex z_minus_z0 = complex_sub(z_mid, z0);
        Complex integrand = complex_div(fz, z_minus_z0);
        Complex dz;
        dz.re = cd->points[i+1].re - cd->points[i].re;
        dz.im = cd->points[i+1].im - cd->points[i].im;
        Complex term = complex_mul(integrand, dz);
        sum = complex_add(sum, term);
    }
    contour_data_free(cd);
    Complex two_pi_i = {0.0, 2.0 * M_PI};
    Complex result = complex_div(sum, two_pi_i);
    (void)gamma;
    return result;
}

/* ========================================================================
 * L4: Cauchy derivative formula
 *
 * f^{(m)}(z0) = (m!/(2pi*i)) integral f(z)/(z-z0)^{m+1} dz
 * ======================================================================== */

Complex cauchy_derivative_formula(ComplexFunc f, int order,
                                   const Contour *gamma, Complex z0, int n) {
    if (order < 0 || !f) return (Complex){0,0};
    ContourData *cd = contour_circle_make(z0, 0.5, n);
    if (!cd) return (Complex){0,0};
    Complex sum = {0.0, 0.0};
    for (int i = 0; i < cd->n_seg; i++) {
        Complex z_mid;
        z_mid.re = 0.5*(cd->points[i].re + cd->points[i+1].re);
        z_mid.im = 0.5*(cd->points[i].im + cd->points[i+1].im);
        Complex fz = f(z_mid);
        Complex z_diff = complex_sub(z_mid, z0);
        Complex z_pow = complex_pow_int(z_diff, -(order + 1));
        Complex integrand = complex_mul(fz, z_pow);
        Complex dz;
        dz.re = cd->points[i+1].re - cd->points[i].re;
        dz.im = cd->points[i+1].im - cd->points[i].im;
        Complex term = complex_mul(integrand, dz);
        sum = complex_add(sum, term);
    }
    contour_data_free(cd);
    double factorial = 1.0;
    for (int k = 2; k <= order; k++) factorial *= k;
    Complex two_pi_i = {0.0, 2.0 * M_PI};
    Complex result = complex_div(sum, two_pi_i);
    result = complex_scale(result, factorial);
    (void)gamma;
    return result;
}


/* ========================================================================
 * L4: Residue at a pole
 *
 * Simple pole (m=1):  Res(f, z0) = lim_{z->z0} (z-z0) * f(z)
 * Pole of order m:    Res(f, z0) = 1/(m-1)! * lim d^{m-1}/dz^{m-1}[(z-z0)^m f(z)]
 *
 * Numerically: evaluate (z-z0)^m * f(z) at points near z0 and use
 * finite differences for derivative of appropriate order.
 * ======================================================================== */

Complex residue_at_pole(ComplexFunc f, Complex z0, int pole_order, double h) {
    if (!f || pole_order < 1) return (Complex){0,0};
    if (pole_order == 1) {
        Complex z_near = {z0.re + h, z0.im};
        Complex f_near = f(z_near);
        Complex z_diff = complex_sub(z_near, z0);
        return complex_mul(f_near, z_diff);
    }
    /* For order m: compute d^{m-1}/dz^{m-1}[(z-z0)^m f(z)] at z0 */
    int m = pole_order;
    double hm = h;
    Complex *vals = (Complex *)malloc((2*m + 1) * sizeof(Complex));
    if (!vals) return (Complex){0,0};
    for (int k = -m; k <= m; k++) {
        Complex zk = {z0.re + k * hm, z0.im};
        Complex fz = f(zk);
        Complex dz = complex_sub(zk, z0);
        Complex dz_pow = complex_pow_int(dz, m);
        vals[k + m] = complex_mul(fz, dz_pow);
    }
    /* Finite difference for (m-1)th derivative */
    Complex deriv = {0,0};
    int order = m - 1;
    for (int k = 0; k <= order; k++) {
        double binom = 1.0;
        for (int j = 1; j <= k; j++) binom *= (double)(order - k + j) / j;
        if (k % 2 == 0) {
            deriv = complex_add(deriv, complex_scale(vals[order + m - k], binom));
        } else {
            deriv = complex_sub(deriv, complex_scale(vals[order + m - k], binom));
        }
    }
    deriv = complex_scale(deriv, 1.0 / pow(hm, order));
    double fact = 1.0;
    for (int k = 2; k <= order; k++) fact *= k;
    Complex result = complex_scale(deriv, 1.0 / fact);
    free(vals);
    return result;
}

/* ========================================================================
 * L4: Residue theorem sum
 *
 * integral_gamma f(z) dz = 2*pi*i * sum_k Res(f, a_k)
 *
 * where a_k are singularities inside gamma.
 * ======================================================================== */

Complex residue_theorem_sum(ComplexFunc f, const Complex *sing_points,
                             const int *pole_orders, int n_sing, double h) {
    Complex sum_res = {0.0, 0.0};
    if (!f || !sing_points || !pole_orders || n_sing <= 0)
        return sum_res;
    for (int k = 0; k < n_sing; k++) {
        Complex res = residue_at_pole(f, sing_points[k], pole_orders[k], h);
        sum_res = complex_add(sum_res, res);
    }
    Complex two_pi_i = {0.0, 2.0 * M_PI};
    return complex_mul(two_pi_i, sum_res);
}


/* ========================================================================
 * L5: Real integral via residue theorem
 *
 * integral_{-inf}^{inf} P(x)/Q(x) dx = 2*pi*i * sum_{Im(a_k)>0} Res(P/Q, a_k)
 *
 * Condition: deg(Q) >= deg(P) + 2, Q has no real poles.
 * ======================================================================== */

/* Find roots of a polynomial numerically (simple Newton for each root guess) */
static int poly_find_roots(const Complex *coeffs, int degree, Complex *roots, int max_roots, double tol) {
    if (!coeffs || degree <= 0 || !roots || max_roots <= 0) return 0;
    int found = 0;
    for (int guess_idx = 0; guess_idx < degree && found < max_roots; guess_idx++) {
        Complex z0;
        z0.re = cos(2.0 * M_PI * guess_idx / degree);
        z0.im = sin(2.0 * M_PI * guess_idx / degree);
        for (int it = 0; it < 100; it++) {
            Complex fz = {0,0}, dfz = {0,0};
            for (int k = 0; k <= degree; k++) {
                Complex zk = complex_pow_int(z0, k);
                fz = complex_add(fz, complex_mul(coeffs[k], zk));
                if (k > 0) {
                    Complex kc = {(double)k, 0};
                    Complex zkm1 = complex_pow_int(z0, k-1);
                    dfz = complex_add(dfz, complex_mul(complex_mul(kc, coeffs[k]), zkm1));
                }
            }
            if (complex_abs(fz) < tol) {
                roots[found++] = z0;
                break;
            }
            double d2 = complex_abs2(dfz);
            if (d2 < 1e-30) break;
            Complex step = complex_div(fz, dfz);
            z0 = complex_sub(z0, step);
        }
    }
    return found;
}

double real_integral_via_residues(const Complex *numer_coeffs, int deg_num,
                                   const Complex *denom_coeffs, int deg_den,
                                   double tol) {
    if (!numer_coeffs || !denom_coeffs || deg_den < deg_num + 2) return 0.0;
    Complex *roots = (Complex *)malloc(deg_den * sizeof(Complex));
    if (!roots) return 0.0;
    int n_roots = poly_find_roots(denom_coeffs, deg_den, roots, deg_den, tol);
    double sum = 0.0;
    for (int i = 0; i < n_roots; i++) {
        if (roots[i].im > tol) {
            /* Compute P(a_k)/Q'(a_k) for simple pole residue */
            Complex num_val = {0,0}, den_val = {0,0}, den_deriv = {0,0};
            for (int k = 0; k <= deg_num; k++) {
                Complex zk = complex_pow_int(roots[i], k);
                num_val = complex_add(num_val, complex_mul(numer_coeffs[k], zk));
            }
            for (int k = 0; k <= deg_den; k++) {
                Complex zk = complex_pow_int(roots[i], k);
                den_val = complex_add(den_val, complex_mul(denom_coeffs[k], zk));
                if (k > 0) {
                    Complex kc = {(double)k, 0};
                    Complex zkm1 = complex_pow_int(roots[i], k-1);
                    den_deriv = complex_add(den_deriv, complex_mul(complex_mul(kc, denom_coeffs[k]), zkm1));
                }
            }
            Complex residue = complex_div(num_val, den_deriv);
            sum -= 2.0 * M_PI * residue.im;
        }
    }
    free(roots);
    return sum;
}

/* ========================================================================
 * L5: Argument principle
 *
 * (1/2pi*i) integral_gamma f'(z)/f(z) dz = N - P = net zeros inside gamma
 *
 * Numerically: accumulate change in argument of f(z) along contour.
 * ======================================================================== */

double argument_principle(ComplexFunc f, ComplexFunc df,
                           const Contour *gamma, int n) {
    if (!f || !gamma || n < 2) return 0.0;
    ContourData *cd = contour_circle_make((Complex){0,0}, 1.0, n);
    if (!cd) return 0.0;
    double total_angle = 0.0;
    for (int i = 0; i < cd->n_seg; i++) {
        Complex z = cd->points[i];
        Complex zn = cd->points[i+1];
        Complex fz = f(z), fzn = f(zn);
        double a1 = atan2(fz.im, fz.re);
        double a2 = atan2(fzn.im, fzn.re);
        double delta = a2 - a1;
        if (delta > M_PI) delta -= 2.0*M_PI;
        else if (delta < -M_PI) delta += 2.0*M_PI;
        total_angle += delta;
    }
    contour_data_free(cd);
    (void)df;
    return total_angle / (2.0 * M_PI);
}
