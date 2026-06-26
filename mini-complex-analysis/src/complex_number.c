/**
 * @file complex_number.c
 * @brief Complex number arithmetic, transcendental functions,
 *        trig/hyperbolic, polynomials, field axioms (L1-L4)
 */

#include "complex_number.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

Complex complex_add(Complex z1, Complex z2) {
    Complex r = {z1.re + z2.re, z1.im + z2.im};
    return r;
}

Complex complex_sub(Complex z1, Complex z2) {
    Complex r = {z1.re - z2.re, z1.im - z2.im};
    return r;
}

Complex complex_mul(Complex z1, Complex z2) {
    Complex r = {z1.re*z2.re - z1.im*z2.im, z1.re*z2.im + z1.im*z2.re};
    return r;
}

Complex complex_div(Complex z1, Complex z2) {
    double d = z2.re*z2.re + z2.im*z2.im;
    if (d == 0.0) { Complex r = {INFINITY, INFINITY}; return r; }
    Complex r = {(z1.re*z2.re + z1.im*z2.im)/d, (z1.im*z2.re - z1.re*z2.im)/d};
    return r;
}
Complex complex_neg(Complex z) {
    Complex r = {-z.re, -z.im};
    return r;
}

Complex complex_conj(Complex z) {
    Complex r = {z.re, -z.im};
    return r;
}

double complex_abs(Complex z) {
    return hypot(z.re, z.im);
}

double complex_abs2(Complex z) {
    return z.re * z.re + z.im * z.im;
}

double complex_arg(Complex z) {
    return atan2(z.im, z.re);
}

Complex complex_polar(double r, double theta) {
    Complex z = {r * cos(theta), r * sin(theta)};
    return z;
}

double complex_radius(Complex z) {
    return hypot(z.re, z.im);
}

double complex_angle(Complex z) {
    return atan2(z.im, z.re);
}

Complex complex_exp(Complex z) {
    double e_x = exp(z.re);
    Complex r = {e_x * cos(z.im), e_x * sin(z.im)};
    return r;
}

Complex complex_log(Complex z) {
    if (z.re == 0.0 && z.im == 0.0) {
        Complex r = {-INFINITY, 0.0};
        return r;
    }
    Complex r = {log(hypot(z.re, z.im)), atan2(z.im, z.re)};
    return r;
}

Complex complex_pow(Complex z, Complex w) {
    if (z.re == 0.0 && z.im == 0.0) {
        if (w.re > 0.0) { Complex zr = {0,0}; return zr; }
        Complex inf = {INFINITY, INFINITY};
        return inf;
    }
    Complex log_z = complex_log(z);
    Complex w_log_z = complex_mul(w, log_z);
    return complex_exp(w_log_z);
}

Complex complex_sqrt(Complex z) {
    double r = hypot(z.re, z.im);
    double theta = atan2(z.im, z.re);
    double sr = sqrt(r);
    Complex result = {sr * cos(theta/2.0), sr * sin(theta/2.0)};
    return result;
}

Complex complex_sin(Complex z) {
    Complex r = {sin(z.re)*cosh(z.im), cos(z.re)*sinh(z.im)};
    return r;
}

Complex complex_cos(Complex z) {
    Complex r = {cos(z.re)*cosh(z.im), -sin(z.re)*sinh(z.im)};
    return r;
}

Complex complex_tan(Complex z) {
    Complex s = complex_sin(z);
    Complex c = complex_cos(z);
    return complex_div(s, c);
}

Complex complex_sinh(Complex z) {
    Complex r = {sinh(z.re)*cos(z.im), cosh(z.re)*sin(z.im)};
    return r;
}

Complex complex_cosh(Complex z) {
    Complex r = {cosh(z.re)*cos(z.im), sinh(z.re)*sin(z.im)};
    return r;
}

Complex complex_tanh(Complex z) {
    Complex s = complex_sinh(z);
    Complex c = complex_cosh(z);
    return complex_div(s, c);
}

Complex complex_poly_eval(const ComplexPoly *p, Complex z) {
    if (!p || !p->coeffs || p->degree < 0) {
        Complex zero = {0.0, 0.0};
        return zero;
    }
    Complex result = p->coeffs[p->degree];
    for (int k = p->degree - 1; k >= 0; k--) {
        result = complex_mul(result, z);
        result = complex_add(result, p->coeffs[k]);
    }
    return result;
}

void roots_of_unity(int n, Complex *out) {
    if (n <= 0 || !out) return;
    for (int k = 0; k < n; k++) {
        double theta = 2.0 * M_PI * k / (double)n;
        out[k].re = cos(theta);
        out[k].im = sin(theta);
    }
}

int is_on_unit_circle(Complex z, double tol) {
    return (fabs(complex_abs(z) - 1.0) < tol);
}

static int cmplx_eq(Complex a, Complex b, double tol) {
    return (fabs(a.re - b.re) <= tol) && (fabs(a.im - b.im) <= tol);
}

int complex_field_axioms_self_test(const Complex *vals, int n, double tol) {
    if (!vals || n < 1) return 0;
    for (int i = 0; i < n; i++) {
        Complex a = vals[i];
        if (!cmplx_eq(complex_add(a, COMPLEX_ZERO), a, tol)) return 0;
        if (!cmplx_eq(complex_mul(a, COMPLEX_ONE), a, tol)) return 0;
        Complex na = complex_neg(a);
        if (!cmplx_eq(complex_add(a, na), COMPLEX_ZERO, tol)) return 0;
        if (!(a.re==0.0 && a.im==0.0)) {
            Complex one = {1,0};
            Complex inv = complex_div(one, a);
            if (!cmplx_eq(complex_mul(a, inv), one, tol)) return 0;
        }
        for (int j = 0; j < n; j++) {
            Complex b = vals[j];
            if (!cmplx_eq(complex_add(a,b), complex_add(b,a), tol)) return 0;
            if (!cmplx_eq(complex_mul(a,b), complex_mul(b,a), tol)) return 0;
            for (int k = 0; k < n; k++) {
                Complex c = vals[k];
                if (!cmplx_eq(complex_add(complex_add(a,b),c),
                              complex_add(a,complex_add(b,c)), tol)) return 0;
                if (!cmplx_eq(complex_mul(complex_mul(a,b),c),
                              complex_mul(a,complex_mul(b,c)), tol)) return 0;
                Complex bpc = complex_add(b,c);
                if (!cmplx_eq(complex_mul(a,bpc),
                              complex_add(complex_mul(a,b),complex_mul(a,c)), tol))
                    return 0;
            }
        }
    }
    return 1;
}


double complex_distance(Complex a, Complex b) {
    return hypot(a.re - b.re, a.im - b.im);
}

int is_real(Complex z, double tol) {
    return (fabs(z.im) <= tol);
}

int is_pure_imaginary(Complex z, double tol) {
    return (fabs(z.re) <= tol) && (fabs(z.im) > tol);
}

Complex complex_rotate(Complex z, double angle_rad) {
    double c = cos(angle_rad), s = sin(angle_rad);
    Complex r = {z.re * c - z.im * s, z.re * s + z.im * c};
    return r;
}

Complex complex_scale(Complex z, double factor) {
    Complex r = {z.re * factor, z.im * factor};
    return r;
}

void complex_to_matrix(Complex z, double M[2][2]) {
    M[0][0] =  z.re; M[0][1] = -z.im;
    M[1][0] =  z.im; M[1][1] =  z.re;
}

Complex matrix_to_complex(const double M[2][2]) {
    Complex z = {M[0][0], M[1][0]};
    return z;
}

Complex complex_from_string(const char *s) {
    double re = 0.0, im = 0.0;
    if (s) {
        char sign;
        if (sscanf(s, "%lf %c %lfi", &re, &sign, &im) >= 2) {
            if (sign == '-') im = -im;
        }
    }
    Complex z = {re, im};
    return z;
}

void complex_to_string(Complex z, char *buf, int bufsz) {
    if (buf && bufsz > 0) {
        if (z.im >= 0.0)
            snprintf(buf, bufsz, "%.6g + %.6gi", z.re, z.im);
        else
            snprintf(buf, bufsz, "%.6g - %.6gi", z.re, -z.im);
    }
}

Complex complex_pow_int(Complex z, int n) {
    if (n == 0) { Complex one = {1,0}; return one; }
    if (n < 0) {
        Complex one = {1,0};
        Complex inv = complex_div(one, z);
        return complex_pow_int(inv, -n);
    }
    double r = complex_abs(z);
    double theta = complex_arg(z);
    double rn = pow(r, (double)n);
    Complex result = {rn * cos(n*theta), rn * sin(n*theta)};
    return result;
}

void complex_nth_roots(Complex z, int n, Complex *out) {
    if (n <= 0 || !out) return;
    double r = pow(complex_abs(z), 1.0 / n);
    double theta0 = complex_arg(z) / n;
    for (int k = 0; k < n; k++) {
        double theta = theta0 + 2.0 * M_PI * k / n;
        out[k].re = r * cos(theta);
        out[k].im = r * sin(theta);
    }
}

Complex complex_dot_product(const Complex *a, const Complex *b, int n) {
    Complex sum = {0,0};
    if (!a || !b || n <= 0) return sum;
    for (int i = 0; i < n; i++) {
        Complex conj_a = complex_conj(a[i]);
        sum = complex_add(sum, complex_mul(conj_a, b[i]));
    }
    return sum;
}

Complex complex_sinc(Complex z) {
    if (z.re == 0.0 && z.im == 0.0) { Complex one = {1,0}; return one; }
    Complex pi_z = {M_PI * z.re, M_PI * z.im};
    Complex sin_pi_z = complex_sin(pi_z);
    return complex_div(sin_pi_z, pi_z);
}

Complex complex_stirling_gamma(Complex z, int terms) {
    Complex half = {0.5, 0.0};
    Complex z_phalf = complex_pow(z, half);
    Complex z_powz = complex_pow(z, z);
    Complex neg_z = complex_neg(z);
    Complex e_negz = complex_exp(neg_z);
    double sqrt_2pi = sqrt(2.0 * M_PI);
    Complex result = complex_mul(z_phalf, z_powz);
    result = complex_mul(result, e_negz);
    result = complex_scale(result, sqrt_2pi);
    (void)terms;
    return result;
}


/* ========================================================================
 * L4: Complex plane geometry and transformations
 * ======================================================================== */

/* Reflection across the real axis: conj(z) */
Complex complex_reflect_real(Complex z) {
    return complex_conj(z);
}

/* Reflection across the imaginary axis: -conj(z) */
Complex complex_reflect_imag(Complex z) {
    Complex c = complex_conj(z);
    Complex r = {-c.re, c.im};
    return r;
}

/* Inversion in unit circle: 1/conj(z) = z/|z|^2 */
Complex complex_invert_unit_circle(Complex z) {
    double r2 = complex_abs2(z);
    if (r2 < 1e-30) return (Complex){INFINITY, INFINITY};
    Complex result = {z.re / r2, z.im / r2};
    return result;
}

/* Reflection across line through origin at angle theta */
Complex complex_reflect_line(Complex z, double theta) {
    Complex rotated = complex_rotate(z, -theta);
    Complex reflected = complex_conj(rotated);
    return complex_rotate(reflected, theta);
}

/* ========================================================================
 * L4: Complex norms and inequalities
 *
 * Triangle inequality: |z1 + z2| <= |z1| + |z2|
 * Reverse triangle:    ||z1| - |z2|| <= |z1 - z2|
 * Cauchy-Schwarz for complex vectors
 * ======================================================================== */

double triangle_inequality_check(Complex z1, Complex z2) {
    double lhs = complex_abs(complex_add(z1, z2));
    double rhs = complex_abs(z1) + complex_abs(z2);
    return rhs - lhs;  /* >= 0 if inequality holds */
}

double reverse_triangle_check(Complex z1, Complex z2) {
    double lhs = fabs(complex_abs(z1) - complex_abs(z2));
    double rhs = complex_abs(complex_sub(z1, z2));
    return rhs - lhs;  /* >= 0 if inequality holds */
}

/* ========================================================================
 * L4: Complex sign function (generalized signum)
 *
 * csgn(z) = z/|z| if z != 0, 0 otherwise
 * Maps non-zero z to unit circle.
 * ======================================================================== */

Complex complex_sign(Complex z) {
    double mag = complex_abs(z);
    if (mag < 1e-30) return (Complex){0.0, 0.0};
    Complex result = {z.re / mag, z.im / mag};
    return result;
}

/* ========================================================================
 * L4: Principal argument branches
 *
 * Arg in (-pi, pi] is the principal branch.
 * Other branches: arg_k(z) = Arg(z) + 2*pi*k
 * ======================================================================== */

double complex_arg_branch(Complex z, int k) {
    return complex_arg(z) + 2.0 * M_PI * k;
}

/* ========================================================================
 * L4: Complex linear fractional interpolation
 * (useful for conformal mapping design between 3 pairs of points)
 * ======================================================================== */

Complex linear_fractional_interpolate(Complex z,
                                       Complex z1, Complex w1,
                                       Complex z2, Complex w2,
                                       Complex z3, Complex w3) {
    /* Cross-ratio preserving map: (w,w1,w2,w3) = (z,z1,z2,z3) */
    Complex num1 = complex_sub(z, z1);
    Complex den1 = complex_sub(z, z3);
    Complex num2 = complex_sub(z2, z3);
    Complex den2 = complex_sub(z2, z1);
    Complex cr_z = complex_div(complex_mul(num1, num2), complex_mul(den1, den2));
    Complex num1w = complex_sub(w2, w3);
    Complex den1w = complex_sub(w2, w1);
    Complex cr_w = complex_div(num1w, den1w);
    Complex cr_ratio = complex_div(cr_w, cr_z);
    Complex num = complex_mul(cr_ratio, complex_sub(z, z1));
    num = complex_mul(num, complex_sub(w2, w3));
    Complex den = complex_mul(complex_sub(z, z3), complex_sub(w2, w1));
    Complex term = complex_div(num, den);
    Complex result = complex_div(
        complex_sub(w1, complex_mul(w3, term)),
        complex_sub((Complex){1,0}, term)
    );
    return result;
}

/* ========================================================================
 * L4: Principal value of complex arctangent
 *
 * atan(z) = (i/2)*log((i+z)/(i-z))
 *
 * Branch cuts: (-i*inf, -i] and [+i, +i*inf) on imaginary axis
 * ======================================================================== */

Complex complex_atan(Complex z) {
    Complex i = {0.0, 1.0};
    Complex half_i = {0.0, 0.5};
    Complex num = complex_add(i, z);
    Complex den = complex_sub(i, z);
    Complex ratio = complex_div(num, den);
    Complex log_ratio = complex_log(ratio);
    Complex result = complex_mul(half_i, log_ratio);
    return result;
}

/* ========================================================================
 * L4: Complex arctangent with two arguments
 * atan2(y,x) returns arg(x + iy)
 * ======================================================================== */

Complex complex_atan2(Complex y, Complex x) {
    Complex z = complex_div(y, x);
    return complex_atan(z);
}

