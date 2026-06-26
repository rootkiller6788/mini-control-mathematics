/**
 * @file complex_number.h
 * @brief Complex number type definitions and elementary operations
 *
 * L1 --- Definitions: Complex number as ordered pair (a + ib),
 * polar form (r--e^{i--}), algebraic structure (field over C)
 *
 * Course alignment:
 *   MIT 18.04 Complex Variables with Applications
 *   Stanford MATH 106 Functions of a Complex Variable
 *   Berkeley MATH 185 Complex Analysis
 *   Cambridge/Michigan/Purdue/ETH/Tsinghua: core complex arithmetic
 *
 * Textbook: Churchill & Brown, "Complex Variables and Applications" (2014)
 *           Ahlfors, "Complex Analysis" (1979)
 */

#ifndef COMPLEX_NUMBER_H
#define COMPLEX_NUMBER_H

#ifdef _WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stddef.h>

/* --------------------------------------------------------------------------
 * L1: Core Type --- complex number as a + ib
 * -------------------------------------------------------------------------- */

/** Complex number in Cartesian form: z = re + i--im */
typedef struct {
    double re;  /**< real part */
    double im;  /**< imaginary part */
} Complex;

/* --------------------------------------------------------------------------
 * L1: Constants
 * -------------------------------------------------------------------------- */

/** The complex zero: 0 + 0i */
static const Complex COMPLEX_ZERO = {0.0, 0.0};

/** The complex one: 1 + 0i */
static const Complex COMPLEX_ONE = {1.0, 0.0};

/** The imaginary unit: 0 + 1i */
static const Complex COMPLEX_I = {0.0, 1.0};

/* --------------------------------------------------------------------------
 * L2: Elementary complex arithmetic (algebraic field operations)
 * -------------------------------------------------------------------------- */

Complex complex_add(Complex z1, Complex z2);
Complex complex_sub(Complex z1, Complex z2);
Complex complex_mul(Complex z1, Complex z2);
Complex complex_div(Complex z1, Complex z2);
Complex complex_neg(Complex z);
Complex complex_conj(Complex z);

/* --------------------------------------------------------------------------
 * L2: magnitude, argument, polar form
 * -------------------------------------------------------------------------- */

double complex_abs(Complex z);
double complex_abs2(Complex z);
double complex_arg(Complex z);
Complex complex_polar(double r, double theta);
double complex_radius(Complex z);
double complex_angle(Complex z);

/* --------------------------------------------------------------------------
 * L2: exp, log, power, sqrt --- transcendental functions
 * -------------------------------------------------------------------------- */

Complex complex_exp(Complex z);
Complex complex_log(Complex z);
Complex complex_pow(Complex z, Complex w);
Complex complex_sqrt(Complex z);

/* --------------------------------------------------------------------------
 * L2: trigonometric and hyperbolic functions
 * -------------------------------------------------------------------------- */

Complex complex_sin(Complex z);
Complex complex_cos(Complex z);
Complex complex_tan(Complex z);
Complex complex_sinh(Complex z);
Complex complex_cosh(Complex z);
Complex complex_tanh(Complex z);

/* --------------------------------------------------------------------------
 * L3: Complex polynomial representation
 * -------------------------------------------------------------------------- */

typedef struct {
    int      degree;
    Complex *coeffs;
} ComplexPoly;

Complex complex_poly_eval(const ComplexPoly *p, Complex z);
void roots_of_unity(int n, Complex *out);
int is_on_unit_circle(Complex z, double tol);

/* --------------------------------------------------------------------------
 * L4: Field axioms verification utilities
 * -------------------------------------------------------------------------- */


/* --------------------------------------------------------------------------
 * Additional geometric and utility operations
 * -------------------------------------------------------------------------- */

double complex_distance(Complex a, Complex b);
int is_real(Complex z, double tol);
int is_pure_imaginary(Complex z, double tol);
Complex complex_rotate(Complex z, double angle_rad);
Complex complex_scale(Complex z, double factor);
void complex_to_matrix(Complex z, double M[2][2]);
Complex matrix_to_complex(const double M[2][2]);
Complex complex_from_string(const char *s);
void complex_to_string(Complex z, char *buf, int bufsz);
Complex complex_pow_int(Complex z, int n);
void complex_nth_roots(Complex z, int n, Complex *out);
Complex complex_dot_product(const Complex *a, const Complex *b, int n);
Complex complex_sinc(Complex z);
Complex complex_stirling_gamma(Complex z, int terms);

/* Additional complex functions */
Complex complex_reflect_real(Complex z);
Complex complex_reflect_imag(Complex z);
Complex complex_invert_unit_circle(Complex z);
Complex complex_reflect_line(Complex z, double theta);
double triangle_inequality_check(Complex z1, Complex z2);
double reverse_triangle_check(Complex z1, Complex z2);
Complex complex_sign(Complex z);
double complex_arg_branch(Complex z, int k);
Complex complex_atan(Complex z);
Complex complex_atan2(Complex y, Complex x);
Complex linear_fractional_interpolate(Complex z, Complex z1, Complex w1,
                                       Complex z2, Complex w2,
                                       Complex z3, Complex w3);

int complex_field_axioms_self_test(const Complex *test_values, int n, double tol);

#endif /* COMPLEX_NUMBER_H */
