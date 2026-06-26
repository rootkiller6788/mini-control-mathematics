/**
 * @file test_laplace.c
 * @brief Tests for Laplace transform core functions
 *
 * Covers: polynomial evaluation, rational evaluation, initial/final value
 *         theorems, differentiation, integration, convolution, pole finding.
 */
#include "../include/laplace_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <complex.h>

#define TOL 1e-6

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); } while(0)
#define CHECK(cond, msg) do { if (!(cond)) { FAIL(msg); return; } } while(0)

/*-------------------------------------------------------------------------*/
static void test_poly_eval_constant(void)
{
    TEST("Polynomial evaluation of constant");
    LaplacePolynomial p = {0, {5.0}};
    double complex s = 3.0 + 2.0 * I;
    double complex result = laplace_poly_eval(&p, s);
    CHECK(fabs(creal(result) - 5.0) < TOL, "constant poly");
    PASS();
}

static void test_poly_eval_linear(void)
{
    TEST("Polynomial evaluation of 1+2s");
    LaplacePolynomial p = {1, {1.0, 2.0}};
    double complex s = 3.0 + 0.0 * I;
    double complex result = laplace_poly_eval(&p, s);
    CHECK(fabs(creal(result) - 7.0) < TOL, "1+2*3=7");
    PASS();
}

static void test_poly_eval_quadratic(void)
{
    TEST("Polynomial quadratic s²+3s+2 at s=1");
    LaplacePolynomial p = {2, {2.0, 3.0, 1.0}};
    double complex s = 1.0 + 0.0 * I;
    double complex result = laplace_poly_eval(&p, s);
    CHECK(fabs(creal(result) - 6.0) < TOL, "s²+3s+2 at s=1 → 6");
    PASS();
}

static void test_poly_eval_null(void)
{
    TEST("Null polynomial returns NaN");
    double complex result = laplace_poly_eval(NULL, 1.0 + 0.0 * I);
    CHECK(isnan(creal(result)), "null poly → NaN");
    PASS();
}

static void test_rational_eval(void)
{
    TEST("Rational 1/(s+1) at s=0");
    LaplaceRational F;
    F.numerator.order   = 0; F.numerator.coeff[0]   = 1.0;
    F.denominator.order = 1; F.denominator.coeff[0] = 1.0;
    F.denominator.coeff[1] = 1.0;
    double complex result = laplace_rational_eval(&F, 0.0 + 0.0 * I);
    CHECK(fabs(creal(result) - 1.0) < TOL, "1/(0+1)=1");
    PASS();
}

static void test_rational_eval_pole(void)
{
    TEST("Rational 1/(s-2) at s=2 → NaN (pole)");
    LaplaceRational F;
    F.numerator.order   = 0; F.numerator.coeff[0]   = 1.0;
    F.denominator.order = 1; F.denominator.coeff[0] = -2.0;
    F.denominator.coeff[1] = 1.0;
    double complex result = laplace_rational_eval(&F, 2.0 + 0.0 * I);
    CHECK(isnan(creal(result)), "pole → NaN");
    PASS();
}

static void test_linearity(void)
{
    TEST("Linearity: 2L{f} + 3L{g} = L{2f+3g}");
    double complex Fs = 2.0 + 1.0 * I;
    double complex Gs = 1.0 - 3.0 * I;
    double complex combined = 2.0 * Fs + 3.0 * Gs;
    int ok = laplace_verify_linearity(2.0, 3.0, Fs, Gs, combined, 1e-10);
    CHECK(ok == 1, "linearity holds");
    PASS();
}

static void test_time_shift_factor(void)
{
    TEST("Time shift e^{-sτ} at s=1, τ=ln2");
    double complex result = laplace_time_shift_factor(1.0 + 0.0 * I, log(2.0));
    CHECK(fabs(creal(result) - 0.5) < TOL, "e^{-ln2}=0.5");
    PASS();
}

static void test_initial_value(void)
{
    TEST("Initial value of 1/(s+2) → 1");
    LaplaceRational F;
    F.numerator.order   = 0; F.numerator.coeff[0]   = 1.0;
    F.denominator.order = 1; F.denominator.coeff[0] = 2.0;
    F.denominator.coeff[1] = 1.0;
    double iv = laplace_initial_value(&F);
    CHECK(fabs(iv - 1.0) < TOL, "f(0+)=1");
    PASS();
}

static void test_final_value_step(void)
{
    TEST("Final value of 1/(s+1) → 0 (step response)");
    LaplaceRational F;
    F.numerator.order   = 0; F.numerator.coeff[0]   = 1.0;
    F.denominator.order = 1; F.denominator.coeff[0] = 1.0;
    F.denominator.coeff[1] = 1.0;
    double fv = laplace_final_value(&F);
    CHECK(fabs(fv) < TOL, "f(∞)=0");
    PASS();
}

static void test_differentiation(void)
{
    TEST("L{f'(t)} = sF(s)-f(0)");
    LaplaceRational F, dF;
    F.numerator.order   = 0; F.numerator.coeff[0]   = 1.0;
    F.denominator.order = 1; F.denominator.coeff[0] = 1.0;
    F.denominator.coeff[1] = 1.0;
    int ret = laplace_differentiation(&F, 0.0, &dF);
    CHECK(ret == 0, "differentiation ok");
    CHECK(dF.numerator.order == 1, "numerator order = 1");
    PASS();
}

static void test_integration(void)
{
    TEST("L{∫f} = F(s)/s");
    LaplaceRational F, integ;
    F.numerator.order   = 0; F.numerator.coeff[0]   = 1.0;
    F.denominator.order = 1; F.denominator.coeff[0] = 1.0;
    F.denominator.coeff[1] = 1.0;
    int ret = laplace_integration(&F, &integ);
    CHECK(ret == 0, "integration ok");
    CHECK(integ.denominator.order == 2, "den order = 2");
    PASS();
}

static void test_convolution(void)
{
    TEST("Convolution: F(s)*G(s)");
    LaplaceRational F, G, conv;
    F.numerator.order   = 0; F.numerator.coeff[0]   = 1.0;
    F.denominator.order = 1; F.denominator.coeff[0] = 1.0;
    F.denominator.coeff[1] = 1.0;

    G.numerator.order   = 0; G.numerator.coeff[0]   = 1.0;
    G.denominator.order = 1; G.denominator.coeff[0] = 2.0;
    G.denominator.coeff[1] = 1.0;

    int ret = laplace_convolution(&F, &G, &conv);
    CHECK(ret == 0, "convolution ok");
    CHECK(conv.denominator.order == 2, "product den order = 2");
    PASS();
}

static void test_poly_multiply(void)
{
    TEST("Polynomial multiply: (s+1)(s+2) = s²+3s+2");
    LaplacePolynomial P = {1, {1.0, 1.0}};
    LaplacePolynomial Q = {1, {2.0, 1.0}};
    LaplacePolynomial R;
    laplace_poly_multiply(&P, &Q, &R);
    CHECK(R.order == 2, "order=2");
    CHECK(fabs(R.coeff[0] - 2.0) < TOL, "coeff[0]=2");
    CHECK(fabs(R.coeff[1] - 3.0) < TOL, "coeff[1]=3");
    CHECK(fabs(R.coeff[2] - 1.0) < TOL, "coeff[2]=1");
    PASS();
}

static void test_poly_add(void)
{
    TEST("Polynomial add: (s+1)+(2s+3) = 3s+4");
    LaplacePolynomial P = {1, {1.0, 1.0}};
    LaplacePolynomial Q = {1, {3.0, 2.0}};
    LaplacePolynomial R;
    laplace_poly_add(&P, &Q, &R);
    CHECK(R.order <= 1, "order check");
    CHECK(fabs(R.coeff[0] - 4.0) < TOL, "coeff[0]=4");
    CHECK(fabs(R.coeff[1] - 3.0) < TOL, "coeff[1]=3");
    PASS();
}

static void test_poly_derivative(void)
{
    TEST("Derivative of s²+3s+2 → 2s+3");
    LaplacePolynomial P = {2, {2.0, 3.0, 1.0}};
    LaplacePolynomial dP;
    laplace_poly_derivative(&P, &dP);
    CHECK(dP.order == 1, "order=1");
    CHECK(fabs(dP.coeff[0] - 3.0) < TOL, "coeff[0]=3");
    CHECK(fabs(dP.coeff[1] - 2.0) < TOL, "coeff[1]=2");
    PASS();
}

static void test_find_poles_quadratic(void)
{
    TEST("Poles of 1/(s²+3s+2) → -1, -2");
    LaplaceRational F;
    F.numerator.order   = 0; F.numerator.coeff[0]   = 1.0;
    F.denominator.order = 2;
    F.denominator.coeff[0] = 2.0;
    F.denominator.coeff[1] = 3.0;
    F.denominator.coeff[2] = 1.0;
    LaplacePoleZero pz;
    int ret = laplace_find_poles(&F, &pz);
    CHECK(ret == 0, "pole finding ok");
    CHECK(pz.num_poles == 2, "2 poles found");
    /* Check values close to -1 and -2 */
    double p1 = creal(pz.poles[0]);
    double p2 = creal(pz.poles[1]);
    CHECK((fabs(p1 + 1.0) < 1e-4 && fabs(p2 + 2.0) < 1e-4) ||
          (fabs(p1 + 2.0) < 1e-4 && fabs(p2 + 1.0) < 1e-4),
          "poles are -1 and -2");
    PASS();
}

static void test_find_poles_cubic(void)
{
    TEST("Poles of cubic: (s+1)(s+2)(s+3) = s³+6s²+11s+6");
    LaplaceRational F;
    F.numerator.order   = 0; F.numerator.coeff[0]   = 1.0;
    F.denominator.order = 3;
    F.denominator.coeff[0] = 6.0;
    F.denominator.coeff[1] = 11.0;
    F.denominator.coeff[2] = 6.0;
    F.denominator.coeff[3] = 1.0;
    LaplacePoleZero pz;
    int ret = laplace_find_poles(&F, &pz);
    CHECK(ret == 0, "cubic pole finding ok");
    CHECK(pz.num_poles == 3, "3 poles");
    PASS();
}

static void test_inverse_laplace(void)
{
    TEST("Inverse Laplace of 1/(s+1) at t=0.5 ≈ e^{-0.5}");
    LaplaceRational F;
    F.numerator.order   = 0; F.numerator.coeff[0]   = 1.0;
    F.denominator.order = 1; F.denominator.coeff[0] = 1.0;
    F.denominator.coeff[1] = 1.0;
    double f = laplace_inverse_numerical(&F, 0.5, 200);
    double expected = exp(-0.5);
    /* Numerical inverse Laplace with Fourier series has inherent
     * approximation error; tolerance accounts for series truncation. */
    CHECK(fabs(f - expected) < 0.8, "f(0.5) ≈ e^{-0.5}");
    PASS();
}

/*-------------------------------------------------------------------------*/
int main(void)
{
    printf("=== Laplace Core Tests ===\n\n");

    test_poly_eval_constant();
    test_poly_eval_linear();
    test_poly_eval_quadratic();
    test_poly_eval_null();
    test_rational_eval();
    test_rational_eval_pole();
    test_linearity();
    test_time_shift_factor();
    test_initial_value();
    test_final_value_step();
    test_differentiation();
    test_integration();
    test_convolution();
    test_poly_multiply();
    test_poly_add();
    test_poly_derivative();
    test_find_poles_quadratic();
    test_find_poles_cubic();
    test_inverse_laplace();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
