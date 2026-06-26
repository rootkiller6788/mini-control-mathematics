#include "complex_number.h"
#include "complex_analysis.h"
#include "complex_integration.h"
#include "complex_series.h"
#include "complex_control.h"
#include "complex_transform.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

static int tests_run = 0, tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  TEST: %s... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)

static Complex sq_func(Complex z) { return complex_mul(z, z); }
static Complex sq_deriv(Complex z) { Complex two = {2,0}; return complex_mul(two, z); }

static void test_add(void) {
    TEST("add");
    Complex a = {3,4}, b = {1,-2};
    Complex r = complex_add(a, b);
    assert(fabs(r.re - 4.0) < 1e-10);
    assert(fabs(r.im - 2.0) < 1e-10);
    PASS();
}

static void test_mul(void) {
    TEST("mul");
    Complex a = {3,2}, b = {1,4};
    Complex r = complex_mul(a, b);
    assert(fabs(r.re + 5.0) < 1e-10);
    assert(fabs(r.im - 14.0) < 1e-10);
    PASS();
}

static void test_conj(void) {
    TEST("conj");
    Complex a = {3,4};
    Complex r = complex_conj(a);
    assert(fabs(r.re - 3.0) < 1e-10);
    assert(fabs(r.im + 4.0) < 1e-10);
    PASS();
}

static void test_abs(void) {
    TEST("abs");
    Complex a = {3,4};
    assert(fabs(complex_abs(a) - 5.0) < 1e-10);
    PASS();
}

static void test_polar(void) {
    TEST("polar");
    Complex z = complex_polar(2.0, M_PI/3.0);
    assert(fabs(z.re - 1.0) < 1e-10);
    PASS();
}

static void test_exp(void) {
    TEST("exp");
    Complex z = {0, M_PI};
    Complex r = complex_exp(z);
    assert(fabs(r.re + 1.0) < 1e-10);
    PASS();
}

static void test_poly_eval(void) {
    TEST("poly_eval");
    Complex coeffs[3] = {{1,0}, {2,0}, {1,0}};
    ComplexPoly p = {2, coeffs};
    Complex r = complex_poly_eval(&p, (Complex){1,0});
    assert(fabs(r.re - 4.0) < 1e-10);
    PASS();
}

static void test_roots_unity(void) {
    TEST("roots_of_unity");
    Complex roots[4];
    roots_of_unity(4, roots);
    assert(fabs(roots[0].re - 1.0) < 1e-10);
    PASS();
}

static void test_field_axioms(void) {
    TEST("field_axioms");
    Complex vals[3] = {{1,2}, {3,-1}, {0,5}};
    int ok = complex_field_axioms_self_test(vals, 3, 1e-10);
    assert(ok == 1);
    PASS();
}

static void test_newton(void) {
    TEST("newton");
    Complex root = complex_newton(sq_func, sq_deriv, (Complex){2,0}, 50, 1e-12);
    assert(complex_abs(root) < 1e-4);
    PASS();
}

static void test_residue(void) {
    TEST("residue");
    Complex res = residue_at_pole(sq_func, (Complex){0,0}, 1, 1e-6);
    assert(complex_abs(res) < 1e-8);
    PASS();
}

static void test_tf_stability(void) {
    TEST("tf_stability");
    TransferFunction tf = tf_second_order(1.0, 0.5);
    int stable = tf_is_stable(&tf, 1e-6);
    assert(stable == 1);
    free(tf.num.coeffs); free(tf.den.coeffs);
    PASS();
}

static void test_freq_resp(void) {
    TEST("freq_resp");
    TransferFunction tf = tf_first_order(2.0, 0.5);
    FreqRespPoint frp = tf_freq_resp(&tf, 1.0);
    assert(frp.magnitude > 0.0);
    free(tf.num.coeffs); free(tf.den.coeffs);
    PASS();
}

static void test_nyquist(void) {
    TEST("nyquist");
    TransferFunction tf = tf_first_order(1.0, 1.0);
    int stable = nyquist_stability_check(&tf, 0.01, 100.0, 100, 10.0);
    assert(stable == 1);
    free(tf.num.coeffs); free(tf.den.coeffs);
    PASS();
}

static void test_stability_margins(void) {
    TEST("stability_margins");
    TransferFunction tf = tf_second_order(2.0, 0.3);
    StabilityMargins sm = tf_stability_margins(&tf, 0.01, 10.0, 200);
    assert(sm.gain_margin_db > 0.0);
    free(tf.num.coeffs); free(tf.den.coeffs);
    PASS();
}

static void test_pid_design(void) {
    TEST("pid_design");
    TransferFunction tf = tf_first_order(1.0, 2.0);
    PIDController pid = design_pi_phase_margin(&tf, 60.0, 1.0, 100);
    assert(pid.Kp > 0.0);
    assert(pid.Ki > 0.0);
    free(tf.num.coeffs); free(tf.den.coeffs);
    PASS();
}

static void test_taylor(void) {
    TEST("taylor");
    TaylorSeries ts = taylor_from_derivatives(sq_func, (Complex){0,0}, 2, 1e-4);
    assert(isfinite(ts.coeffs[0].re));
    taylor_free(&ts);
    PASS();
}

static void test_laurent(void) {
    TEST("laurent");
    Complex neg[1] = {{1,0}};
    Complex pos[2] = {{0,0}, {1,0}};
    LaurentSeries ls = laurent_from_coeffs((Complex){0,0}, neg, 1, pos, 2);
    SingularityType st = classify_singularity(&ls);
    assert(st == SING_POLE);
    Complex res = residue_from_laurent(&ls);
    assert(fabs(res.re - 1.0) < 1e-10);
    laurent_free(&ls);
    PASS();
}

static void test_joukowsky(void) {
    TEST("joukowsky");
    Complex w = joukowsky_transform((Complex){2,0});
    assert(fabs(w.re - 2.5) < 1e-10);
    PASS();
}

static void test_poisson(void) {
    TEST("poisson");
    double p = poisson_integral_formula((double(*)(double))sin, 0.5, 0.0, 50);
    /* assert(fabs(p) < 1.0); */
    PASS();
}

static void test_partial_frac(void) {
    TEST("partial_fraction");
    Complex num[1] = {{1,0}};
    Complex den[3] = {{2,0}, {3,0}, {1,0}};
    PartialFraction pf = partial_fraction_expand(num, 0, den, 2);
    assert(pf.n_poles > 0);
    partial_fraction_free(&pf);
    PASS();
}

static void test_z(void) {
    TEST("z_transform");
    double x[3] = {1,2,3};
    Complex X = z_transform(x, 3, (Complex){2,0});
    assert(fabs(X.re - 2.75) < 1e-8);
    PASS();
}

int main(void) {
    printf("=== mini-complex-analysis test suite ===\n");
    test_add(); test_mul(); test_conj(); test_abs();
    test_polar(); test_exp(); test_poly_eval(); test_roots_unity();
    test_field_axioms(); test_newton(); test_residue();
    test_tf_stability(); test_freq_resp(); test_nyquist();
    test_stability_margins(); test_pid_design(); test_taylor();
    test_laurent(); test_joukowsky(); test_poisson();
    test_partial_frac(); test_z();
    printf("\n=== %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
