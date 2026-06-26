/**
 * @file test_stability.c
 * @brief Tests for Routh-Hurwitz and Jury stability criteria
 */
#include "../include/stability.h"
#include "../include/laplace_core.h"
#include "../include/z_transform.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define TOL 1e-6

static int tests_run = 0, tests_passed = 0;
#define TEST(n) do { tests_run++; printf("  TEST: %s ... ", n); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(m) do { printf("FAIL: %s\n", m); } while(0)
#define CHECK(c, m) do { if (!(c)) { FAIL(m); return; } } while(0)

static void test_routh_stable(void)
{
    TEST("Routh-Hurwitz: s²+2s+1 is stable");
    LaplacePolynomial den = {2, {1.0, 2.0, 1.0}};
    RouthResult result;
    int ret = routh_stability(&den, &result);
    CHECK(ret == 0, "routh ok");
    CHECK(result.is_stable == 1, "stable");
    PASS();
}

static void test_routh_unstable(void)
{
    TEST("Routh-Hurwitz: s²-s+1 is unstable");
    LaplacePolynomial den = {2, {1.0, -1.0, 1.0}};
    RouthResult result;
    routh_stability(&den, &result);
    CHECK(result.is_stable == 0, "unstable detected");
    PASS();
}

static void test_routh_third_order(void)
{
    TEST("Routh-Hurwitz: s³+3s²+3s+1 = (s+1)³ is stable");
    LaplacePolynomial den = {3, {1.0, 3.0, 3.0, 1.0}};
    RouthResult result;
    routh_stability(&den, &result);
    CHECK(result.is_stable == 1, "3rd order stable");
    PASS();
}

static void test_jury_stable(void)
{
    TEST("Jury: 1 - 0.5z⁻¹ is stable (pole at z=0.5)");
    ZPoly poly = {1, {1.0, -0.5}};
    JuryResult result;
    int ret = jury_stability(&poly, &result);
    CHECK(ret == 0, "jury ok");
    CHECK(result.is_stable == 1, "stable (|0.5|<1)");
    PASS();
}

static void test_jury_unstable(void)
{
    TEST("Jury: 1 - 2z⁻¹ is unstable (pole at z=2)");
    ZPoly poly = {1, {1.0, -2.0}};
    JuryResult result;
    jury_stability(&poly, &result);
    CHECK(result.is_stable == 0, "unstable (|2|>1)");
    PASS();
}

int main(void)
{
    printf("=== Stability Tests ===\n\n");
    test_routh_stable();
    test_routh_unstable();
    test_routh_third_order();
    test_jury_stable();
    test_jury_unstable();
    printf("\n=== %d/%d passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
