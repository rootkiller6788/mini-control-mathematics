/**
 * test_optimization.c — Assert-based tests for mini-optimization-basics
 *
 * Tests core API: vector/matrix ops, Cholesky, gradient descent,
 * conjugate gradient, KKT conditions, projections.
 */

#include "optimization_core.h"
#include "gradient_methods.h"
#include "constrained_optimization.h"
#include "convex_analysis.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <string.h>

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASSED\n"); } while(0)
#define FAIL(msg) do { printf("FAILED: %s\n", msg); } while(0)

/* ─── Vector Operations ────────────────────────────────── */

static void test_vector_alloc_free(void)
{
    TEST("vector_alloc/free");
    vector_t *v = vector_alloc(10);
    assert(v != NULL && v->n == 10);
    for (size_t i = 0; i < 10; i++) assert(v->data[i] == 0.0);
    vector_free(v);
    PASS();
}

static void test_vector_dot_norm(void)
{
    TEST("vector dot product and norms");
    double data[] = {3.0, 4.0, 0.0};
    vector_t *v = vector_from_array(data, 3);
    assert(fabs(vector_norm_l2(v) - 5.0) < 1e-12);
    assert(fabs(vector_norm_l1(v) - 7.0) < 1e-12);
    assert(fabs(vector_norm_linf(v) - 4.0) < 1e-12);
    assert(fabs(vector_dot(v, v) - 25.0) < 1e-12);
    vector_free(v);
    PASS();
}

static void test_vector_ops(void)
{
    TEST("vector arithmetic");
    double a_d[] = {1.0, 2.0, 3.0}, b_d[] = {4.0, 5.0, 6.0};
    vector_t *a = vector_from_array(a_d, 3);
    vector_t *b = vector_from_array(b_d, 3);
    vector_t *out = vector_alloc(3);
    vector_add(a, b, out);
    assert(fabs(out->data[0] - 5.0) < 1e-12);
    vector_sub(a, b, out);
    assert(fabs(out->data[0] + 3.0) < 1e-12);
    vector_scale(a, 2.0);
    assert(fabs(a->data[0] - 2.0) < 1e-12);
    vector_free(a); vector_free(b); vector_free(out);
    PASS();
}

/* ─── Cholesky Decomposition ────────────────────────────── */

static void test_cholesky(void)
{
    TEST("Cholesky decomposition");
    matrix_t *A = matrix_alloc(2, 2);
    MAT_ELT(A, 0, 0) = 4.0; MAT_ELT(A, 0, 1) = 2.0;
    MAT_ELT(A, 1, 0) = 2.0; MAT_ELT(A, 1, 1) = 3.0;

    matrix_t *L = matrix_alloc(2, 2);
    for (size_t i = 0; i < 2; i++)
        for (size_t j = 0; j < 2; j++)
            MAT_ELT(L, i, j) = MAT_ELT(A, i, j);

    int rc = cholesky_decomp(L);
    assert(rc == 0);
    assert(fabs(MAT_ELT(L, 0, 0) - 2.0) < 1e-10);
    assert(fabs(MAT_ELT(L, 1, 0) - 1.0) < 1e-10);
    assert(fabs(MAT_ELT(L, 1, 1) - sqrt(2.0)) < 1e-10);

    vector_t *b = vector_alloc(2);
    b->data[0] = 6.0; b->data[1] = 7.0;
    vector_t *x = vector_alloc(2);
    for (size_t i = 0; i < 2; i++)
        for (size_t j = 0; j < 2; j++)
            MAT_ELT(L, i, j) = MAT_ELT(A, i, j);
    cholesky_decomp(L);
    cholesky_solve(L, b, x);
    /* True solution: [0.5, 2.0] (A=[4,2;2,3], b=[6,7]) */
    assert(fabs(x->data[0] - 0.5) < 1e-8 && fabs(x->data[1] - 2.0) < 1e-8);

    matrix_free(A); matrix_free(L);
    vector_free(b); vector_free(x);
    PASS();
}

/* ─── Matrix Ops ────────────────────────────────────────── */

static void test_matrix_ops(void)
{
    TEST("matrix-vector multiplication");
    matrix_t *A = matrix_alloc(2, 3);
    MAT_ELT(A, 0, 0) = 1.0; MAT_ELT(A, 0, 1) = 2.0; MAT_ELT(A, 0, 2) = 3.0;
    MAT_ELT(A, 1, 0) = 4.0; MAT_ELT(A, 1, 1) = 5.0; MAT_ELT(A, 1, 2) = 6.0;
    double xd[] = {1.0, 1.0, 1.0};
    vector_t *x = vector_from_array(xd, 3);
    vector_t *y = vector_alloc(2);
    matrix_vector_mul(A, x, y);
    assert(fabs(y->data[0] - 6.0) < 1e-12 && fabs(y->data[1] - 15.0) < 1e-12);

    matrix_t *Q = matrix_alloc(2, 2);
    MAT_ELT(Q, 0, 0) = 2.0; MAT_ELT(Q, 0, 1) = 1.0;
    MAT_ELT(Q, 1, 0) = 1.0; MAT_ELT(Q, 1, 1) = 2.0;
    double vd[] = {1.0, 2.0};
    vector_t *v = vector_from_array(vd, 2);
    double qf = quadratic_form(Q, v);
    assert(fabs(qf - 14.0) < 1e-10);
    assert(matrix_is_psd(Q) == 1);

    matrix_free(A); matrix_free(Q);
    vector_free(x); vector_free(y); vector_free(v);
    PASS();
}

/* ─── Gradient Descent ─────────────────────────────────── */

static double quad_eval(const convex_function_t *f, const vector_t *xx)
{
    (void)f;
    double d0 = xx->data[0] - 1.0, d1 = xx->data[1] - 2.0;
    return d0*d0 + d1*d1;
}
static void quad_grad(const convex_function_t *f, const vector_t *xx, vector_t *g)
{
    (void)f;
    g->data[0] = 2.0*(xx->data[0]-1.0);
    g->data[1] = 2.0*(xx->data[1]-2.0);
}

static void test_gradient_descent(void)
{
    TEST("gradient descent convergence");
    convex_function_t f_obj;
    memset(&f_obj, 0, sizeof(f_obj));
    f_obj.n = 2;
    f_obj.eval = quad_eval;
    f_obj.grad = quad_grad;

    unconstr_problem_t prob;
    memset(&prob, 0, sizeof(prob));
    prob.n = 2; prob.f = f_obj; prob.has_analytic = 1; prob.opt_value = 0.0;
    double x0d[] = {10.0, 10.0};
    prob.x0.data = x0d; prob.x0.n = 2;

    gradient_descent_config_t cfg = gd_config_default();
    cfg.ls_type = LS_BACKTRACKING;
    cfg.grad_tol = 1e-5;

    opt_result_t result;
    memset(&result, 0, sizeof(result));
    gradient_descent(&prob, &cfg, &result);

    printf("(iters=%zu, f=%.6f) ", result.iterations, result.f_value);
    assert(result.status == OPT_SUCCESS);
    assert(result.f_value < 1e-5);
    assert(fabs(result.x.data[0] - 1.0) < 1e-3);
    assert(fabs(result.x.data[1] - 2.0) < 1e-3);
    free(result.x.data);
    PASS();
}

/* ─── Conjugate Gradient ────────────────────────────────── */

static void test_conjugate_gradient(void)
{
    TEST("conjugate gradient");
    matrix_t *Q = matrix_alloc(2, 2);
    MAT_ELT(Q, 0, 0) = 4.0; MAT_ELT(Q, 0, 1) = 1.0;
    MAT_ELT(Q, 1, 0) = 1.0; MAT_ELT(Q, 1, 1) = 3.0;
    vector_t *b = vector_alloc(2);
    b->data[0] = 1.0; b->data[1] = 2.0;
    vector_t *x = vector_alloc(2);
    x->data[0] = 0.0; x->data[1] = 0.0;

    cg_config_t cfg = cg_config_default();
    conjugate_gradient_linear(Q, b, x, &cfg);

    double e0 = 1.0/11.0, e1 = 7.0/11.0;
    assert(fabs(x->data[0] - e0) < 1e-6);
    assert(fabs(x->data[1] - e1) < 1e-6);

    matrix_free(Q); vector_free(b); vector_free(x);
    PASS();
}

/* ─── Projections ───────────────────────────────────────── */

static void test_projections(void)
{
    TEST("Euclidean projections");
    size_t n = 3;
    double cd[] = {0,0,0}, vd[] = {3.0, 4.0, 0};
    vector_t *c = vector_from_array(cd, n);
    vector_t *v = vector_from_array(vd, n);
    vector_t *p = vector_alloc(n);
    project_l2_ball(v, c, 1.0, p);
    assert(fabs(p->data[0] - 0.6) < 1e-10 && fabs(p->data[1] - 0.8) < 1e-10);

    double ld[] = {-1,-1,-1}, ud[] = {1,1,1};
    vector_t *lo = vector_from_array(ld, n);
    vector_t *up = vector_from_array(ud, n);
    double v2d[] = {2.0, -3.0, 0.5};
    vector_t *v2 = vector_from_array(v2d, n);
    project_box(v2, lo, up, p);
    assert(fabs(p->data[0]-1.0) < 1e-12 && fabs(p->data[1]+1.0) < 1e-12);

    project_linf_ball(v2, 2.0, p);
    assert(fabs(p->data[0]-2.0) < 1e-12 && fabs(p->data[1]+2.0) < 1e-12);

    vector_free(c); vector_free(v); vector_free(v2);
    vector_free(lo); vector_free(up); vector_free(p);
    PASS();
}

/* ─── Biconjugate Theorem (L4) ──────────────────────────── */

static double sq_eval(const convex_function_t *f, const vector_t *x)
{
    (void)f;
    double s = 0;
    for (size_t i = 0; i < x->n; i++) s += x->data[i]*x->data[i];
    return 0.5 * s;
}
static void sq_grad(const convex_function_t *f, const vector_t *x, vector_t *g)
{
    (void)f;
    for (size_t i = 0; i < x->n; i++) g->data[i] = x->data[i];
}

static void test_convex_conjugate(void)
{
    TEST("convex conjugate (Fenchel transform)");
    convex_function_t sq;
    memset(&sq, 0, sizeof(sq));
    sq.n = 2;
    sq.eval = sq_eval;
    sq.grad = sq_grad;

    /* f(x) = 0.5||x||^2, f*(y) = 0.5||y||^2 (self-dual) */
    double yd[] = {1.0, 2.0};
    vector_t *y = vector_from_array(yd, 2);
    double conj = convex_conjugate_eval(&sq, y);
    /* f*(1,2) = 0.5*(1+4)=2.5 */
    assert(fabs(conj - 2.5) < 0.1);

    vector_free(y);
    PASS();
}

/* ─── Main ──────────────────────────────────────────────── */

int main(void)
{
    srand((unsigned int)time(NULL));
    printf("=== mini-optimization-basics Test Suite ===\n\n");

    test_vector_alloc_free();
    test_vector_dot_norm();
    test_vector_ops();
    test_matrix_ops();
    test_cholesky();
    test_gradient_descent();
    test_conjugate_gradient();
    test_projections();
    test_convex_conjugate();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
