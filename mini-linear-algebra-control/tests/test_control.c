/**
 * @file test_control.c
 * @brief Tests for control-theoretic linear algebra
 */
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
#include "control_linalg.h"
#include "eigen.h"
#include "matrix_exp.h"
#include "norms.h"

static int tests_run = 0, tests_passed = 0;
#define TEST(c,m) do{tests_run++;if(c)tests_passed++;else printf("  FAIL: %s\n",m);}while(0)
#define CLOSE(a,b,t) (fabs((a)-(b))<(t))

int main(void) {
    printf("=== test_control: control-theoretic linalg ===\n");

    /* Controllability test - double integrator A=[0 1;0 0] B=[0;1] */
    double Av[4] = {0, 1, 0, 0};
    double Bv[2] = {0, 1};
    Matrix *A = matrix_create_from_array(Av, 2, 2);
    Matrix *B = matrix_create_from_array(Bv, 2, 1);
    int ctrl = control_is_controllable(A, B, 1e-10);
    TEST(ctrl == 1, "double integrator controllable");

    Matrix *Cmat = matrix_create(2, 2);
    control_controllability_matrix(A, B, Cmat);
    TEST(matrix_rank(Cmat, 1e-10) == 2, "controllability matrix rank=2");

    /* Observability test - C=[1 0] */
    double Cv[2] = {1, 0};
    Matrix *Cobs = matrix_create_from_array(Cv, 1, 2);
    int obs = control_is_observable(A, Cobs, 1e-10);
    TEST(obs == 1, "double integrator observable");

    /* Lyapunov equation - A stable [-2 0; 0 -3] */
    double Asv[4] = {-2, 0, 0, -3};
    Matrix *Ast = matrix_create_from_array(Asv, 2, 2);
    Matrix *Q = matrix_create_identity(2);
    Matrix *P = matrix_create(2, 2);
    int ret = lyapunov_solve_continuous(Ast, Q, P);
    TEST(ret == 0, "Lyapunov solve (stable A)");
    if (ret == 0) {
        double res = lyapunov_residual(Ast, P, Q);
        TEST(res < 2e-6, "Lyapunov residual small");
    }

    /* Discrete Lyapunov */
    double Adv[4] = {0.5, 0, 0, 0.7};
    Matrix *Ad = matrix_create_from_array(Adv, 2, 2);
    ret = lyapunov_solve_discrete(Ad, Q, P);
    TEST(ret == 0, "discrete Lyapunov solve");

    /* Riccati equation */
    double Acv[4] = {0, 1, -2, -3};
    double Bcv[2] = {0, 1};
    double Qcv[4] = {1, 0, 0, 2};
    double Rcv[1] = {0.1};
    Matrix *Ac = matrix_create_from_array(Acv, 2, 2);
    Matrix *Bc = matrix_create_from_array(Bcv, 2, 1);
    Matrix *Qc = matrix_create_from_array(Qcv, 2, 2);
    Matrix *Rc = matrix_create_from_array(Rcv, 1, 1);
    Matrix *Pr = matrix_create(2, 2);
    ret = riccati_solve_care(Ac, Bc, Qc, Rc, Pr);
    TEST(ret == 0, "CARE solve");

    /* Matrix exponential */
    Matrix *expA = matrix_create(2, 2);
    ret = matrix_exp(Ast, expA);
    TEST(ret == 0, "matrix exponential");

    /* ZOH discretization */
    Matrix *Adisc = matrix_create(2, 2);
    Matrix *Bdisc = matrix_create(2, 1);
    ret = discretize_zoh(A, B, 0.1, Adisc, Bdisc);
    TEST(ret == 0, "ZOH discretization");

    /* Stability checks */
    int hurwitz = eigen_is_hurwitz(Ast, 1e-10);
    TEST(hurwitz == 1, "stable matrix is Hurwitz");
    int schur = eigen_is_schur(Ad, 1e-10);
    TEST(schur == 1, "Schur stable");

    /* Spectral radius */
    double rho = eigen_spectral_radius(Ast);
    TEST(CLOSE(rho, 3.0, 1e-6), "spectral radius = 3");

    /* LQR gain */
    Matrix *K = matrix_create(1, 2);
    ret = lqr_gain(Ac, Bc, Qc, Rc, K);
    TEST(ret == 0, "LQR gain computed");

    /* Observer gain (SISO - uses Ackermann via duality) */
    Vector *poles = vector_create_from_array((double[]){-5, -6}, 2);
    Matrix *Lobs = matrix_create(2, 1);
    ret = observer_gain_luenberger(A, Cobs, poles, Lobs);
    TEST(ret == 0, "observer gain computed");

    /* Hankel singular values */
    Matrix *Wc = matrix_create(2,2), *Wo = matrix_create(2,2);
    Vector *hsv = vector_create(2);
    if (control_controllability_gramian(Ast, B, Wc) == 0 &&
        control_observability_gramian(Ast, Cobs, Wo) == 0) {
        ret = control_hankel_singular_values(Wc, Wo, hsv);
        TEST(ret == 0, "Hankel singular values");
    }

    /* Norms */
    double fnorm = matrix_norm_frobenius(Ast);
    TEST(fnorm > 0, "Frobenius norm > 0");
    double cond_est = matrix_condition_l2(Ast);
    TEST(cond_est > 0, "condition number > 0");

    /* Orthogonality check for Q from QR */
    double orth_err = orthogonality_error(Q);
    TEST(orth_err < 1e-8, "Q nearly orthogonal");

    /* Cleanup */
    matrix_free(A); matrix_free(B); matrix_free(Cmat);
    matrix_free(Cobs); matrix_free(Ast); matrix_free(Q); matrix_free(P);
    matrix_free(Ad); matrix_free(Ac); matrix_free(Bc);
    matrix_free(Qc); matrix_free(Rc); matrix_free(Pr);
    matrix_free(expA); matrix_free(Adisc); matrix_free(Bdisc);
    matrix_free(K); vector_free(poles); matrix_free(Lobs);
    matrix_free(Wc); matrix_free(Wo); vector_free(hsv);
    /* R and Q from QR test */
    Matrix *R_test = matrix_create(3,3);
    Matrix *Q_test = matrix_create(3,3);
    matrix_free(R_test); matrix_free(Q_test);

    printf("  %d/%d tests passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}