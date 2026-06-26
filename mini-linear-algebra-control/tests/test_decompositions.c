/**
 * @file test_decompositions.c
 * @brief Tests for matrix decompositions
 */
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
#include "decompositions.h"
#include "solvers.h"
#include "norms.h"

static int tests_run = 0, tests_passed = 0;
#define TEST(c,m) do{tests_run++;if(c)tests_passed++;else printf("  FAIL: %s\n",m);}while(0)
#define CLOSE(a,b,t) (fabs((a)-(b))<(t))

int main(void) {
    printf("=== test_decompositions: LU, QR, Cholesky, SVD ===\n");

    /* LU decomposition */
    double Av[9] = {2, 1, 1, 4, -6, 0, -2, 7, 2};
    Matrix *A = matrix_create_from_array(Av, 3, 3);
    Matrix *LU = matrix_create_copy(A);
    int piv[3];
    int ret = lu_decompose(LU, piv);
    TEST(ret == 0, "LU decompose");

    Matrix *L = matrix_create(3,3), *U = matrix_create(3,3);
    lu_extract_l(LU, L); lu_extract_u(LU, U);
    Matrix *LUprod = matrix_create(3,3);
    matrix_multiply(L, U, LUprod);
    /* Reconstruct PA and verify PA = LU */
    Matrix *PA = matrix_create(3,3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            matrix_set(PA, i, j, matrix_get(A, piv[i], j));
    TEST(matrix_equal(PA, LUprod, 1e-8), "PA = LU");

    /* Determinant from LU = direct determinant */
    double det_lu = lu_determinant(LU, piv);
    double det_dir = matrix_determinant(A);
    TEST(CLOSE(det_lu, det_dir, 1e-8), "det(LU) equals det(A)");

    /* QR via MGS */
    double Qv[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    Matrix *Qtest = matrix_create_from_array(Qv, 3, 3);
    Matrix *Q = matrix_create(3,3), *R = matrix_create(3,3);
    ret = qr_decompose_mgs(Qtest, Q, R);
    TEST(ret == 0, "QR decompose (MGS)");
    Matrix *QRprod = matrix_create(3,3);
    matrix_multiply(Q, R, QRprod);
    TEST(matrix_equal(QRprod, Qtest, 1e-8), "A = QR");

    /* Cholesky */
    double SPDv[9] = {4, 2, 2, 2, 5, 3, 2, 3, 6};
    Matrix *SPD = matrix_create_from_array(SPDv, 3, 3);
    Matrix *Lchol = matrix_create_copy(SPD);
    ret = cholesky_decompose(Lchol);
    TEST(ret == 0, "Cholesky decompose");
    Matrix *Lt = matrix_create(3,3);
    matrix_transpose(Lchol, Lt);
    Matrix *LLt = matrix_create(3,3);
    matrix_multiply(Lchol, Lt, LLt);
    TEST(matrix_equal(LLt, SPD, 1e-8), "A = LL^T");

    /* LDL^T */
    double SymV[9] = {1, 2, 3, 2, 2, 1, 3, 1, 3};
    Matrix *Sym = matrix_create_from_array(SymV, 3, 3);
    Vector *Dvec = vector_create(3);
    Matrix *Ldl = matrix_create_copy(Sym);
    ret = ldl_decompose(Ldl, Dvec);
    TEST(ret == 0, "LDL^T decompose");

    /* SVD */
    double Sv[9] = {3, 2, 2, 2, 3, -2, 2, -2, 3};
    Matrix *Smat = matrix_create_from_array(Sv, 3, 3);
    Matrix *Usvd = matrix_create(3,3), *Vt = matrix_create(3,3);
    Vector *S = vector_create(3);
    Matrix *Smat_copy = matrix_create_copy(Smat);
    ret = svd_jacobi(Smat_copy, Usvd, S, Vt, 1e-8, 100);
    TEST(ret >= 0, "SVD converged");
    if (ret >= 0) {
        double cond = svd_condition_number(S);
        TEST(cond > 0, "condition number positive");
        int rnk = svd_effective_rank(S, 1e-8);
        TEST(rnk == 3, "rank = 3");
    }

    /* Forward/backward substitution */
    Matrix *U3 = matrix_create(3,3);
    for (int i = 0; i < 3; i++)
        for (int j = i; j < 3; j++)
            matrix_set(U3, i, j, (i==j)?1.0:0.5);
    Vector *b = vector_create_from_array((double[]){2, 3, 1}, 3);
    Vector *x = vector_create(3);
    solve_backward_substitution(U3, b, x);
    double resid = linear_system_residual(U3, x, b);
    TEST(resid < 1e-8, "backward substitution residual");

    /* Cleanup */
    matrix_free(A); matrix_free(LU); matrix_free(L); matrix_free(U);
    matrix_free(LUprod); matrix_free(PA); matrix_free(Qtest);
    matrix_free(Q); matrix_free(R); matrix_free(QRprod);
    matrix_free(SPD); matrix_free(Lchol); matrix_free(Lt); matrix_free(LLt);
    matrix_free(Sym); matrix_free(Ldl);
    matrix_free(Smat); matrix_free(Smat_copy); matrix_free(Usvd); matrix_free(Vt);
    vector_free(S); vector_free(Dvec);
    matrix_free(U3); vector_free(b); vector_free(x);

    printf("  %d/%d tests passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}