/**
 * @file test_basic.c
 * @brief Tests for matrix/vector operations
 */
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
#include "solvers.h"
#include "norms.h"

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(cond, msg) do { tests_run++; if (cond) { tests_passed++; } else { printf("  FAIL: %s\n", msg); } } while(0)

#define CLOSE(a, b, tol) (fabs((a)-(b)) < (tol))

int main(void) {
    printf("=== test_basic: matrix/vector operations ===\n");

    /* Matrix create/free */
    Matrix *A = matrix_create(3, 3);
    TEST(A != NULL, "matrix_create");
    TEST(A->rows == 3 && A->cols == 3, "matrix dimensions");
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            matrix_set(A, i, j, (double)(i * 3 + j + 1));

    TEST(CLOSE(matrix_get(A, 0, 0), 1.0, 1e-10), "matrix_set/get (0,0)");
    TEST(CLOSE(matrix_get(A, 2, 2), 9.0, 1e-10), "matrix_set/get (2,2)");

    /* Matrix multiply by identity */
    Matrix *B = matrix_create(3, 3);
    for (int i = 0; i < 3; i++) matrix_set(B, i, i, 1.0);
    Matrix *C = matrix_create(3, 3);
    matrix_multiply(A, B, C);
    TEST(CLOSE(matrix_get(C, 0, 0), 1.0, 1e-10), "matrix multiply by identity");

    /* Transpose */
    Matrix *At = matrix_create(3, 3);
    matrix_transpose(A, At);
    TEST(CLOSE(matrix_get(At, 0, 1), 4.0, 1e-10), "matrix transpose");

    /* Trace = 1+5+9 = 15 */
    double tr = matrix_trace(A);
    TEST(CLOSE(tr, 15.0, 1e-10), "matrix trace");

    /* Determinant of tridiagonal */
    double Av2[9] = {2, -1, 0, -1, 2, -1, 0, -1, 2};
    Matrix *Tmat = matrix_create_from_array(Av2, 3, 3);
    double det = matrix_determinant(Tmat);
    TEST(CLOSE(det, 4.0, 1e-8), "determinant of tridiagonal (det=4)");

    /* Identity determinant = 1 */
    Matrix *I3 = matrix_create_identity(3);
    TEST(CLOSE(matrix_determinant(I3), 1.0, 1e-10), "det(I)=1");

    /* Vector dot product */
    double v1v[3] = {1, 2, 3};
    double v2v[3] = {4, 5, 6};
    Vector *v1 = vector_create_from_array(v1v, 3);
    Vector *v2 = vector_create_from_array(v2v, 3);
    TEST(CLOSE(vector_dot(v1, v2), 32.0, 1e-10), "dot product 1*4+2*5+3*6=32");
    TEST(CLOSE(vector_norm_l2(v1), sqrt(14.0), 1e-10), "L2 norm");

    /* Cross product */
    Vector *cross = vector_create(3);
    vector_cross(v1, v2, cross);
    TEST(CLOSE(vector_get(cross, 0), -3.0, 1e-10), "cross product x");
    TEST(CLOSE(vector_get(cross, 1), 6.0, 1e-10), "cross product y");
    TEST(CLOSE(vector_get(cross, 2), -3.0, 1e-10), "cross product z");

    /* Gram-Schmidt */
    Vector *gs1 = vector_create_from_array((double[]){1,0,0}, 3);
    Vector *gs2 = vector_create_from_array((double[]){1,1,0}, 3);
    Vector *gs3 = vector_create_from_array((double[]){1,1,1}, 3);
    Vector *gs_vecs[3] = {gs1, gs2, gs3};
    int ret = vector_gram_schmidt(gs_vecs, 3, 1e-10);
    TEST(ret == 0, "Gram-Schmidt success");
    if (ret == 0) {
        TEST(CLOSE(vector_dot(gs_vecs[0], gs_vecs[1]), 0.0, 1e-10), "GS orth 0-1");
        TEST(CLOSE(vector_dot(gs_vecs[0], gs_vecs[2]), 0.0, 1e-10), "GS orth 0-2");
        TEST(CLOSE(vector_dot(gs_vecs[1], gs_vecs[2]), 0.0, 1e-10), "GS orth 1-2");
    }

    /* Linear system solve */
    double Av3[9] = {2, 1, 0, 1, 2, 1, 0, 1, 2};
    Matrix *A2 = matrix_create_from_array(Av3, 3, 3);
    Vector *b = vector_create_from_array((double[]){1, 2, 3}, 3);
    Vector *x = vector_create(3);
    ret = solve_gaussian_elimination(A2, b, x);
    TEST(ret == 0, "Gaussian elimination");
    if (ret == 0) {
        double res = linear_system_residual(A2, x, b);
        TEST(res < 1e-8, "linear system residual small");
    }

    /* Cleanup */
    matrix_free(A); matrix_free(B); matrix_free(C); matrix_free(At);
    matrix_free(Tmat); matrix_free(I3); matrix_free(A2);
    vector_free(v1); vector_free(v2); vector_free(cross);
    vector_free(gs1); vector_free(gs2); vector_free(gs3);
    vector_free(b); vector_free(x);

    printf("  %d/%d tests passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}