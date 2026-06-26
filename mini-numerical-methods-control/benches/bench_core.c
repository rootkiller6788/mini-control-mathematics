#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../include/numerical_core.h"
#include "../include/linear_systems.h"

int main(void) {
    printf("=== Benchmark: Core Numerical Operations ===

");
    size_t n = 200;
    Matrix *A = matrix_create(n, n);
    Matrix *B = matrix_create(n, n);
    for (size_t i = 0; i < n; i++) {
        A->data[i*n+i] = 2.0;
        B->data[i*n+i] = 1.0;
        if (i > 0) A->data[i*n+(i-1)] = -1.0;
        if (i < n-1) A->data[i*n+(i+1)] = -1.0;
    }

    clock_t start, end;
    printf("Matrix multiply (%zux%zu x %zux%zu)... ", n, n, n, n);
    fflush(stdout);
    start = clock();
    Matrix *C = matrix_multiply(A, B);
    end = clock();
    double t_mult = (double)(end - start) / CLOCKS_PER_SEC;
    printf("%.4f s
", t_mult);

    printf("LU decomposition (%zux%zu)... ", n, n);
    fflush(stdout);
    start = clock();
    LUResult *lu = lu_decompose(A);
    end = clock();
    double t_lu = (double)(end - start) / CLOCKS_PER_SEC;
    printf("%.4f s
", t_lu);

    printf("Determinant... ");
    fflush(stdout);
    start = clock();
    double det = matrix_determinant(A);
    end = clock();
    double t_det = (double)(end - start) / CLOCKS_PER_SEC;
    printf("%.4f s (det=%.6e)
", t_det, det);

    printf("
Summary: multiply=%.4fs, LU=%.4fs, det=%.4fs
", t_mult, t_lu, t_det);
    printf("FLOP/s (LU): ~%.2e
", (2.0*n*n*n/3.0) / fmax(t_lu, 1e-6));

    matrix_free(A); matrix_free(B); matrix_free(C);
    lu_free(lu);
    return 0;
}
