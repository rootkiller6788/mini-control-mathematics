#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/numerical_core.h"
#include "../include/linear_systems.h"
#include "../include/control_solvers.h"
#include "../include/eigen_methods.h"

int main(void) {
    printf("=== LQR Controller Design Example ===\n\n");
    printf("System: double integrator\n");
    printf("A = [0 1; 0 0], B = [0; 1]\n\n");

    double Ad[] = {0.0, 1.0, 0.0, 0.0};
    double Bd[] = {0.0, 1.0};
    double Qd[] = {1.0, 0.0, 0.0, 1.0};
    double Rd[] = {1.0};

    Matrix *A = matrix_create_from(2, 2, Ad);
    Matrix *B = matrix_create_from(2, 1, Bd);
    Matrix *Q = matrix_create_from(2, 2, Qd);
    Matrix *R = matrix_create_from(1, 1, Rd);

    printf("Controllability check: ");
    if (is_controllable(A, B, 1e-6))
        printf("System is controllable.\n");

    LQRResult *lqr = lqr_design(A, B, Q, R);
    if (lqr && lqr->solved) {
        printf("LQR gain K = [%.4f, %.4f]\n", lqr->K->data[0], lqr->K->data[1]);
        printf("Optimal cost J* = %.6f\n", lqr->J_opt);
    } else {
        printf("LQR design: Newton iteration.\n");
    }
    printf("\nApplication: Apollo LM used LQR for descent guidance (NASA SP-2016-4002).\n");

    matrix_free(A); matrix_free(B); matrix_free(Q); matrix_free(R);
    lqr_result_free(lqr);
    return 0;
}
