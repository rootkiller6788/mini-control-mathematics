#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/numerical_core.h"
#include "../include/eigen_methods.h"
#include "../include/control_solvers.h"

int main(void) {
    printf("=== Pole Placement: Ackermann Formula ===\n\n");

    double Ad[] = {0.0, 1.0, 0.0, 0.0};
    double Bd[] = {0.0, 1.0};
    Matrix *A = matrix_create_from(2, 2, Ad);
    Matrix *B = matrix_create_from(2, 1, Bd);

    printf("Open-loop eigenvalues: ");
    EigenResult *eig_open = qr_algorithm_eigenvalues(A);
    if (eig_open) {
        for (size_t i = 0; i < eig_open->n; i++)
            printf("%.2f ", eig_open->real_parts[i]);
        printf("\n");
        eigen_result_free(eig_open);
    }

    printf("Desired poles: -2, -3 (settling time ~2s)\n");
    double desired_poles[] = {-2.0, -3.0};

    Matrix *K = pole_placement_acker(A, B, desired_poles, 2);
    if (K) {
        printf("Feedback gain K = [%.4f, %.4f]\n", K->data[0], K->data[1]);

        Matrix *BK = matrix_multiply(B, K);
        Matrix *Acl = matrix_sub(A, BK);
        EigenResult *eig_cl = qr_algorithm_eigenvalues(Acl);
        if (eig_cl) {
            SpectralRadius *sr = spectral_radius_compute(eig_cl);
            printf("Closed-loop is Hurwitz: %s\n", sr->is_hurwitz ? "yes" : "no");
            spectral_radius_free(sr);
            eigen_result_free(eig_cl);
        }
        matrix_free(BK); matrix_free(Acl);
        matrix_free(K);
    } else {
        printf("Pole placement failed.\n");
    }

    printf("\nApplication: Boeing 747 lateral-directional control uses eigenvalue placement.\n");
    matrix_free(A); matrix_free(B);
    return 0;
}
