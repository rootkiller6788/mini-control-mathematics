/**
 * @example pole_placement_demo.c
 * @brief Pole placement via Ackermann and Bass-Gura
 * Demonstrates L6: Canonical pole placement problem
 */
#include <stdio.h>
#include "matrix.h"
#include "vector.h"
#include "control_linalg.h"
#include "eigen.h"

int main(void) {
    printf("=== Pole Placement Demo ===\n\n");
    Matrix *A = matrix_create_from_array((double[]){0,1,0,0},2,2);
    Matrix *B = matrix_create_from_array((double[]){0,1},2,1);
    Matrix *K = matrix_create(1,2);
    Vector *poles = vector_create_from_array((double[]){-2,-3},2);

    printf("System: double integrator (open-loop poles: 0,0)\n");
    printf("Desired poles: -2, -3\n\n");

    if (pole_placement_ackermann(A,B,poles,K)==0) {
        printf("Ackermann K = [%.4f %.4f]\n",matrix_get(K,0,0),matrix_get(K,0,1));
        Matrix *Acl=matrix_create(2,2), *BK=matrix_create(2,2);
        matrix_multiply(B,K,BK); matrix_subtract(A,BK,Acl);
        Vector *re=vector_create(2),*im=vector_create(2);
        eigen_qr_general(Acl,re,im,1e-10,100);
        printf("Achieved poles: %.4f, %.4f\n",vector_get(re,0),vector_get(re,1));
        matrix_free(Acl); matrix_free(BK); vector_free(re); vector_free(im);
    }

    if (pole_placement_bass_gura(A,B,poles,K)==0)
        printf("Bass-Gura  K = [%.4f %.4f]\n",matrix_get(K,0,0),matrix_get(K,0,1));

    printf("\nAggressive poles: -5, -6\n");
    poles->data[0]=-5; poles->data[1]=-6;
    if (pole_placement_ackermann(A,B,poles,K)==0)
        printf("Ackermann K = [%.4f %.4f]\n",matrix_get(K,0,0),matrix_get(K,0,1));

    matrix_free(A); matrix_free(B); matrix_free(K); vector_free(poles);
    printf("\nDone.\n");
    return 0;
}