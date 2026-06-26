/**
 * @example controllability_check.c
 * @brief Check controllability and observability of linear systems
 * Demonstrates L6: Canonical controllability/observability problem
 */
#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "control_linalg.h"
#include "norms.h"

int main(void) {
    printf("=== Controllability & Observability Analysis ===\n\n");
    /* Double Integrator */
    printf("1. Double Integrator: A=[0 1;0 0], B=[0;1]\n");
    Matrix *A = matrix_create_from_array((double[]){0,1,0,0},2,2);
    Matrix *B = matrix_create_from_array((double[]){0,1},2,1);
    Matrix *Cmat = matrix_create(2,2);
    control_controllability_matrix(A,B,Cmat);
    printf("   Rank(C)=%d, Controllable=%s\n\n",
           matrix_rank(Cmat,1e-10),
           control_is_controllable(A,B,1e-10)?"YES":"NO");

    /* Inverted Pendulum */
    printf("2. Inverted Pendulum: A=[0 1;g/l 0], B=[0;1/ml^2]\n");
    double l=0.5, m=0.2;
    matrix_free(A); matrix_free(B);
    A=matrix_create_from_array((double[]){0,1,9.81/l,0},2,2);
    B=matrix_create_from_array((double[]){0,1.0/(m*l*l)},2,1);
    control_controllability_matrix(A,B,Cmat);
    printf("   Rank(C)=%d, Controllable=%s, Cond=%.2e\n\n",
           matrix_rank(Cmat,1e-10),
           control_is_controllable(A,B,1e-10)?"YES":"NO",
           controllability_condition_number(A,B));

    /* Uncontrollable */
    printf("3. Uncontrollable: A=[-1 0;0 -2], B=[1;0]\n");
    matrix_free(A); matrix_free(B);
    A=matrix_create_from_array((double[]){-1,0,0,-2},2,2);
    B=matrix_create_from_array((double[]){1,0},2,1);
    printf("   Controllable=%s (expected: NO)\n\n",
           control_is_controllable(A,B,1e-10)?"YES":"NO");

    printf("4. PBH test at lambda=-2: %s (expected: FAIL)\n",
           control_pbh_controllability_test(A,B,-2.0,0.0,1e-10)?
           "PASS" : "FAIL (rank<n)");

    matrix_free(A); matrix_free(B); matrix_free(Cmat);
    printf("\nDone.\n");
    return 0;
}