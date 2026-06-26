#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/numerical_core.h"
#include "../include/eigen_methods.h"

int main() {
    double Ad[] = {2, -1, -1, 2};
    Matrix *A = matrix_create_from(2, 2, Ad);
    Vector *v = vector_create(2);
    v->data[0] = 1.0; v->data[1] = -0.3;
    double lambda;
    int ret = power_iteration(A, &lambda, v, 100, 1e-10);
    printf("Return: %d, lambda = %.15f\n", ret, lambda);
    printf("Eigenvector: [%.10f, %.10f]\n", v->data[0], v->data[1]);
    printf("Check: A*v = [%.10f, %.10f]\n",
           2*v->data[0] - v->data[1],
           -v->data[0] + 2*v->data[1]);
    printf("Check: lambda*v = [%.10f, %.10f]\n",
           lambda*v->data[0], lambda*v->data[1]);
    matrix_free(A); vector_free(v);
    return 0;
}
