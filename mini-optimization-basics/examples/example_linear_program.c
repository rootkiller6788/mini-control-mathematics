/**
 * example_linear_program.c — LP Example: Resource Allocation
 *
 * Demonstrates solving a linear program using the simplex method.
 *
 * Problem: A factory produces two products (x1, x2).
 *   Profit: $40 per x1, $30 per x2
 *   Machine A: 2*x1 + x2 <= 100 hours
 *   Machine B: x1 + x2 <= 80 hours
 *   Machine C: x1 <= 40 units
 *   Non-negativity: x1, x2 >= 0
 *
 * Optimal: x1 = 20, x2 = 60, profit = 40*20 + 30*60 = 2600
 */

#include "constrained_optimization.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

static double lin_eval(const convex_function_t *f, const vector_t *x)
{
    (void)f;
    return x->data[0] + x->data[1];
}

static void lin_grad(const convex_function_t *f, const vector_t *x, vector_t *g)
{
    (void)f; (void)x;
    g->data[0] = 1.0;
    g->data[1] = 1.0;
}

int main(void)
{
    printf("=== Linear Programming: Resource Allocation ===\n\n");

    /* Set up LP in standard form:
     * minimize -40*x1 - 30*x2
     * s.t. 2*x1 + x2 + s1 = 100
     *       x1 + x2 + s2 = 80
     *       x1 + s3 = 40
     *       x1, x2, s1, s2, s3 >= 0
     */

    linear_program_t lp;
    memset(&lp, 0, sizeof(lp));
    lp.n = 5;  /* x1, x2, s1, s2, s3 */
    lp.m = 3;

    /* Cost vector (negative for maximization) */
    lp.c.n = 5;
    double cdat[] = {-40.0, -30.0, 0.0, 0.0, 0.0};
    lp.c.data = cdat;

    /* Constraint matrix A (3x5) */
    lp.A = *matrix_alloc(3, 5);
    double adat[] = {
        2.0, 1.0, 1.0, 0.0, 0.0,
        1.0, 1.0, 0.0, 1.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 1.0
    };
    memcpy(lp.A.data, adat, 15 * sizeof(double));

    /* RHS */
    lp.b.n = 3;
    double bdat[] = {100.0, 80.0, 40.0};
    lp.b.data = bdat;

    /* Initialize result */
    opt_result_t result;
    memset(&result, 0, sizeof(result));

    simplex_config_t cfg = simplex_config_default();
    cfg.verbose = 0;

    simplex_solve(&lp, &cfg, &result);

    printf("Simplex result:\n");
    printf("  Status: %d (0=OPTIMAL)\n", result.status);
    printf("  Iterations: %zu\n", result.iterations);
    printf("  Optimal value: %.4f (profit = %.4f)\n", result.f_value, -result.f_value);
    printf("  x1 = %.4f, x2 = %.4f\n", result.x.data[0], result.x.data[1]);
    printf("  Slack: s1=%.4f, s2=%.4f, s3=%.4f\n",
           result.x.data[2], result.x.data[3], result.x.data[4]);

    /* Expected: x1=20, x2=60, profit=2600 */
    double profit = -result.f_value;
    printf("\nExpected: x1=20, x2=60, profit=2600\n");
    printf("Actual:   x1=%.1f, x2=%.1f, profit=%.0f\n",
           result.x.data[0], result.x.data[1], profit);

    free(result.x.data);
    matrix_free(&lp.A);

    printf("\n=== LP Example Complete ===\n");
    return 0;
}
