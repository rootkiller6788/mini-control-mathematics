#include "complex_number.h"
#include "complex_integration.h"
#include "complex_series.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== Real Integral via Residue Theorem ===\n\n");
    printf("Problem: I = integral_0^{2pi} dtheta / (5 + 4*cos(theta))\n\n");
    printf("Method: Substitute z = e^{i*theta}\n");
    printf("After substitution, poles at z = -1/2 and z = -2\n");
    printf("Only z = -1/2 is inside unit circle.\n\n");

    /* Analytical residue calculation */
    double residue = 1.0 / 3.0;
    double I_analytical = 2.0 * M_PI / 3.0;
    double I_residue = 2.0 * M_PI * residue;

    printf("Residue at z=-1/2: %.6f\n", residue);
    printf("I (residue theorem) = 2*pi * %.6f = %.10f\n", residue, I_residue);
    printf("I (analytical)      = 2*pi/3 = %.10f\n", I_analytical);
    printf("Error: %.2e\n\n", fabs(I_residue - I_analytical));

    /* Numeric verification via Simpson integration */
    int n = 1000;
    double dtheta = 2.0 * M_PI / n;
    double num_sum = 0.0;
    for (int i = 0; i < n; i++) {
        double theta = i * dtheta;
        double f_val = 1.0 / (5.0 + 4.0 * cos(theta));
        if (i == 0 || i == n-1) num_sum += f_val;
        else if (i % 2 == 1) num_sum += 4.0 * f_val;
        else num_sum += 2.0 * f_val;
    }
    double I_simpson = (dtheta / 3.0) * num_sum;
    printf("I (Simpson, n=%d) = %.10f\n", n, I_simpson);
    printf("Error vs analytic: %.2e\n", fabs(I_simpson - I_analytical));
    printf("\n=== Residue Theorem Application Complete ===\n");
    return 0;
}
