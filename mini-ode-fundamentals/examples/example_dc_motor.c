/**
 * @example example_dc_motor.c
 * @brief  End-to-end example: DC motor step response and transfer function.
 *
 * Demonstrates:
 *   - DC motor state-space model simulation
 *   - Transfer function extraction for control design
 *   - Steady-state analysis
 *
 * Reference: Ogata "Modern Control Engineering" (2010), Chap. 4
 *            Franklin, Powell, Emami-Naeini "Feedback Control" (2019)
 *
 * Application: Servo motor control (L7).
 */

#include "ode_types.h"
#include "ode_ivp.h"
#include "ode_numerical.h"
#include "ode_nonlinear.h"
#include "ode_applications.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(void) {
    printf("=== DC Motor Step Response Analysis ===\n\n");

    /* Typical small DC motor parameters (Pittman 9000 series approx.) */
    DCMotorParams motor = {
        .J = 1.0e-5,     /* Rotor inertia [kg·m²] */
        .b = 1.0e-6,     /* Viscous friction [N·m·s] */
        .K = 0.02,       /* Torque/EMF constant [N·m/A] */
        .L = 2.0e-4,     /* Armature inductance [H] */
        .R = 0.5,        /* Armature resistance [Ω] */
        .V = 12.0        /* Applied voltage [V] */
    };

    printf("Motor parameters:\n");
    printf("  J = %.2e kg·m²,  b = %.2e N·m·s\n", motor.J, motor.b);
    printf("  K = %.2f N·m/A,  L = %.2e H,  R = %.2f Ω\n", motor.K, motor.L, motor.R);
    printf("  V = %.1f V\n\n", motor.V);

    /* Transfer function analysis */
    double num[1], den[3];
    ode_dc_motor_transfer_function(&motor, num, den);
    printf("Transfer function G(s) = omega(s)/V(s):\n");
    printf("  G(s) = %.4e / (%.4e·s² + %.4e·s + %.4e)\n\n", num[0], den[0], den[1], den[2]);

    /* Characteristic equation: den[0]s² + den[1]s + den[2] = 0 */
    double a = den[0], b = den[1], c_val = den[2];
    double disc = b*b - 4.0*a*c_val;
    printf("Poles of the system:\n");
    if (disc >= 0) {
        double p1 = (-b + sqrt(disc)) / (2.0*a);
        double p2 = (-b - sqrt(disc)) / (2.0*a);
        printf("  s₁ = %.2f,  s₂ = %.2f (both real, stable)\n", p1, p2);
    } else {
        double re = -b / (2.0*a);
        double im = sqrt(-disc) / (2.0*a);
        printf("  s₁,₂ = %.2f ± j%.2f (complex, stable)\n", re, im);
    }

    /* Time domain simulation */
    double t_end = 0.05;
    int n_steps = 200;
    double *theta = (double *)malloc((size_t)n_steps * sizeof(double));
    double *omega = (double *)malloc((size_t)n_steps * sizeof(double));
    double *current = (double *)malloc((size_t)n_steps * sizeof(double));
    double *t_vals = (double *)malloc((size_t)n_steps * sizeof(double));

    int ret = ode_dc_motor_simulate(&motor, t_end, n_steps, theta, omega, current, t_vals);
    if (ret != 0) {
        printf("Simulation failed!\n");
        free(theta); free(omega); free(current); free(t_vals);
        return 1;
    }

    /* Steady-state values (from DC analysis) */
    double omega_ss = motor.K * motor.V / (motor.b * motor.R + motor.K * motor.K);
    double i_ss = motor.V / motor.R - motor.K * omega_ss / motor.R;
    printf("\nSteady-state (t → ∞):\n");
    printf("  omega_ss = %.4f rad/s  (%.2f RPM)\n",
           omega_ss, omega_ss * 60.0 / (2.0 * M_PI));
    printf("  i_ss     = %.4f A\n", i_ss);

    /* Print simulation trajectory */
    printf("\nTime response:\n");
    printf("  t (ms)    theta(rad)   omega(rad/s)  current(A)\n");
    int step = n_steps / 10;
    for (int i = 0; i < n_steps; i += step) {
        printf("  %6.2f    %+.6f   %+.6f   %+.6f\n",
               t_vals[i]*1000.0, theta[i], omega[i], current[i]);
    }
    printf("  ... (%d total points)\n", n_steps);

    /* Verify steady state */
    printf("\nFinal values: omega=%.4f (expected %.4f), current=%.4f (expected %.4f)\n",
           omega[n_steps-1], omega_ss, current[n_steps-1], i_ss);

    free(theta); free(omega); free(current); free(t_vals);
    return 0;
}
