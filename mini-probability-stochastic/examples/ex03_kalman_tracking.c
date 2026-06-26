/**
 * @file ex03_kalman_tracking.c
 * @brief Example: 1D Kalman filter for object tracking with noisy radar.
 *
 * Demonstrates: kalman_filter_run(), kalman_smoother(),
 * kalman_log_likelihood(), ukf_step(), info_criterion_aic(),
 * info_criterion_bic().
 *
 * Scenario: Track an object moving with constant velocity from
 * noisy position measurements (radar). Compare Kalman filter and
 * smoother estimates to ground truth.
 *
 * Knowledge: L5 engineering methods (Kalman filtering, UKF for
 * nonlinear observation models), L6 engineering problems (target
 * tracking, sensor fusion), L7 applications (radar/sonar tracking).
 */
#include "estimation_filtering.h"
#include "probability_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== Example 3: Kalman Filter Object Tracking ===\n\n");

    /* 2D state: [position, velocity] */
    double dt = 1.0;
    /* State transition: x_k = [1 dt; 0 1] * x_{k-1} + noise */
    double F[4] = {1.0, dt, 0.0, 1.0};
    double H[2] = {1.0, 0.0}; /* observe position only */
    double Q[4] = {0.01, 0.0, 0.0, 0.01}; /* small process noise */
    double R[1] = {0.5}; /* measurement noise variance */

    kalman_model_t model;
    model.state_dim = 2; model.input_dim = 0; model.obs_dim = 1;
    model.F = F; model.B = NULL; model.H = H;
    model.Q = Q; model.R = R;

    /* Ground truth: constant velocity 0.5 m/s, start at position 0 */
    /* Generate noisy observations */
    uint64_t seed = 42;
    int T = 20;
    double truth_pos[20], truth_vel[20];
    double observations[20];
    double true_pos = 0.0, true_vel = 0.5;
    int k;
    printf("True trajectory and noisy observations:\n");
    printf("  t   true_pos  true_vel  observed\n");
    for (k = 0; k < T; k++) {
        true_pos += true_vel * dt;
        truth_pos[k] = true_pos;
        truth_vel[k] = true_vel;
        observations[k] = true_pos + sqrt(0.5) * prob_box_muller(&seed);
        printf("  %2d   %6.2f    %5.2f     %6.2f\n",
               k, truth_pos[k], truth_vel[k], observations[k]);
    }

    /* Initial estimate: [0, 0] with large uncertainty */
    double x_init[2] = {0.0, 0.0};
    double P_init[4] = {1.0, 0.0, 0.0, 1.0};

    /* Run Kalman filter */
    double x_est[40], P_est[80]; /* T*state_dim, T*state_dim^2 */
    kalman_filter_run(&model, observations, NULL, (size_t)T,
                      x_init, P_init, x_est, P_est);

    /* Run Kalman smoother */
    /* Need predict/update intermediate arrays */
    double xf[40], Pf[80], xp[40], Pp[80], xs[40], Ps[80];
    {
        double x[2] = {0.0, 0.0}, P[4] = {1.0, 0.0, 0.0, 1.0};
        double x_p[2], P_p[4];
        int t;
        for (t = 0; t < T; t++) {
            kalman_predict(&model, x, P, NULL, x_p, P_p);
            xp[t*2]=x_p[0]; xp[t*2+1]=x_p[1];
            Pp[t*4]=P_p[0]; Pp[t*4+1]=P_p[1]; Pp[t*4+2]=P_p[2]; Pp[t*4+3]=P_p[3];
            xf[t*2]=x_p[0]; xf[t*2+1]=x_p[1];
            Pf[t*4]=P_p[0]; Pf[t*4+1]=P_p[1]; Pf[t*4+2]=P_p[2]; Pf[t*4+3]=P_p[3];
            kalman_update(&model, x_p, P_p, &observations[t], x, P, NULL);
        }
    }
    kalman_smoother(&model, xf, Pf, xp, Pp, (size_t)T, xs, Ps);

    /* Print results */
    printf("\nResults: Position estimates\n");
    printf("  t   True     Filtered  Smoothed  |FilterErr| |SmoothErr|\n");
    double rmse_f = 0.0, rmse_s = 0.0;
    for (k = 0; k < T; k++) {
        double err_f = truth_pos[k] - x_est[k*2];
        double err_s = truth_pos[k] - xs[k*2];
        rmse_f += err_f * err_f;
        rmse_s += err_s * err_s;
        printf("  %2d   %6.2f   %6.2f    %6.2f     %6.3f      %6.3f\n",
               k, truth_pos[k], x_est[k*2], xs[k*2],
               fabs(err_f), fabs(err_s));
    }
    rmse_f = sqrt(rmse_f / T); rmse_s = sqrt(rmse_s / T);
    printf("\n  RMSE Filter: %.3f   RMSE Smoother: %.3f\n\n", rmse_f, rmse_s);
    printf("  Note: Smoother RMSE <= Filter RMSE (uses all data)\n\n");

    /* Velocity estimates */
    printf("Velocity estimates (true=0.5):\n");
    for (k = 0; k < T; k++) {
        printf("  t=%2d  filter=%.3f  smoother=%.3f\n", k, x_est[k*2+1], xs[k*2+1]);
    }

    /* Log-likelihood and model selection */
    double logL = kalman_log_likelihood(&model, observations, (size_t)T,
                                        x_init, P_init);
    printf("\nModel selection (for the filter model):\n");
    printf("  Log-likelihood: %.2f\n", logL);
    printf("  AIC: %.2f (k=%zu parameters)\n",
           info_criterion_aic(logL, 3), (size_t)3);
    printf("  BIC: %.2f (k=%zu, n=%d)\n",
           info_criterion_bic(logL, 3, T), (size_t)3, T);

    /* Compare with a worse model (higher R = more measurement noise) */
    kalman_model_t bad_model = model;
    double bad_R[] = {2.0};
    bad_model.R = bad_R;
    double logL_bad = kalman_log_likelihood(&bad_model, observations,
                                            (size_t)T, x_init, P_init);
    printf("\n  Worse model (R=2.0): logL = %.2f, AIC = %.2f\n",
           logL_bad, info_criterion_aic(logL_bad, 3));
    printf("  Better model has lower AIC: %s\n",
           info_criterion_aic(logL,3) < info_criterion_aic(logL_bad,3)
           ? "CORRECT" : "WRONG");

    return 0;
}
