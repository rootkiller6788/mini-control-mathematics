#include "complex_control.h"
#include "complex_number.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    printf("=== PID Controller Design for DC Motor ===\n\n");
    TransferFunction plant;
    plant.gain = 10.0;
    plant.num.degree = 0;
    plant.num.coeffs = (Complex *)calloc(1, sizeof(Complex));
    plant.num.coeffs[0] = (Complex){1,0};
    plant.den.degree = 3;
    plant.den.coeffs = (Complex *)calloc(4, sizeof(Complex));
    plant.den.coeffs[0] = (Complex){0,0};
    plant.den.coeffs[1] = (Complex){10,0};
    plant.den.coeffs[2] = (Complex){7,0};
    plant.den.coeffs[3] = (Complex){1,0};
    printf("Plant: G(s) = 10 / (s*(s+2)*(s+5))\n\n");
    FreqRespPoint frp_ol = tf_freq_resp(&plant, 2.0);
    printf("|G(j*2)| = %.4f (%.1f dB), arg = %.1f deg\n\n",
           frp_ol.magnitude, frp_ol.magnitude_db, frp_ol.phase_deg);
    printf("--- PI Controller Design (PM=60deg, wc=2.0) ---\n");
    PIDController pid = design_pi_phase_margin(&plant, 60.0, 2.0, 100);
    printf("Designed PI: Kp=%.4f, Ki=%.4f\n\n", pid.Kp, pid.Ki);
    TransferFunction pid_tf = pid_to_tf(&pid);
    TransferFunction loop_tf = tf_series(&pid_tf, &plant);
    TransferFunction unity_h;
    unity_h.gain = 1.0;
    unity_h.num.degree = 0;
    unity_h.num.coeffs = (Complex *)calloc(1, sizeof(Complex));
    unity_h.num.coeffs[0] = (Complex){1,0};
    unity_h.den.degree = 0;
    unity_h.den.coeffs = (Complex *)calloc(1, sizeof(Complex));
    unity_h.den.coeffs[0] = (Complex){1,0};
    TransferFunction cl_tf = tf_feedback(&loop_tf, &unity_h);
    StabilityMargins sm = tf_stability_margins(&loop_tf, 0.1, 10.0, 200);
    printf("--- Closed-loop Performance ---\n");
    printf("Gain Margin: %.2f dB\n", sm.gain_margin_db);
    printf("Phase Margin: %.2f deg\n", sm.phase_margin_deg);
    printf("Gain Crossover: %.4f rad/s\n", sm.gain_crossover);
    int stable = tf_is_stable(&cl_tf, 1e-6);
    printf("Closed-loop stable: %s\n", stable ? "YES" : "NO");
    double dc_gain = tf_dc_gain(&cl_tf);
    printf("DC gain: %.4f\n", dc_gain);
    double bw = tf_bandwidth(&cl_tf, 0.1, 10.0, 100);
    printf("-3dB Bandwidth: %.4f rad/s\n", bw);
    free(plant.num.coeffs); free(plant.den.coeffs);
    free(pid_tf.num.coeffs); free(pid_tf.den.coeffs);
    free(loop_tf.num.coeffs); free(loop_tf.den.coeffs);
    free(cl_tf.num.coeffs); free(cl_tf.den.coeffs);
    free(unity_h.num.coeffs); free(unity_h.den.coeffs);
    printf("\n=== Design Complete ===\n");
    return 0;
}
