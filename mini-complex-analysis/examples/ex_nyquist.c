#include "complex_control.h"
#include "complex_number.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    printf("=== Nyquist Stability Analysis ===\n\n");
    TransferFunction plant;
    plant.gain = 2.0;
    plant.num.degree = 0;
    plant.num.coeffs = (Complex *)calloc(1, sizeof(Complex));
    plant.num.coeffs[0] = (Complex){1,0};
    plant.den.degree = 3;
    plant.den.coeffs = (Complex *)calloc(4, sizeof(Complex));
    plant.den.coeffs[0] = (Complex){1,0};
    plant.den.coeffs[1] = (Complex){3,0};
    plant.den.coeffs[2] = (Complex){3,0};
    plant.den.coeffs[3] = (Complex){1,0};
    printf("Plant: G(s) = 2 / (s+1)^3\n");
    printf("Relative degree: %d\n\n", tf_relative_degree(&plant));
    int ol_stable = tf_is_stable(&plant, 1e-6);
    printf("Open-loop stable: %s\n", ol_stable ? "YES" : "NO");
    Complex poles[3];
    int n_poles = tf_find_poles(&plant, poles, 3, 1e-8, 200);
    printf("Open-loop poles (%d found):\n", n_poles);
    for (int i = 0; i < n_poles; i++)
        printf("  p%d = %.4f + %.4fi\n", i+1, poles[i].re, poles[i].im);
    printf("\n--- Nyquist Analysis ---\n");
    int cl_stable = nyquist_stability_check(&plant, 0.01, 20.0, 200, 10.0);
    printf("Closed-loop stable (Nyquist): %s\n", cl_stable ? "YES" : "NO");
    StabilityMargins sm = tf_stability_margins(&plant, 0.01, 10.0, 200);
    printf("\n--- Stability Margins ---\n");
    printf("Gain Margin: %.2f dB\n", sm.gain_margin_db);
    printf("Phase Margin: %.2f deg\n", sm.phase_margin_deg);
    printf("Phase Crossover: %.4f rad/s\n", sm.phase_crossover);
    printf("Gain Crossover: %.4f rad/s\n", sm.gain_crossover);
    printf("\n--- Frequency Response ---\n");
    double freqs[4] = {0.1, 1.0, 3.0, 10.0};
    for (int i = 0; i < 4; i++) {
        FreqRespPoint frp = tf_freq_resp(&plant, freqs[i]);
        printf("  w=%.1f: |G|=%.4f (%.1f dB), ph=%.1f deg\n",
               freqs[i], frp.magnitude, frp.magnitude_db, frp.phase_deg);
    }
    free(plant.num.coeffs); free(plant.den.coeffs);
    printf("\n=== Analysis Complete ===\n");
    return 0;
}
