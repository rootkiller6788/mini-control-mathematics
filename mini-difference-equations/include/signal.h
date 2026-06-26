/**
 * signal.h — Discrete-Time Signal Types and Generators
 *
 * Domain: Fundamental discrete-time signals used as building blocks
 * for difference equation analysis and digital signal processing.
 *
 * Knowledge Level Coverage:
 *   L1: Standard signals (impulse, step, ramp, exponential, sinusoid)
 *   L3: Signal metrics (energy, power, RMS, peak)
 *
 * Reference: Oppenheim & Schafer "Discrete-Time Signal Processing" (2010)
 * Reference: Proakis & Manolakis "Digital Signal Processing" (2013)
 *
 * Course Mapping:
 *   MIT 6.003 — Signals and Systems (discrete signals)
 *   Stanford EE 264 — Digital Signal Processing
 *   Berkeley EE 123 — Digital Signal Processing
 */

#ifndef SIGNAL_H
#define SIGNAL_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── L1: Signal Types ────────────────────────────────────────────────── */

/**
 * Window function types for FIR filter design.
 */
typedef enum {
    WINDOW_RECTANGULAR,    /* w[n] = 1 */
    WINDOW_HANN,           /* w[n] = 0.5(1 - cos(2πn/(N-1))) */
    WINDOW_HAMMING,        /* w[n] = 0.54 - 0.46 cos(2πn/(N-1)) */
    WINDOW_BLACKMAN,       /* w[n] = 0.42 - 0.5cos(2πn/(N-1)) + 0.08cos(4πn/(N-1)) */
    WINDOW_BARTLETT,       /* Triangular window */
    WINDOW_KAISER          /* Kaiser window with β parameter */
} WindowType;

/**
 * Filter type for design.
 */
typedef enum {
    FILTER_LOWPASS,
    FILTER_HIGHPASS,
    FILTER_BANDPASS,
    FILTER_BANDSTOP
} FilterType;

/**
 * FIR filter specification.
 */
typedef struct {
    FilterType type;
    int        order;           /* Filter order N */
    double     f_cutoff_low;    /* Lower cutoff (normalized, 0 to 0.5) */
    double     f_cutoff_high;   /* Upper cutoff (for bandpass/bandstop) */
    WindowType window;          /* Window function */
    double     kaiser_beta;     /* Kaiser β (only if WINDOW_KAISER) */
    double    *coefficients;    /* Filter coefficients h[0..N] */
} FIRFilter;

/**
 * IIR filter specification (second-order sections / direct form).
 */
typedef struct {
    FilterType type;
    int        order;           /* Filter order */
    double     f_cutoff;        /* Cutoff frequency (normalized) */
    double    *b_coeffs;        /* Feedforward coeffs b[0..order] */
    double    *a_coeffs;        /* Feedback coeffs a[1..order] (a[0]=1) */
    double     sampling_freq;   /* Sampling frequency in Hz */
} IIRFilter;

/**
 * Moving average filter specification.
 */
typedef struct {
    int     window_size;   /* Number of samples to average */
    double *buffer;        /* Circular buffer of past samples */
    int     buffer_index;  /* Current write position in buffer */
    double  sum;           /* Running sum for O(1) update */
} MovingAvgFilter;

/**
 * Exponential moving average (EMA) / first-order IIR low-pass.
 * y[n] = α·x[n] + (1-α)·y[n-1]
 * Time constant τ ≈ -T/ln(1-α) samples.
 */
typedef struct {
    double  alpha;         /* Smoothing factor 0 < α ≤ 1 */
    double  y_prev;        /* Previous output y[n-1] */
    double  sampling_rate; /* Hz (for time constant computation) */
} ExpMovingAvg;

/**
 * Digital resonator / second-order oscillator.
 * y[n] = 2r·cos(ω₀)·y[n-1] - r²·y[n-2] + x[n]
 */
typedef struct {
    double  r;             /* Pole radius (0 < r < 1 for stable) */
    double  omega0;        /* Resonant frequency (rad/sample) */
    double  y1;            /* y[n-1] */
    double  y2;            /* y[n-2] */
} DigitalResonator;

/* ── API: L1 — Signal Generation ─────────────────────────────────────── */

/**
 * Generate unit impulse δ[n]:
 *   δ[n] = 1 for n = 0, 0 otherwise.
 * Complexity: O(N)
 */
void signal_impulse(double *out, int N, int delay);

/**
 * Generate unit step u[n]:
 *   u[n] = 0 for n < 0, 1 for n ≥ 0.
 * Complexity: O(N)
 */
void signal_step(double *out, int N, int start);

/**
 * Generate ramp r[n] = n for n ≥ 0.
 * Complexity: O(N)
 */
void signal_ramp(double *out, int N);

/**
 * Generate real exponential a^n · u[n].
 * Complexity: O(N)
 */
void signal_exponential(double *out, int N, double base);

/**
 * Generate complex exponential e^{jωn} (stored as interleaved real/imag).
 * Complexity: O(N)
 */
void signal_complex_exponential(double *out_real, double *out_imag,
                                int N, double omega);

/**
 * Generate sinusoidal signal A·cos(ωn + φ).
 * Complexity: O(N)
 */
void signal_sinusoid(double *out, int N,
                     double amplitude, double omega, double phase);

/**
 * Generate white Gaussian noise (Box-Muller method).
 * Complexity: O(N)
 */
void signal_white_noise(double *out, int N, double mean, double stddev);

/**
 * Generate chirp signal: A·cos(ω₀ n + 0.5·k·n²).
 * Linear frequency sweep from ω₀ to ω₀ + k·N.
 * Complexity: O(N)
 */
void signal_chirp(double *out, int N,
                  double amplitude, double omega0, double chirp_rate);

/* ── API: L3 — Signal Metrics ────────────────────────────────────────── */

/**
 * Compute signal energy: E = Σ |x[n]|²
 * Complexity: O(N)
 */
double signal_energy(const double *x, int N);

/**
 * Compute signal average power: P = (1/N) Σ |x[n]|²
 * Complexity: O(N)
 */
double signal_power(const double *x, int N);

/**
 * Compute RMS value: sqrt((1/N) Σ x[n]²)
 * Complexity: O(N)
 */
double signal_rms(const double *x, int N);

/**
 * Find peak value and its index.
 * Complexity: O(N)
 */
double signal_peak(const double *x, int N, int *peak_index);

/**
 * Compute signal mean (DC component).
 * Complexity: O(N)
 */
double signal_mean(const double *x, int N);

/**
 * Compute signal variance.
 * Complexity: O(N)
 */
double signal_variance(const double *x, int N);

/* ── API: L1/L2 — Window Functions ───────────────────────────────────── */

/**
 * Generate window function coefficients w[0..N-1].
 * Complexity: O(N)
 */
void window_generate(double *w, int N, WindowType type, double kaiser_beta);

/* ── API: L6 — Filter Building Blocks ────────────────────────────────── */

/**
 * Initialize moving average filter.
 * Complexity: O(window_size)
 */
void moving_avg_init(MovingAvgFilter *maf, int window_size);

/**
 * Process one sample through moving average filter.
 * Complexity: O(1) — uses running sum
 */
double moving_avg_update(MovingAvgFilter *maf, double x);

/**
 * Initialize exponential moving average.
 * Time constant τ = -T/ln(1-α) → α = 1 - e^{-T/τ}
 * Complexity: O(1)
 */
void exp_moving_avg_init(ExpMovingAvg *ema, double alpha, double fs);

/**
 * Process one sample through EMA filter.
 * Complexity: O(1)
 */
double exp_moving_avg_update(ExpMovingAvg *ema, double x);

/**
 * Initialize digital resonator.
 * Complexity: O(1)
 */
void digital_resonator_init(DigitalResonator *res,
                            double r, double omega0);

/**
 * Process one sample through digital resonator.
 * Complexity: O(1)
 */
double digital_resonator_update(DigitalResonator *res, double x);

/**
 * Apply FIR filter to entire signal (offline).
 * y[n] = Σ_{k=0}^{M} h[k] · x[n-k]
 * Complexity: O(N·M)
 */
void fir_filter_apply(const double *h, int M,
                      const double *x, int N,
                      double *y);

/**
 * Apply IIR filter to entire signal (offline, direct form I).
 * y[n] = Σ b_k x[n-k] - Σ a_k y[n-k]
 * Complexity: O(N·max(order_num, order_den))
 */
void iir_filter_apply(const double *b, int num_order,
                      const double *a, int den_order,
                      const double *x, int N,
                      double *y);

/**
 * Free FIR/IIR filter resources.
 */
/* L6: FIR filter design */
int design_iir_lowpass(double cutoff_freq, double sampling_freq, int order, IIRFilter *filter);
int design_fir_filter(FilterType type, int order, double fc_low, double fc_high, WindowType win, FIRFilter *fir);
void fir_filter_free(FIRFilter *fir);
void iir_filter_free(IIRFilter *iir);
void moving_avg_free(MovingAvgFilter *maf);

#ifdef __cplusplus
}
#endif

#endif /* SIGNAL_H */
