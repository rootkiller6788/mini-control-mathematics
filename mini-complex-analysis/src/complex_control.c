/**
 * @file complex_control.c
 * @brief Control theory via complex analysis: transfer functions,
 *        Nyquist stability, root locus, Bode analysis, PID design
 *
 * L1-L7: From transfer functions through Nyquist to PID controller design
 *
 * Key theorems and methods:
 *   Nyquist Stability Criterion (1932): Z = N + P
 *   Argument Principle applied to feedback stability
 *   Root Locus method (Evans, 1948)
 *   Bode analysis (gain/phase margins)
 *   PID tuning via frequency-domain methods
 *
 * Course alignment:
 *   MIT 6.302 -- Feedback Systems
 *   Stanford AA 278 -- Control System Design
 *   Cambridge 3F2 -- Control Systems Design
 *   Tsinghua -- Automatic Control Theory
 *
 * Textbook: Ogata (2010), Dorf & Bishop (2016), Astrom & Murray (2010)
 */

#include "complex_control.h"
#include "complex_integration.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

Complex tf_eval(const TransferFunction *tf, Complex s) {
    if (!tf) return (Complex){0,0};
    Complex num_val = complex_poly_eval(&tf->num, s);
    Complex den_val = complex_poly_eval(&tf->den, s);
    Complex ratio = complex_div(num_val, den_val);
    Complex result = complex_scale(ratio, tf->gain);
    return result;
}

/* Polynomial root finding via companion matrix eigenvalues (simplified Newton) */
static int poly_roots_newton(const ComplexPoly *p, Complex *roots,
                              int max_roots, double tol, int max_iter) {
    if (!p || !roots || max_roots <= 0) return 0;
    int found = 0;
    int deg = p->degree;
    if (deg <= 0) return 0;
    /* Try initial guesses on circles of increasing radius */
    for (int r_idx = 0; r_idx < 3 && found < max_roots && found < deg; r_idx++) {
        double radius = 0.5 + r_idx * 1.0;
        for (int ang = 0; ang < deg && found < max_roots && found < deg; ang++) {
            double theta = 2.0 * M_PI * ang / deg;
            Complex z = {radius * cos(theta), radius * sin(theta)};
            for (int it = 0; it < max_iter; it++) {
                Complex fz = complex_poly_eval(p, z);
                double abs_fz = complex_abs(fz);
                if (abs_fz < tol) {
                    /* Check if this root is already found */
                    int dup = 0;
                    for (int j = 0; j < found; j++) {
                        if (complex_distance(z, roots[j]) < tol * 10) { dup = 1; break; }
                    }
                    if (!dup) { roots[found++] = z; }
                    break;
                }
                /* Compute derivative */
                Complex dfz = {0,0};
                for (int k = 1; k <= deg; k++) {
                    Complex kc = {(double)k, 0};
                    Complex zkm1 = complex_pow_int(z, k-1);
                    dfz = complex_add(dfz, complex_mul(complex_mul(kc, p->coeffs[k]), zkm1));
                }
                double d2 = complex_abs2(dfz);
                if (d2 < 1e-30) break;
                z = complex_sub(z, complex_div(fz, dfz));
            }
        }
    }
    return found;
}

int tf_find_poles(const TransferFunction *tf, Complex *poles,
                  int max_poles, double tol, int max_iter) {
    if (!tf || !poles || max_poles <= 0) return 0;
    return poly_roots_newton(&tf->den, poles, max_poles, tol, max_iter);
}

int tf_find_zeros(const TransferFunction *tf, Complex *zeros,
                  int max_zeros, double tol, int max_iter) {
    if (!tf || !zeros || max_zeros <= 0) return 0;
    return poly_roots_newton(&tf->num, zeros, max_zeros, tol, max_iter);
}

int tf_is_stable(const TransferFunction *tf, double stability_margin) {
    if (!tf) return 0;
    int deg = tf->den.degree;
    if (deg <= 0) return 1;
    Complex *poles = (Complex *)calloc(deg, sizeof(Complex));
    if (!poles) return 0;
    int n_poles = poly_roots_newton(&tf->den, poles, deg, 1e-8, 200);
    int stable = 1;
    for (int i = 0; i < n_poles; i++) {
        if (poles[i].re > -stability_margin) { stable = 0; break; }
    }
    free(poles);
    return stable;
}

int tf_relative_degree(const TransferFunction *tf) {
    if (!tf) return 0;
    return tf->den.degree - tf->num.degree;
}


/* Polynomial multiplication: C(z) = A(z) * B(z) */
static ComplexPoly poly_multiply(const ComplexPoly *a, const ComplexPoly *b) {
    ComplexPoly result;
    result.degree = a->degree + b->degree;
    result.coeffs = (Complex *)calloc(result.degree + 1, sizeof(Complex));
    if (!result.coeffs) { result.degree = -1; return result; }
    for (int i = 0; i <= a->degree; i++) {
        for (int j = 0; j <= b->degree; j++) {
            Complex prod = complex_mul(a->coeffs[i], b->coeffs[j]);
            result.coeffs[i + j] = complex_add(result.coeffs[i + j], prod);
        }
    }
    return result;
}

/* Polynomial addition */
static ComplexPoly poly_add(const ComplexPoly *a, const ComplexPoly *b) {
    int max_deg = (a->degree > b->degree) ? a->degree : b->degree;
    ComplexPoly result;
    result.degree = max_deg;
    result.coeffs = (Complex *)calloc(max_deg + 1, sizeof(Complex));
    if (!result.coeffs) { result.degree = -1; return result; }
    for (int i = 0; i <= a->degree; i++)
        result.coeffs[i] = complex_add(result.coeffs[i], a->coeffs[i]);
    for (int i = 0; i <= b->degree; i++)
        result.coeffs[i] = complex_add(result.coeffs[i], b->coeffs[i]);
    return result;
}

/* System interconnections */
TransferFunction tf_series(const TransferFunction *g1, const TransferFunction *g2) {
    TransferFunction result;
    result.gain = g1->gain * g2->gain;
    result.num = poly_multiply(&g1->num, &g2->num);
    result.den = poly_multiply(&g1->den, &g2->den);
    return result;
}

TransferFunction tf_parallel(const TransferFunction *g1, const TransferFunction *g2) {
    TransferFunction result;
    ComplexPoly num1_den2 = poly_multiply(&g1->num, &g2->den);
    ComplexPoly num2_den1 = poly_multiply(&g2->num, &g1->den);
    result.den = poly_multiply(&g1->den, &g2->den);
    /* Scale numerators by gains */
    for (int i = 0; i <= num1_den2.degree; i++)
        num1_den2.coeffs[i] = complex_scale(num1_den2.coeffs[i], g1->gain);
    for (int i = 0; i <= num2_den1.degree; i++)
        num2_den1.coeffs[i] = complex_scale(num2_den1.coeffs[i], g2->gain);
    result.num = poly_add(&num1_den2, &num2_den1);
    result.gain = 1.0;
    free(num1_den2.coeffs);
    free(num2_den1.coeffs);
    return result;
}

TransferFunction tf_feedback(const TransferFunction *g, const TransferFunction *h) {
    TransferFunction result;
    ComplexPoly num_g_h = poly_multiply(&g->num, &h->num);
    ComplexPoly den_g_h = poly_multiply(&g->den, &h->den);
    /* Characteristic denominator: 1 + G*H = (den_g*den_h + num_g*num_h) / (den_g*den_h) */
    /* Scale num_g_h by g->gain * h->gain */
    double gh_gain = g->gain * h->gain;
    for (int i = 0; i <= num_g_h.degree; i++)
        num_g_h.coeffs[i] = complex_scale(num_g_h.coeffs[i], gh_gain);
    result.den = poly_add(&den_g_h, &num_g_h);
    result.num = poly_multiply(&g->num, &h->den);
    for (int i = 0; i <= result.num.degree; i++)
        result.num.coeffs[i] = complex_scale(result.num.coeffs[i], g->gain);
    result.gain = 1.0;
    free(num_g_h.coeffs);
    free(den_g_h.coeffs);
    return result;
}


/* Frequency response: evaluate G(j*omega) */
FreqRespPoint tf_freq_resp(const TransferFunction *tf, double omega) {
    FreqRespPoint frp;
    frp.omega = omega;
    if (!tf) {
        frp.magnitude = 0; frp.magnitude_db = -INFINITY;
        frp.phase_rad = 0; frp.phase_deg = 0;
        return frp;
    }
    Complex jw = {0.0, omega};
    Complex G_jw = tf_eval(tf, jw);
    frp.magnitude = complex_abs(G_jw);
    frp.magnitude_db = 20.0 * log10(frp.magnitude > 1e-30 ? frp.magnitude : 1e-30);
    frp.phase_rad = atan2(G_jw.im, G_jw.re);
    frp.phase_deg = frp.phase_rad * 180.0 / M_PI;
    return frp;
}

void tf_bode_data(const TransferFunction *tf, const double *omegas,
                  int n_pts, FreqRespPoint *out) {
    if (!tf || !omegas || !out || n_pts <= 0) return;
    for (int i = 0; i < n_pts; i++) {
        out[i] = tf_freq_resp(tf, omegas[i]);
    }
}

double tf_bandwidth(const TransferFunction *tf, double omega_min,
                    double omega_max, int n_search) {
    if (!tf || n_search < 2) return 0.0;
    FreqRespPoint dc = tf_freq_resp(tf, 0.0);
    double target_mag = dc.magnitude / sqrt(2.0);
    double bw = omega_min;
    for (int i = 0; i < n_search; i++) {
        double w = omega_min + (omega_max - omega_min) * i / (n_search - 1);
        FreqRespPoint frp = tf_freq_resp(tf, w);
        if (frp.magnitude <= target_mag) { bw = w; break; }
    }
    return bw;
}

/* ========================================================================
 * L5: Nyquist Stability Criterion
 *
 * Nyquist contour: imaginary axis from -jR to +jR plus semicircle.
 * Evaluate L(j*omega) for omega in [omega_min, omega_max].
 * Count clockwise encirclements of -1 point.
 *
 * Theorem (Nyquist, 1932):
 *   Z = N + P
 *   Z = # unstable closed-loop poles
 *   N = # clockwise encirclements of -1
 *   P = # unstable open-loop poles
 * ======================================================================== */

Contour nyquist_contour(double R, int n_segments) {
    Contour c;
    c.t0 = 0.0;
    c.t1 = (double)n_segments;
    c.closed = 1;
    c.param = NULL;
    (void)R;
    return c;
}

void nyquist_evaluate(const TransferFunction *loop_tf,
                      const double *omegas, int n_pts,
                      Complex *nyquist_pts) {
    if (!loop_tf || !omegas || !nyquist_pts || n_pts <= 0) return;
    for (int i = 0; i < n_pts; i++) {
        Complex jw = {0.0, omegas[i]};
        nyquist_pts[i] = tf_eval(loop_tf, jw);
    }
}

int nyquist_encirclements(const Complex *nyquist_pts, int n_pts, Complex p) {
    if (!nyquist_pts || n_pts < 2) return 0;
    double total_angle = 0.0;
    for (int i = 0; i < n_pts - 1; i++) {
        Complex diff_a = complex_sub(nyquist_pts[i], p);
        Complex diff_b = complex_sub(nyquist_pts[i+1], p);
        double a1 = atan2(diff_a.im, diff_a.re);
        double a2 = atan2(diff_b.im, diff_b.re);
        double delta = a2 - a1;
        if (delta > M_PI) delta -= 2.0 * M_PI;
        else if (delta < -M_PI) delta += 2.0 * M_PI;
        total_angle += delta;
    }
    return (int)round(total_angle / (2.0 * M_PI));
}

int nyquist_stability_check(const TransferFunction *loop_tf,
                             double omega_min, double omega_max,
                             int n_pts, double R) {
    if (!loop_tf || n_pts < 2) return 0;
    /* Count open-loop unstable poles (Re > 0) */
    Complex *poles = (Complex *)calloc(loop_tf->den.degree, sizeof(Complex));
    int n_poles = poly_roots_newton(&loop_tf->den, poles,
                                     loop_tf->den.degree, 1e-8, 200);
    int P = 0;
    for (int i = 0; i < n_poles; i++) {
        if (poles[i].re > 0.0) P++;
    }
    free(poles);
    /* Evaluate Nyquist plot */
    Complex *nyq = (Complex *)malloc(n_pts * sizeof(Complex));
    if (!nyq) return 0;
    double *omegas = (double *)malloc(n_pts * sizeof(double));
    for (int i = 0; i < n_pts; i++)
        omegas[i] = omega_min + (omega_max - omega_min) * i / (n_pts - 1);
    nyquist_evaluate(loop_tf, omegas, n_pts, nyq);
    /* Encirclements of -1 */
    Complex minus_one = {-1.0, 0.0};
    int N = nyquist_encirclements(nyq, n_pts, minus_one);
    free(nyq);
    free(omegas);
    (void)R;
    /* Z = N + P; stable iff Z = 0 */
    return (N + P == 0) ? 1 : 0;
}


/* ========================================================================
 * L5: Gain and Phase Margins
 *
 * Gain margin (GM): 20*log10(1/|L(j*omega_pc)|) dB
 *   where omega_pc is phase crossover (phase = -180 deg)
 *
 * Phase margin (PM): 180 + arg(L(j*omega_gc)) degrees
 *   where omega_gc is gain crossover (|L| = 1 = 0 dB)
 *
 * Rule of thumb: GM >= 6 dB, PM >= 45 deg for robust stability.
 * ======================================================================== */

StabilityMargins tf_stability_margins(const TransferFunction *loop_tf,
                                       double omega_min, double omega_max,
                                       int n_search) {
    StabilityMargins sm;
    sm.gain_margin_db = INFINITY;
    sm.phase_margin_deg = INFINITY;
    sm.phase_crossover = 0.0;
    sm.gain_crossover = 0.0;
    if (!loop_tf || n_search < 2) return sm;

    double *omegas = (double *)malloc(n_search * sizeof(double));
    FreqRespPoint *frp = (FreqRespPoint *)malloc(n_search * sizeof(FreqRespPoint));
    for (int i = 0; i < n_search; i++)
        omegas[i] = omega_min + (omega_max - omega_min) * i / (n_search - 1);
    tf_bode_data(loop_tf, omegas, n_search, frp);

    /* Find phase crossover: phase closest to -180 deg */
    for (int i = 0; i < n_search; i++) {
        double phase_dist = fabs(frp[i].phase_deg + 180.0);
        if (phase_dist < 5.0 || (i > 0 &&
            (frp[i-1].phase_deg + 180.0) * (frp[i].phase_deg + 180.0) < 0)) {
            sm.phase_crossover = frp[i].omega;
            sm.gain_margin_db = -frp[i].magnitude_db;
            break;
        }
    }

    /* Find gain crossover: |L| closest to 1 (0 dB) */
    for (int i = 0; i < n_search; i++) {
        if (fabs(frp[i].magnitude_db) < 1.0 || (i > 0 &&
            frp[i-1].magnitude_db * frp[i].magnitude_db < 0)) {
            sm.gain_crossover = frp[i].omega;
            sm.phase_margin_deg = frp[i].phase_deg + 180.0;
            if (sm.phase_margin_deg > 180.0) sm.phase_margin_deg -= 360.0;
            break;
        }
    }

    free(omegas);
    free(frp);
    return sm;
}

/* ========================================================================
 * L5: Root Locus
 *
 * Characteristic equation: 1 + K*G(s)*H(s) = 0
 * As K varies 0 -> inf, closed-loop poles trace the root locus.
 *
 * Angle condition: arg(G(s)H(s)) = 180 + k*360 degrees
 * Magnitude condition: K = 1/|G(s)H(s)|
 *
 * This implementation traces branches starting from open-loop poles.
 * ======================================================================== */

int root_locus_branch(const TransferFunction *ol_tf, Complex start_pole,
                       double k_max, int n_steps, RootLocusPoint *out) {
    if (!ol_tf || !out || n_steps <= 0) return 0;
    Complex z = start_pole;
    double k_step = k_max / n_steps;
    for (int step = 0; step < n_steps; step++) {
        double K = (step + 1) * k_step;
        /* Solve 1 + K*G(s) = 0 => find s s.t. K*G(s) = -1 */
        /* Use Newton on f(s) = K*G(s) + 1 = 0 */
        for (int it = 0; it < 50; it++) {
            Complex G = tf_eval(ol_tf, z);
            Complex f_val;
            f_val.re = K * G.re + 1.0;
            f_val.im = K * G.im;
            if (complex_abs(f_val) < 1e-8) break;
            /* Approximate derivative */
            Complex zp = {z.re + 1e-6, z.im};
            Complex Gp = tf_eval(ol_tf, zp);
            Complex df_val;
            df_val.re = K * (Gp.re - G.re) / 1e-6;
            df_val.im = K * (Gp.im - G.im) / 1e-6;
            double d2 = complex_abs2(df_val);
            if (d2 < 1e-30) break;
            Complex step_c = complex_div(f_val, df_val);
            z = complex_sub(z, step_c);
        }
        out[step].pole = z;
        out[step].gain = K;
    }
    return n_steps;
}

int root_locus_full(const TransferFunction *ol_tf, double k_max,
                     int n_steps, RootLocusPoint *out, int max_pts) {
    if (!ol_tf || !out || max_pts <= 0) return 0;
    int deg = ol_tf->den.degree;
    if (deg <= 0) return 0;
    Complex *poles = (Complex *)calloc(deg, sizeof(Complex));
    int n_poles = poly_roots_newton(&ol_tf->den, poles, deg, 1e-8, 200);
    int total_pts = 0;
    for (int p = 0; p < n_poles && total_pts + n_steps <= max_pts; p++) {
        int n = root_locus_branch(ol_tf, poles[p], k_max, n_steps,
                                   out + total_pts);
        total_pts += n;
    }
    free(poles);
    return total_pts;
}


/* ========================================================================
 * L7: PID Controller Design
 *
 * PID: C(s) = Kp + Ki/s + Kd*s
 *
 * With derivative filter for properness: C(s) = Kp + Ki/s + Kd*s/(1 + s/N)
 *
 * PI design for desired phase margin:
 *   desired phase boost = target_PM - (180 + arg(G(j*wc)))
 *   Ki/wc = tan(phase_boost) * Kp
 *   |C(j*wc)*G(j*wc)| = 1 gives Kp
 * ======================================================================== */

TransferFunction pid_to_tf(const PIDController *pid) {
    TransferFunction tf;
    /* C(s) = Kp + Ki/s + Kd*s/(1 + s/N) */
    ComplexPoly num, den;
    /* Denominator: s*(1 + s/N) = s + s^2/N */
    den.degree = 2;
    den.coeffs = (Complex *)calloc(3, sizeof(Complex));
    den.coeffs[0] = (Complex){0,0};      /* constant */
    den.coeffs[1] = (Complex){1,0};      /* s */
    den.coeffs[2] = (Complex){1.0/pid->N, 0}; /* s^2/N */
    /* Numerator: Kp*s + Ki*(1+s/N) + Kd*s^2 */
    num.degree = 2;
    num.coeffs = (Complex *)calloc(3, sizeof(Complex));
    num.coeffs[0] = (Complex){pid->Ki, 0};                         /* Ki */
    num.coeffs[1] = (Complex){pid->Kp + pid->Ki/pid->N, 0};       /* Kp + Ki/N */
    num.coeffs[2] = (Complex){pid->Kp/pid->N + pid->Kd, 0};       /* Kp/N + Kd */
    tf.num = num;
    tf.den = den;
    tf.gain = 1.0;
    return tf;
}

/* PI controller design for target phase margin */
PIDController design_pi_phase_margin(const TransferFunction *plant,
                                      double target_pm_deg,
                                      double target_wc,
                                      int n_search) {
    PIDController pid = {0, 0, 0, 100.0};
    if (!plant) return pid;

    /* Evaluate plant at target crossover */
    FreqRespPoint plant_frp = tf_freq_resp(plant, target_wc);
    double plant_mag = plant_frp.magnitude;
    double plant_phase_deg = plant_frp.phase_deg;

    /* Required phase boost from PI controller */
    double pm_rad = target_pm_deg * M_PI / 180.0;
    double plant_phase_rad = plant_phase_deg * M_PI / 180.0;
    double phase_boost = pm_rad - M_PI - plant_phase_rad;
    while (phase_boost < -M_PI) phase_boost += 2.0*M_PI;
    while (phase_boost > M_PI) phase_boost -= 2.0*M_PI;
    /* PI phase boost: arg(1 + Ki/(j*wc*Kp)) = atan(-Ki/(wc*Kp)) */
    double tan_boost = tan(phase_boost);
    /* Ki/Kp = -wc * tan(phase_boost) */

    /* Magnitude condition: |C(j*wc)| * |G(j*wc)| = 1 */
    /* |C(j*wc)| = sqrt(Kp^2 + (Ki/wc)^2) */
    /* Kp = 1 / (|G| * sqrt(1 + (Ki/(wc*Kp))^2)) */
    double ratio = -target_wc * tan_boost; /* Ki/Kp */
    if (ratio < 0) ratio = 0.1; /* ensure positive Ki for stability */
    double denom_factor = sqrt(1.0 + ratio * ratio / (target_wc * target_wc));
    pid.Kp = 1.0 / (plant_mag * denom_factor);
    pid.Ki = ratio * pid.Kp;
    pid.Kd = 0.0;
    pid.N = 100.0;
    (void)n_search;
    return pid;
}

/* ========================================================================
 * L7: Additional control utilities
 * ======================================================================== */

/* Step response characteristics via inverse Laplace (simplified) */
double step_response_steady_state(const TransferFunction *cl_tf) {
    if (!cl_tf) return 0.0;
    /* Final value theorem: lim_{t->inf} y(t) = lim_{s->0} s * (1/s) * G_cl(s) */
    /* = G_cl(0) = tf_eval(cl_tf, 0) */
    Complex dc = tf_eval(cl_tf, (Complex){0,0});
    return dc.re;
}

/* DC gain of a transfer function */
double tf_dc_gain(const TransferFunction *tf) {
    if (!tf) return 0.0;
    Complex s_zero = {0.0, 0.0};
    Complex val = tf_eval(tf, s_zero);
    return val.re;
}

/* Check if a system is minimum-phase (all zeros in left half-plane) */
int tf_is_minimum_phase(const TransferFunction *tf, double tol) {
    if (!tf) return 0;
    int deg = tf->num.degree;
    if (deg <= 0) return 1;
    Complex *zeros = (Complex *)calloc(deg, sizeof(Complex));
    int n_zeros = poly_roots_newton(&tf->num, zeros, deg, 1e-8, 200);
    int min_phase = 1;
    for (int i = 0; i < n_zeros; i++) {
        if (zeros[i].re > tol) { min_phase = 0; break; }
    }
    free(zeros);
    return min_phase;
}

/* ========================================================================
 * L7: Typical system models for application examples
 * ======================================================================== */

/* First-order system: K / (tau*s + 1) */
TransferFunction tf_first_order(double K, double tau) {
    TransferFunction tf;
    tf.gain = K;
    tf.num.degree = 0;
    tf.num.coeffs = (Complex *)calloc(1, sizeof(Complex));
    tf.num.coeffs[0] = (Complex){1.0, 0};
    tf.den.degree = 1;
    tf.den.coeffs = (Complex *)calloc(2, sizeof(Complex));
    tf.den.coeffs[0] = (Complex){1.0, 0};
    tf.den.coeffs[1] = (Complex){tau, 0};
    return tf;
}

/* Second-order system: omega_n^2 / (s^2 + 2*zeta*omega_n*s + omega_n^2) */
TransferFunction tf_second_order(double omega_n, double zeta) {
    TransferFunction tf;
    tf.gain = 1.0;
    tf.num.degree = 0;
    tf.num.coeffs = (Complex *)calloc(1, sizeof(Complex));
    tf.num.coeffs[0] = (Complex){omega_n * omega_n, 0};
    tf.den.degree = 2;
    tf.den.coeffs = (Complex *)calloc(3, sizeof(Complex));
    tf.den.coeffs[0] = (Complex){omega_n * omega_n, 0};
    tf.den.coeffs[1] = (Complex){2.0 * zeta * omega_n, 0};
    tf.den.coeffs[2] = (Complex){1.0, 0};
    return tf;
}

/* Integrator: 1/s */
TransferFunction tf_integrator(void) {
    TransferFunction tf;
    tf.gain = 1.0;
    tf.num.degree = 0;
    tf.num.coeffs = (Complex *)calloc(1, sizeof(Complex));
    tf.num.coeffs[0] = (Complex){1.0, 0};
    tf.den.degree = 1;
    tf.den.coeffs = (Complex *)calloc(2, sizeof(Complex));
    tf.den.coeffs[0] = (Complex){0.0, 0};
    tf.den.coeffs[1] = (Complex){1.0, 0};
    return tf;
}

/* Time delay approximation (1st-order Pade): e^{-Td*s} ~ (1 - Td*s/2)/(1 + Td*s/2) */
TransferFunction tf_pade_delay(double Td) {
    TransferFunction tf;
    tf.gain = 1.0;
    tf.num.degree = 1;
    tf.num.coeffs = (Complex *)calloc(2, sizeof(Complex));
    tf.num.coeffs[0] = (Complex){1.0, 0};
    tf.num.coeffs[1] = (Complex){-Td / 2.0, 0};
    tf.den.degree = 1;
    tf.den.coeffs = (Complex *)calloc(2, sizeof(Complex));
    tf.den.coeffs[0] = (Complex){1.0, 0};
    tf.den.coeffs[1] = (Complex){Td / 2.0, 0};
    return tf;
}


/* ========================================================================
 * L7: Lead-Lag Compensator Design
 *
 * Lead: C(s) = K * (T*s + 1)/(alpha*T*s + 1), alpha < 1 (phase advance)
 * Lag:  C(s) = K * (T*s + 1)/(beta*T*s + 1), beta > 1 (gain reduction)
 *
 * Lead compensator provides phase boost at target frequency.
 * Lag compensator provides low-frequency gain increase.
 * ======================================================================== */

TransferFunction tf_lead_compensator(double K, double T, double alpha) {
    TransferFunction tf;
    tf.gain = K;
    tf.num.degree = 1;
    tf.num.coeffs = (Complex *)calloc(2, sizeof(Complex));
    tf.num.coeffs[0] = (Complex){1.0, 0};
    tf.num.coeffs[1] = (Complex){T, 0};
    tf.den.degree = 1;
    tf.den.coeffs = (Complex *)calloc(2, sizeof(Complex));
    tf.den.coeffs[0] = (Complex){1.0, 0};
    tf.den.coeffs[1] = (Complex){alpha * T, 0};
    return tf;
}

TransferFunction tf_lag_compensator(double K, double T, double beta) {
    TransferFunction tf;
    tf.gain = K;
    tf.num.degree = 1;
    tf.num.coeffs = (Complex *)calloc(2, sizeof(Complex));
    tf.num.coeffs[0] = (Complex){1.0, 0};
    tf.num.coeffs[1] = (Complex){T, 0};
    tf.den.degree = 1;
    tf.den.coeffs = (Complex *)calloc(2, sizeof(Complex));
    tf.den.coeffs[0] = (Complex){1.0, 0};
    tf.den.coeffs[1] = (Complex){beta * T, 0};
    return tf;
}

/* Lead compensator design for desired phase margin */
PIDController design_lead_phase_margin(const TransferFunction *plant,
                                        double target_pm_deg,
                                        double target_wc) {
    PIDController pid = {0, 0, 0, 100.0};
    if (!plant) return pid;
    FreqRespPoint frp = tf_freq_resp(plant, target_wc);
    double needed_boost = target_pm_deg - (frp.phase_deg + 180.0);
    while (needed_boost < -360) needed_boost += 360;
    while (needed_boost > 360) needed_boost -= 360;
    double alpha = (1.0 - sin(needed_boost * M_PI / 180.0))
                 / (1.0 + sin(needed_boost * M_PI / 180.0));
    if (alpha < 0.01) alpha = 0.01;
    if (alpha > 1.0) alpha = 1.0;
    double T = 1.0 / (target_wc * sqrt(alpha));
    double K = 1.0 / (frp.magnitude * sqrt(alpha));
    pid.Kp = K;
    pid.Ki = 0.0;
    pid.Kd = K * T;
    pid.N = 1.0 / (alpha * T);
    return pid;
}

