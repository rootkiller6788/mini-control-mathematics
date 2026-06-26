/**
 * @file complex_control.h
 * @brief Complex analysis for control theory
 *
 * L5-L7: Engineering methods and applications in automatic control
 */

#ifndef COMPLEX_CONTROL_H
#define COMPLEX_CONTROL_H

#include "complex_number.h"
#include "complex_integration.h"
#include <stddef.h>

typedef struct {
    ComplexPoly num;
    ComplexPoly den;
    double      gain;
} TransferFunction;

Complex tf_eval(const TransferFunction *tf, Complex s);
int tf_find_poles(const TransferFunction *tf, Complex *poles,
                  int max_poles, double tol, int max_iter);
int tf_find_zeros(const TransferFunction *tf, Complex *zeros,
                  int max_zeros, double tol, int max_iter);
int tf_is_stable(const TransferFunction *tf, double stability_margin);
int tf_relative_degree(const TransferFunction *tf);

TransferFunction tf_series(const TransferFunction *g1,
                            const TransferFunction *g2);
TransferFunction tf_parallel(const TransferFunction *g1,
                              const TransferFunction *g2);
TransferFunction tf_feedback(const TransferFunction *g,
                              const TransferFunction *h);

typedef struct {
    double omega;
    double magnitude;
    double magnitude_db;
    double phase_deg;
    double phase_rad;
} FreqRespPoint;

FreqRespPoint tf_freq_resp(const TransferFunction *tf, double omega);
void tf_bode_data(const TransferFunction *tf, const double *omegas,
                  int n_pts, FreqRespPoint *out);
double tf_bandwidth(const TransferFunction *tf, double omega_min,
                    double omega_max, int n_search);

Contour nyquist_contour(double R, int n_segments);
void nyquist_evaluate(const TransferFunction *loop_tf,
                      const double *omegas, int n_pts,
                      Complex *nyquist_pts);
int nyquist_encirclements(const Complex *nyquist_pts, int n_pts, Complex p);
int nyquist_stability_check(const TransferFunction *loop_tf,
                             double omega_min, double omega_max,
                             int n_pts, double R);

typedef struct {
    double gain_margin_db;
    double phase_margin_deg;
    double phase_crossover;
    double gain_crossover;
} StabilityMargins;

StabilityMargins tf_stability_margins(const TransferFunction *loop_tf,
                                       double omega_min, double omega_max,
                                       int n_search);

typedef struct {
    Complex pole;
    double  gain;
} RootLocusPoint;

int root_locus_branch(const TransferFunction *ol_tf, Complex start_pole,
                       double k_max, int n_steps, RootLocusPoint *out);
int root_locus_full(const TransferFunction *ol_tf, double k_max,
                     int n_steps, RootLocusPoint *out, int max_pts);

typedef struct {
    double Kp;
    double Ki;
    double Kd;
    double N;
} PIDController;

TransferFunction pid_to_tf(const PIDController *pid);
PIDController design_pi_phase_margin(const TransferFunction *plant,
                                      double target_pm_deg,
                                      double target_wc,
                                      int n_search);


/* Additional transfer function constructors */
TransferFunction tf_first_order(double K, double tau);
TransferFunction tf_second_order(double omega_n, double zeta);
TransferFunction tf_integrator(void);
TransferFunction tf_pade_delay(double Td);
TransferFunction tf_lead_compensator(double K, double T, double alpha);
TransferFunction tf_lag_compensator(double K, double T, double beta);

/* Additional design and analysis */
double step_response_steady_state(const TransferFunction *cl_tf);
double tf_dc_gain(const TransferFunction *tf);
int tf_is_minimum_phase(const TransferFunction *tf, double tol);
PIDController design_lead_phase_margin(const TransferFunction *plant, double target_pm_deg, double target_wc);

#endif /* COMPLEX_CONTROL_H */
