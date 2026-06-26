/**
 * @file    ode_applications.c
 * @brief   Real-world ODE applications: DC motor, RLC circuit,
 *          spring-mass-damper, population models, Lorenz system,
 *          damped pendulum. Each implements an independent L7 application.
 *
 * L7: Applications – engineering case studies with real parameters.
 *
 * Knowledge points:
 *   ode_rhs_dc_motor                  → L7: DC motor state-space model
 *   ode_dc_motor_simulate             → L7: DC motor step response (servo control)
 *   ode_dc_motor_transfer_function    → L7: Motor TF for control design
 *   ode_rhs_rlc_circuit               → L7: RLC circuit model
 *   ode_rlc_simulate                  → L7: RLC transient analysis
 *   ode_rlc_damping_analysis          → L7: RLC damping classification
 *   ode_spring_mass_simulate          → L7: Spring-mass-damper
 *   ode_oscillator_frequency_response → L7: Frequency response (Bode precursor)
 *   ode_find_resonance                → L7: Resonant frequency
 *   ode_predator_prey_simulate        → L7: Population dynamics
 *   ode_lotka_volterra_equilibrium    → L7: LV equilibrium
 *   ode_lotka_volterra_integral       → L7: LV conserved quantity
 *   ode_lorenz_simulate               → L7: Lorenz attractor
 *   ode_lorenz_equilibria             → L7: Lorenz equilibrium analysis
 *   ode_pendulum_simulate             → L7: Pendulum dynamics
 *   ode_pendulum_period               → L7: Pendulum period (small/large angle)
 */

#include "ode_applications.h"
#include "ode_nonlinear.h"
#include "ode_numerical.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ──────────────────────── Helpers ──────────────────────────────────── */

static void vec_copy_n(const double *s, double *d, int n) {
    memcpy(d, s, (size_t)n * sizeof(double));
}

/* ──────────────────────── L7: DC Motor ──────────────────────────────── */

int ode_rhs_dc_motor(double t, const double *y, double *f_out,
                       int dim, void *ctx) {
    (void)t;
    if (!y || !f_out || !ctx || dim != 3) return -1;
    DCMotorParams *p = (DCMotorParams *)ctx;
    /* States: x₁=θ, x₂=ω, x₃=i */
    f_out[0] = y[1];                                  /* θ̇ = ω */
    f_out[1] = -(p->b / p->J) * y[1] + (p->K / p->J) * y[2];  /* ω̇ */
    f_out[2] = -(p->K / p->L) * y[1] - (p->R / p->L) * y[2] + p->V / p->L; /* i̇ */
    return 0;
}

int ode_dc_motor_simulate(const DCMotorParams *params, double t_end,
                            int n_steps, double *theta, double *omega,
                            double *current, double *t_vals) {
    if (!params || !theta || !omega || !current || !t_vals || n_steps < 2)
        return -1;

    double h = t_end / (n_steps - 1);
    double y[3] = {0.0, 0.0, 0.0}; /* Start from rest */

    for (int i = 0; i < n_steps; i++) {
        t_vals[i] = i * h;
        theta[i] = y[0];
        omega[i] = y[1];
        current[i] = y[2];

        /* Forward Euler */
        double f[3];
        ode_rhs_dc_motor(t_vals[i], y, f, 3, (void *)params);
        for (int j = 0; j < 3; j++) y[j] += h * f[j];
    }

    return 0;
}

void ode_dc_motor_transfer_function(const DCMotorParams *params,
                                      double *num, double *den) {
    if (!params || !num || !den) return;

    /* G(s) = ω(s)/V(s) = K / ((Js+b)(Ls+R) + K²) */
    /* = K / (JLs² + (JR+bL)s + (bR + K²)) */
    double a2 = params->J * params->L;
    double a1 = params->J * params->R + params->b * params->L;
    double a0 = params->b * params->R + params->K * params->K;

    num[0] = params->K;
    den[0] = a2;
    den[1] = a1;
    den[2] = a0;
}

/* ──────────────────────── L7: RLC Circuit ──────────────────────────── */

int ode_rhs_rlc_circuit(double t, const double *y, double *f_out,
                          int dim, void *ctx) {
    if (!y || !f_out || !ctx || dim != 2) return -1;
    RLCCircuitParams *p = (RLCCircuitParams *)ctx;

    /* States: x₁ = v_c, x₂ = i_L */
    /* C·v̇_c = i_L  →  v̇_c = i_L / C */
    /* L·i̇_L = v(t) - R·i_L - v_c  →  i̇_L = (v - R·i_L - v_c) / L */
    double v_source = p->v0 * cos(p->omega * t);

    f_out[0] = y[1] / p->C;                          /* v̇_c */
    f_out[1] = (v_source - p->R * y[1] - y[0]) / p->L;  /* i̇_L */
    return 0;
}

int ode_rlc_simulate(const RLCCircuitParams *params, double t_end,
                       int n_steps, double *v_c, double *i_l,
                       double *t_vals) {
    if (!params || !v_c || !i_l || !t_vals || n_steps < 2) return -1;

    double h = t_end / (n_steps - 1);
    double y[2] = {0.0, 0.0};

    for (int i = 0; i < n_steps; i++) {
        t_vals[i] = i * h;
        v_c[i] = y[0];
        i_l[i] = y[1];

        double f[2];
        ode_rhs_rlc_circuit(t_vals[i], y, f, 2, (void *)params);
        for (int j = 0; j < 2; j++) y[j] += h * f[j];
    }
    return 0;
}

int ode_rlc_damping_analysis(const RLCCircuitParams *params,
                               double *omega0, double *zeta) {
    if (!params || !omega0 || !zeta) return -1;

    *omega0 = 1.0 / sqrt(params->L * params->C);
    *zeta = (params->R / 2.0) * sqrt(params->C / params->L);

    /* Classification */
    if (*zeta < 1e-10) return 0;  /* Undamped */
    if (*zeta < 1.0) return 1;    /* Underdamped */
    if (fabs(*zeta - 1.0) < 1e-10) return 2;  /* Critically damped */
    return 3;  /* Overdamped */
}

/* ──────────────────────── L7: Spring-Mass-Damper ────────────────────── */

int ode_spring_mass_simulate(const OscillatorParams *params, double t_end,
                               int n_steps, double *x, double *v,
                               double *t_vals) {
    if (!params || !x || !v || !t_vals || n_steps < 2) return -1;

    double h = t_end / (n_steps - 1);
    double y[2] = {0.01, 0.0}; /* Small initial displacement */

    for (int i = 0; i < n_steps; i++) {
        t_vals[i] = i * h;
        x[i] = y[0];
        v[i] = y[1];

        double f[2];
        ode_rhs_harmonic_oscillator(t_vals[i], y, f, 2, (void *)params);
        for (int j = 0; j < 2; j++) y[j] += h * f[j];
    }
    return 0;
}

void ode_oscillator_frequency_response(const OscillatorParams *params,
                                         const double *omega, int n_freqs,
                                         double *M, double *phase) {
    if (!params || !omega || !M || !phase) return;

    double m = params->m, c = params->c, k = params->k;
    double omega_n = sqrt(k / m);
    double zeta = c / (2.0 * sqrt(m * k));

    for (int i = 0; i < n_freqs; i++) {
        double r = omega[i] / omega_n;
        double denom = sqrt((1.0 - r*r)*(1.0 - r*r) + 4.0*zeta*zeta*r*r);
        M[i] = (denom > 1e-14) ? 1.0 / denom : 0.0;
        phase[i] = atan2(-2.0 * zeta * r, 1.0 - r*r);
    }
}

bool ode_find_resonance(const OscillatorParams *params,
                          double *omega_res, double *M_max) {
    if (!params || !omega_res || !M_max) return false;

    double m = params->m, c = params->c, k = params->k;
    double omega_n = sqrt(k / m);
    double zeta = c / (2.0 * sqrt(m * k));

    if (zeta >= 1.0 / sqrt(2.0)) {
        *omega_res = 0.0;
        *M_max = 1.0;
        return false;
    }

    *omega_res = omega_n * sqrt(1.0 - 2.0 * zeta * zeta);
    double r_res = *omega_res / omega_n;
    double denom = sqrt((1.0 - r_res*r_res)*(1.0 - r_res*r_res) + 4.0*zeta*zeta*r_res*r_res);
    *M_max = 1.0 / denom;
    return true;
}

/* ──────────────────────── L7: Population Dynamics ───────────────────── */

int ode_predator_prey_simulate(const LotkaVolterraParams *params,
                                 double x0, double y0, double t_end,
                                 int n_steps, double *prey, double *pred,
                                 double *t_vals) {
    if (!params || !prey || !pred || !t_vals || n_steps < 2) return -1;

    double h = t_end / (n_steps - 1);
    double y[2] = {x0, y0};

    for (int i = 0; i < n_steps; i++) {
        t_vals[i] = i * h;
        prey[i] = y[0];
        pred[i] = y[1];

        double f[2];
        ode_rhs_lotka_volterra(t_vals[i], y, f, 2, (void *)params);
        for (int j = 0; j < 2; j++) y[j] += h * f[j];
    }
    return 0;
}

void ode_lotka_volterra_equilibrium(const LotkaVolterraParams *params,
                                      double *x_eq, double *y_eq) {
    if (!params || !x_eq || !y_eq) return;
    *x_eq = params->gamma / params->delta;
    *y_eq = params->alpha / params->beta;
}

double ode_lotka_volterra_integral(const LotkaVolterraParams *params,
                                     double x, double y) {
    if (!params || x <= 0.0 || y <= 0.0) return 0.0;
    return params->delta * x - params->gamma * log(x)
         + params->beta * y - params->alpha * log(y);
}

/* ──────────────────────── L7: Lorenz System ─────────────────────────── */

int ode_lorenz_simulate(const LorenzParams *params,
                          double x0, double y0, double z0,
                          double t_end, int n_steps,
                          double *x, double *y, double *z, double *t_vals) {
    if (!params || !x || !y || !z || !t_vals || n_steps < 2) return -1;

    double h = t_end / (n_steps - 1);
    double state[3] = {x0, y0, z0};

    for (int i = 0; i < n_steps; i++) {
        t_vals[i] = i * h;
        x[i] = state[0];
        y[i] = state[1];
        z[i] = state[2];

        double f[3];
        ode_rhs_lorenz(t_vals[i], state, f, 3, (void *)params);
        for (int j = 0; j < 3; j++) state[j] += h * f[j];
    }
    return 0;
}

void ode_lorenz_equilibria(const LorenzParams *params,
                             double *eq_points, int *n_eq) {
    if (!params || !eq_points || !n_eq) return;

    /* C₀ = (0, 0, 0) */
    eq_points[0] = 0.0; eq_points[1] = 0.0; eq_points[2] = 0.0;

    if (params->rho > 1.0) {
        double beta = sqrt(params->beta * (params->rho - 1.0));
        eq_points[3] = beta;
        eq_points[4] = beta;
        eq_points[5] = params->rho - 1.0;
        eq_points[6] = -beta;
        eq_points[7] = -beta;
        eq_points[8] = params->rho - 1.0;
        *n_eq = 3;
    } else {
        *n_eq = 1;
    }
}

/* ──────────────────────── L7: Pendulum ──────────────────────────────── */

int ode_pendulum_simulate(const PendulumParams *params,
                            double theta0, double omega0,
                            double t_end, int n_steps,
                            double *theta, double *omega, double *t_vals) {
    if (!params || !theta || !omega || !t_vals || n_steps < 2) return -1;

    double h = t_end / (n_steps - 1);
    double y[2] = {theta0, omega0};

    for (int i = 0; i < n_steps; i++) {
        t_vals[i] = i * h;
        theta[i] = y[0];
        omega[i] = y[1];

        double f[2];
        ode_rhs_pendulum(t_vals[i], y, f, 2, (void *)params);
        for (int j = 0; j < 2; j++) y[j] += h * f[j];
    }
    return 0;
}

void ode_pendulum_period(const PendulumParams *params, double theta0,
                           double *T_small, double *T_exact) {
    if (!params || !T_small || !T_exact) return;

    double omega0_sq = params->g / params->L;
    *T_small = 2.0 * M_PI / sqrt(omega0_sq);

    /* Large-angle correction via arithmetic-geometric mean */
    double k = sin(fabs(theta0) / 2.0);
    /* AGM-based complete elliptic integral K(k) approximation */
    /* Using polynomial approximation (Abramowitz & Stegun) */
    double k2 = k * k;
    double K_approx = M_PI / 2.0 * (1.0 + k2/4.0 + 9.0*k2*k2/64.0 + 25.0*k2*k2*k2/256.0);
    *T_exact = *T_small * (2.0 / M_PI) * K_approx;
}
