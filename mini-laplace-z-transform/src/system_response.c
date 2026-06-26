/**
 * @file system_response.c
 * @brief Time-Domain System Response — Step, Impulse, Ramp Analysis
 *
 * L5: Response computation via numerical inverse Laplace
 * L6: Transient spec extraction, 2nd-order closed-form solutions,
 *     steady-state error computation.
 */

#include "system_response.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*----------------------------------------------------------------------
 * L5 — Compute Time Response (Numerical Inverse Laplace)
 *----------------------------------------------------------------------*/

int compute_time_response(const TransferFunction *tf,
                           ResponseType resp_type,
                           double param,
                           double T_end, int n_samples,
                           TimeResponse *response)
{
    if (!tf || !response || n_samples < 2 || T_end <= 0.0) {
        errno = EINVAL;
        return -1;
    }
    memset(response, 0, sizeof(*response));

    /* Build Y(s) = G(s)·U(s) depending on input type */
    LaplaceRational Y;
    memset(&Y, 0, sizeof(Y));
    Y.numerator.order   = tf->num_order;
    Y.denominator.order = tf->den_order;
    for (int i = 0; i <= tf->num_order; i++) Y.numerator.coeff[i]   = tf->num_coeff[i];
    for (int i = 0; i <= tf->den_order; i++) Y.denominator.coeff[i] = tf->den_coeff[i];

    /* Multiply by input transform U(s) */
    /* Step: U(s) = 1/s → multiply denominator by s */
    /* Impulse: U(s) = 1 → no change */
    /* Ramp: U(s) = 1/s² → multiply denominator by s² */
    /* Sinusoid: U(s) = ω/(s²+ω²) */

    LaplacePolynomial U_num = {0, {1.0}};
    LaplacePolynomial U_den = {0, {1.0}};

    switch (resp_type) {
    case RESPONSE_STEP:
        U_den.order = 1;
        U_den.coeff[0] = 0.0;
        U_den.coeff[1] = 1.0;  /* D_U(s) = s */
        break;
    case RESPONSE_IMPULSE:
        /* U(s) = 1, no change */
        break;
    case RESPONSE_RAMP:
        U_den.order = 2;
        U_den.coeff[0] = 0.0;
        U_den.coeff[1] = 0.0;
        U_den.coeff[2] = 1.0;  /* D_U(s) = s² */
        break;
    case RESPONSE_SINUSOID:
        U_num.order = 0;
        U_num.coeff[0] = param;  /* ω */
        U_den.order = 2;
        U_den.coeff[0] = param * param;
        U_den.coeff[1] = 0.0;
        U_den.coeff[2] = 1.0;  /* D_U(s) = s² + ω² */
        break;
    default:
        errno = EINVAL;
        return -1;
    }

    /* Y(s) = G(s)·U(s) = N_G/N_U / (D_G·D_U) */
    LaplacePolynomial Y_num, Y_den;
    laplace_poly_multiply(&Y.numerator,   &U_num, &Y_num);
    laplace_poly_multiply(&Y.denominator, &U_den, &Y_den);

    LaplaceRational Y_rational;
    Y_rational.numerator   = Y_num;
    Y_rational.denominator = Y_den;

    /* Allocate response arrays */
    response->n_samples = n_samples;
    response->t_start   = 0.0;
    response->t_step    = T_end / (n_samples - 1);
    response->time      = (double*)malloc(n_samples * sizeof(double));
    response->output    = (double*)malloc(n_samples * sizeof(double));
    response->input     = (double*)malloc(n_samples * sizeof(double));
    if (!response->time || !response->output || !response->input) {
        time_response_free(response);
        return -1;
    }

    for (int i = 0; i < n_samples; i++) {
        double t = i * response->t_step;
        response->time[i] = t;

        /* Input signal at time t */
        switch (resp_type) {
        case RESPONSE_STEP:    response->input[i] = (t >= 0) ? 1.0 : 0.0; break;
        case RESPONSE_IMPULSE: response->input[i] = (i == 0) ? 1e10 : 0.0; break;
        case RESPONSE_RAMP:    response->input[i] = t; break;
        case RESPONSE_SINUSOID: response->input[i] = sin(param * t); break;
        default: response->input[i] = 0.0; break;
        }

        /* Output: numerical inverse Laplace */
        response->output[i] = laplace_inverse_numerical(&Y_rational, t, 100);
    }

    return 0;
}

/*----------------------------------------------------------------------
 * L5 — Extract Transient Specifications from Step Response
 *----------------------------------------------------------------------*/

int extract_transient_specs(const TimeResponse *response,
                             TransientSpecs *specs)
{
    if (!response || !specs || response->n_samples < 2) {
        errno = EINVAL;
        return -1;
    }

    memset(specs, 0, sizeof(*specs));

    int N = response->n_samples;
    const double *t = response->time;
    const double *y = response->output;

    /* Steady-state value (last sample) */
    double y_ss = y[N - 1];
    specs->steady_state_value = y_ss;
    specs->steady_state_error = 1.0 - y_ss;  /* For unity step reference */

    /* Find peak (overshoot) */
    double y_max = y[0];
    int peak_idx = 0;
    for (int i = 1; i < N; i++) {
        if (y[i] > y_max) {
            y_max = y[i];
            peak_idx = i;
        }
    }
    specs->peak_time = t[peak_idx];
    if (y_ss > 0.0 && y_max > y_ss) {
        specs->overshoot_percent = 100.0 * (y_max - y_ss) / y_ss;
    } else {
        specs->overshoot_percent = 0.0;
    }

    /* Rise time (10% → 90%) */
    double y10 = 0.1 * y_ss;
    double y90 = 0.9 * y_ss;
    double t10 = 0.0, t90 = 0.0;

    for (int i = 1; i < N; i++) {
        if (y[i-1] <= y10 && y[i] >= y10)
            t10 = t[i-1] + (t[i] - t[i-1]) * (y10 - y[i-1]) / (y[i] - y[i-1]);
        if (y[i-1] <= y90 && y[i] >= y90)
            t90 = t[i-1] + (t[i] - t[i-1]) * (y90 - y[i-1]) / (y[i] - y[i-1]);
    }
    specs->rise_time = t90 - t10;
    if (specs->rise_time < 0.0) specs->rise_time = 0.0;

    /* Delay time (50% rise) */
    double y50 = 0.5 * y_ss;
    for (int i = 1; i < N; i++) {
        if (y[i-1] <= y50 && y[i] >= y50) {
            specs->delay_time = t[i-1] + (t[i] - t[i-1]) * (y50 - y[i-1]) / (y[i] - y[i-1]);
            break;
        }
    }

    /* Settling time (last time response enters and stays within ±2% of y_ss) */
    double band = 0.02 * fabs(y_ss);
    if (band < 1e-9) band = 0.02;
    double t_settle = 0.0;
    for (int i = N - 1; i >= 0; i--) {
        if (fabs(y[i] - y_ss) > band) {
            if (i < N - 1) t_settle = t[i + 1];
            else t_settle = t[N - 1];
            break;
        }
    }
    specs->settling_time = t_settle;

    specs->is_valid = 1;
    return 0;
}

/*----------------------------------------------------------------------
 * L5 — Closed-Form 2nd-Order Specifications
 *----------------------------------------------------------------------*/

int second_order_specs_closed_form(double omega_n, double zeta,
                                    TransientSpecs *specs)
{
    if (!specs || omega_n <= 0.0 || zeta < 0.0) { errno = EINVAL; return -1; }
    memset(specs, 0, sizeof(*specs));

    specs->steady_state_value = 1.0;
    specs->steady_state_error = 0.0;

    if (zeta == 0.0) {
        /* Undamped */
        specs->rise_time  = 1.8 / omega_n;
        specs->peak_time  = M_PI / (2.0 * omega_n);
        specs->overshoot_percent = 100.0;
        specs->settling_time = 1e300;  /* Never settles */
        specs->delay_time = 0.7 / omega_n;
    } else if (zeta < 1.0) {
        /* Underdamped */
        double wd = omega_n * sqrt(1.0 - zeta * zeta);
        double sigma = zeta * omega_n;
        double phi = acos(zeta);

        specs->rise_time  = (M_PI - phi) / wd;  /* 0-100% rise */
        specs->peak_time  = M_PI / wd;
        specs->overshoot_percent = 100.0 * exp(-M_PI * zeta / sqrt(1.0 - zeta * zeta));
        specs->settling_time = 4.0 / sigma;  /* 2% criterion */
        specs->delay_time = (1.0 + 0.7 * zeta) / omega_n;  /* approximation */
    } else if (zeta == 1.0) {
        /* Critically damped */
        specs->rise_time  = 3.358 / omega_n;
        specs->peak_time  = 1e300;  /* No overshoot */
        specs->overshoot_percent = 0.0;
        specs->settling_time = 4.75 / omega_n;
        specs->delay_time = 1.1 / omega_n;
    } else {
        /* Overdamped */
        double s2 = omega_n * (-zeta + sqrt(zeta * zeta - 1.0));
        specs->rise_time  = (1.0 - 0.4167 * zeta + 2.917 * zeta * zeta) / omega_n;
        specs->peak_time  = 1e300;
        specs->overshoot_percent = 0.0;
        specs->settling_time = 4.0 / (-s2);  /* Dominant slow pole */
        specs->delay_time = (1.0 + 0.6 * zeta + 0.15 * zeta * zeta) / omega_n;
    }

    specs->is_valid = 1;
    return 0;
}

/*----------------------------------------------------------------------
 * L6 — Second-Order Step Response at Time t (Analytical)
 *----------------------------------------------------------------------*/

double second_order_step_response_at(double omega_n, double zeta, double t)
{
    if (omega_n <= 0.0 || zeta < 0.0 || t < 0.0) {
        errno = EINVAL;
        return NAN;
    }

    if (t == 0.0) return 0.0;

    double wn = omega_n;

    if (zeta == 0.0) {
        /* Undamped: y(t) = 1 - cos(ωₙt) */
        return 1.0 - cos(wn * t);
    } else if (zeta < 1.0) {
        /* Underdamped */
        double wd = wn * sqrt(1.0 - zeta * zeta);
        double sigma = zeta * wn;
        double phi = acos(zeta);
        return 1.0 - (exp(-sigma * t) / sqrt(1.0 - zeta * zeta)) * sin(wd * t + phi);
    } else if (zeta == 1.0) {
        /* Critically damped: y(t) = 1 - e^{-ωₙt}·(1 + ωₙt) */
        return 1.0 - exp(-wn * t) * (1.0 + wn * t);
    } else {
        /* Overdamped: y(t) = 1 + (1/(s₂-s₁))·(s₂·e^{s₁t} - s₁·e^{s₂t}) */
        double disc = sqrt(zeta * zeta - 1.0);
        double s1 = wn * (-zeta - disc);
        double s2 = wn * (-zeta + disc);
        return 1.0 + (s2 * exp(s1 * t) - s1 * exp(s2 * t)) / (s2 - s1);
    }
}

/*----------------------------------------------------------------------
 * L6 — Impulse Response via Analytic Residues
 *----------------------------------------------------------------------*/

int impulse_response_analytic(const TransferFunction *tf,
                               double T_end, int n_samples,
                               TimeResponse *response)
{
    if (!tf || !response || n_samples < 2 || T_end <= 0.0) {
        errno = EINVAL;
        return -1;
    }

    /* For strictly proper rational G(s), h(t) = Σ A_k·e^{p_k·t} */
    LaplaceRational L;
    memset(&L, 0, sizeof(L));
    L.numerator.order   = tf->num_order;
    L.denominator.order = tf->den_order;
    for (int i = 0; i <= tf->num_order; i++) L.numerator.coeff[i]   = tf->num_coeff[i];
    for (int i = 0; i <= tf->den_order; i++) L.denominator.coeff[i] = tf->den_coeff[i];

    LaplacePoleZero pz;
    memset(&pz, 0, sizeof(pz));
    if (laplace_find_poles(&L, &pz) != 0) return -1;

    /* Compute residues */
    int n_poles = pz.num_poles;
    double residues[LAPLACE_MAX_ORDER];
    double complex p[LAPLACE_MAX_ORDER];

    /* For simple poles: residue A_k = N(p_k)/D'(p_k) */
    LaplacePolynomial dD;
    laplace_poly_derivative(&L.denominator, &dD);

    for (int k = 0; k < n_poles; k++) {
        p[k] = pz.poles[k];
        double complex N_val = laplace_poly_eval(&L.numerator, p[k]);
        double complex Dp_val = laplace_poly_eval(&dD, p[k]);
        if (cabs(Dp_val) > 1e-12)
            residues[k] = creal(N_val / Dp_val) * tf->gain;
        else
            residues[k] = 0.0;
    }

    /* Build response */
    response->n_samples = n_samples;
    response->t_start   = 0.0;
    response->t_step    = T_end / (n_samples - 1);
    response->time      = (double*)malloc(n_samples * sizeof(double));
    response->output    = (double*)malloc(n_samples * sizeof(double));
    response->input     = (double*)malloc(n_samples * sizeof(double));

    if (!response->time || !response->output || !response->input) {
        time_response_free(response);
        return -1;
    }

    for (int i = 0; i < n_samples; i++) {
        double t = i * response->t_step;
        response->time[i] = t;
        response->input[i] = (i == 0) ? 1.0 / response->t_step : 0.0; /* approx impulse */

        double h = 0.0;
        for (int k = 0; k < n_poles; k++) {
            h += residues[k] * creal(cexp(p[k] * t));
        }
        response->output[i] = h;
    }
    return 0;
}

/*----------------------------------------------------------------------
 * L6 — Steady-State Error
 *----------------------------------------------------------------------*/

int steady_state_error(const TransferFunction *G,
                        ResponseType input_type,
                        double *ess)
{
    if (!G || !ess) { errno = EINVAL; return -1; }

    /* Error constants for unity feedback */
    LaplaceRational L;
    memset(&L, 0, sizeof(L));
    L.numerator.order   = G->num_order;
    L.denominator.order = G->den_order;
    for (int i = 0; i <= G->num_order; i++) L.numerator.coeff[i]   = G->num_coeff[i];
    for (int i = 0; i <= G->den_order; i++) L.denominator.coeff[i] = G->den_coeff[i];

    /* Evaluate at s=0 */
    double complex s0 = 0.0;
    double complex G0 = laplace_rational_eval(&L, s0);

    switch (input_type) {
    case RESPONSE_STEP: {
        /* K_p = lim_{s→0} G(s), e_ss = 1/(1+K_p) */
        if (cabs(G0) > 1e-12) {
            double Kp = creal(G0) * G->gain; /* eval includes gain scaling? we add it */
            /* Actually laplace_rational_eval doesn't multiply by gain */
            /* TransferFunction has gain field; LaplaceRational doesn't. */
            /* G0 is N(0)/D(0) without gain. */
            Kp = creal(G0);
            if (G->sampling_period <= 0.0) {
                /* CT: add gain */
                Kp *= G->gain;
            }
            *ess = 1.0 / (1.0 + Kp);
        } else {
            /* Type ≥ 1: zero steady-state error */
            *ess = 0.0;
        }
        break;
    }
    case RESPONSE_RAMP: {
        /* K_v = lim_{s→0} s·G(s) */
        /* s·G(s) = s·N(s)/D(s) */
        /* If D(0) ≠ 0: s·G(0) = 0 → K_v = 0 → e_ss = ∞ */
        /* If D(0) = 0 (type ≥ 1): K_v = G'(0) */
        if (cabs(G0) < 1e-12) {
            /* D(0)=0, use l'Hôpital */
            LaplacePolynomial dD;
            laplace_poly_derivative(&L.denominator, &dD);
            double complex N0 = laplace_poly_eval(&L.numerator, s0);
            double complex Dp0 = laplace_poly_eval(&dD, s0);
            if (cabs(Dp0) > 1e-12) {
                double Kv = creal(N0 / Dp0) * G->gain;
                *ess = (Kv > 0.0) ? 1.0 / Kv : 1e300;
            } else {
                *ess = 0.0;  /* Type ≥ 2: zero ramp error */
            }
        } else {
            *ess = 1e300;  /* Type 0: infinite ramp error */
        }
        break;
    }
    case RESPONSE_IMPULSE:
    default:
        /* Impulse: no steady-state (approaches 0 for stable systems) */
        *ess = 0.0;
        break;
    }
    return 0;
}

/*----------------------------------------------------------------------
 * Cleanup
 *----------------------------------------------------------------------*/

void time_response_free(TimeResponse *response)
{
    if (response) {
        free(response->time);
        free(response->output);
        free(response->input);
        memset(response, 0, sizeof(*response));
    }
}
