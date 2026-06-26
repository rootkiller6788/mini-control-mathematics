/**
 * @file transform_pairs.c
 * @brief Common Laplace and Z-Transform Pairs — Lookup Table Implementation
 *
 * Implements standard transform pair lookups, custom pair construction,
 * time signal generation for pairs.
 */

#include "transform_pairs.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*----------------------------------------------------------------------
 * L3 — Laplace Pair Lookup
 *----------------------------------------------------------------------*/

int laplace_pair_get(LaplacePairID id, double param, LaplacePairDef *pair)
{
    if (!pair) { errno = EINVAL; return -1; }
    memset(pair, 0, sizeof(*pair));
    pair->id = id;

    LaplaceRational *F = &pair->transform;

    switch (id) {
    case LAPLACE_PAIR_IMPULSE:
        pair->time_name   = "δ(t)";
        pair->s_domain_name = "1";
        F->numerator.order = 0; F->numerator.coeff[0] = 1.0;
        F->denominator.order = 0; F->denominator.coeff[0] = 1.0;
        break;
    case LAPLACE_PAIR_STEP:
        pair->time_name   = "u(t)";
        pair->s_domain_name = "1/s";
        F->numerator.order = 0; F->numerator.coeff[0] = 1.0;
        F->denominator.order = 1; F->denominator.coeff[0] = 0.0;
        F->denominator.coeff[1] = 1.0;
        break;
    case LAPLACE_PAIR_RAMP:
        pair->time_name   = "t·u(t)";
        pair->s_domain_name = "1/s²";
        F->numerator.order = 0; F->numerator.coeff[0] = 1.0;
        F->denominator.order = 2; F->denominator.coeff[0] = 0.0;
        F->denominator.coeff[1] = 0.0; F->denominator.coeff[2] = 1.0;
        break;
    case LAPLACE_PAIR_PARABOLA:
        pair->time_name   = "t²·u(t)";
        pair->s_domain_name = "2/s³";
        F->numerator.order = 0; F->numerator.coeff[0] = 2.0;
        F->denominator.order = 3; F->denominator.coeff[0] = 0.0;
        F->denominator.coeff[1] = 0.0; F->denominator.coeff[2] = 0.0;
        F->denominator.coeff[3] = 1.0;
        break;
    case LAPLACE_PAIR_EXP_DECAY:
        if (param <= 0.0) return -1;
        pair->time_name   = "e^(-at)·u(t)";
        pair->s_domain_name = "1/(s+a)";
        F->numerator.order = 0; F->numerator.coeff[0] = 1.0;
        F->denominator.order = 1; F->denominator.coeff[0] = param;
        F->denominator.coeff[1] = 1.0;
        break;
    case LAPLACE_PAIR_EXP_RAMP:
        if (param <= 0.0) return -1;
        pair->time_name   = "t·e^(-at)·u(t)";
        pair->s_domain_name = "1/(s+a)²";
        F->numerator.order = 0; F->numerator.coeff[0] = 1.0;
        F->denominator.order = 2; F->denominator.coeff[0] = param * param;
        F->denominator.coeff[1] = 2.0 * param; F->denominator.coeff[2] = 1.0;
        break;
    case LAPLACE_PAIR_SINE:
        pair->time_name   = "sin(ωt)·u(t)";
        pair->s_domain_name = "ω/(s²+ω²)";
        F->numerator.order = 0; F->numerator.coeff[0] = param;
        F->denominator.order = 2; F->denominator.coeff[0] = param * param;
        F->denominator.coeff[1] = 0.0; F->denominator.coeff[2] = 1.0;
        break;
    case LAPLACE_PAIR_COSINE:
        pair->time_name   = "cos(ωt)·u(t)";
        pair->s_domain_name = "s/(s²+ω²)";
        F->numerator.order = 1; F->numerator.coeff[0] = 0.0;
        F->numerator.coeff[1] = 1.0;
        F->denominator.order = 2; F->denominator.coeff[0] = param * param;
        F->denominator.coeff[1] = 0.0; F->denominator.coeff[2] = 1.0;
        break;
    case LAPLACE_PAIR_DAMPED_SINE:
        if (param <= 0.0) return -1;
        /* e^(-at)sin(ωt). param interpretation: a used directly, ω from... */
        /* We reuse param only for a; ω is hardcoded to 1.0 for simplicity, */
        /* or use secondary parameter interpretation: param = a, ω = 1.0 */
        pair->time_name   = "e^(-at)sin(t)·u(t)";
        pair->s_domain_name = "1/((s+a)²+1)";
        F->numerator.order = 0; F->numerator.coeff[0] = 1.0;
        F->denominator.order = 2;
        F->denominator.coeff[0] = param * param + 1.0;
        F->denominator.coeff[1] = 2.0 * param;
        F->denominator.coeff[2] = 1.0;
        break;
    case LAPLACE_PAIR_DAMPED_COSINE:
        if (param <= 0.0) return -1;
        pair->time_name   = "e^(-at)cos(t)·u(t)";
        pair->s_domain_name = "(s+a)/((s+a)²+1)";
        F->numerator.order = 1; F->numerator.coeff[0] = param;
        F->numerator.coeff[1] = 1.0;
        F->denominator.order = 2;
        F->denominator.coeff[0] = param * param + 1.0;
        F->denominator.coeff[1] = 2.0 * param;
        F->denominator.coeff[2] = 1.0;
        break;
    case LAPLACE_PAIR_SINH:
        pair->time_name   = "sinh(at)·u(t)";
        pair->s_domain_name = "a/(s²-a²)";
        F->numerator.order = 0; F->numerator.coeff[0] = param;
        F->denominator.order = 2; F->denominator.coeff[0] = -param * param;
        F->denominator.coeff[1] = 0.0; F->denominator.coeff[2] = 1.0;
        break;
    case LAPLACE_PAIR_COSH:
        pair->time_name   = "cosh(at)·u(t)";
        pair->s_domain_name = "s/(s²-a²)";
        F->numerator.order = 1; F->numerator.coeff[0] = 0.0;
        F->numerator.coeff[1] = 1.0;
        F->denominator.order = 2; F->denominator.coeff[0] = -param * param;
        F->denominator.coeff[1] = 0.0; F->denominator.coeff[2] = 1.0;
        break;
    default:
        errno = EINVAL;
        return -1;
    }
    return 0;
}

/*----------------------------------------------------------------------
 * L3 — Z-Transform Pair Lookup
 *----------------------------------------------------------------------*/

int z_pair_get(ZPairID id, double param, double T, ZPairDef *pair)
{
    if (!pair || T <= 0.0) { errno = EINVAL; return -1; }
    memset(pair, 0, sizeof(*pair));
    pair->id = id;

    ZRational *X = &pair->transform;
    X->sampling_period = T;

    switch (id) {
    case Z_PAIR_IMPULSE:
        pair->time_name = "δ[n]";
        pair->z_domain_name = "1";
        X->numerator.order = 0; X->numerator.coeff[0] = 1.0;
        X->denominator.order = 0; X->denominator.coeff[0] = 1.0;
        break;
    case Z_PAIR_STEP:
        pair->time_name = "u[n]";
        pair->z_domain_name = "1/(1-z⁻¹)";
        X->numerator.order = 0; X->numerator.coeff[0] = 1.0;
        X->denominator.order = 1; X->denominator.coeff[0] = 1.0;
        X->denominator.coeff[1] = -1.0;
        break;
    case Z_PAIR_RAMP:
        pair->time_name = "n·u[n]";
        pair->z_domain_name = "z⁻¹/(1-z⁻¹)²";
        X->numerator.order = 1; X->numerator.coeff[0] = 0.0;
        X->numerator.coeff[1] = 1.0;
        X->denominator.order = 2; X->denominator.coeff[0] = 1.0;
        X->denominator.coeff[1] = -2.0; X->denominator.coeff[2] = 1.0;
        break;
    case Z_PAIR_PARABOLA:
        pair->time_name = "n²·u[n]";
        pair->z_domain_name = "z⁻¹(1+z⁻¹)/(1-z⁻¹)³";
        X->numerator.order = 2; X->numerator.coeff[0] = 0.0;
        X->numerator.coeff[1] = 1.0; X->numerator.coeff[2] = 1.0;
        X->denominator.order = 3; X->denominator.coeff[0] = 1.0;
        X->denominator.coeff[1] = -3.0; X->denominator.coeff[2] = 3.0;
        X->denominator.coeff[3] = -1.0;
        break;
    case Z_PAIR_EXP_DECAY:
        if (param <= 0.0) return -1;
        pair->time_name = "aⁿ·u[n]";
        pair->z_domain_name = "1/(1-a·z⁻¹)";
        X->numerator.order = 0; X->numerator.coeff[0] = 1.0;
        X->denominator.order = 1; X->denominator.coeff[0] = 1.0;
        X->denominator.coeff[1] = -param;
        break;
    case Z_PAIR_EXP_RAMP:
        if (param <= 0.0) return -1;
        pair->time_name = "n·aⁿ·u[n]";
        pair->z_domain_name = "a·z⁻¹/(1-a·z⁻¹)²";
        X->numerator.order = 1; X->numerator.coeff[0] = 0.0;
        X->numerator.coeff[1] = param;
        X->denominator.order = 2; X->denominator.coeff[0] = 1.0;
        X->denominator.coeff[1] = -2.0 * param;
        X->denominator.coeff[2] = param * param;
        break;
    case Z_PAIR_SINE: {
        double w = (param > 0) ? param : 1.0;
        pair->time_name = "sin(ωn)·u[n]";
        pair->z_domain_name = "z⁻¹sin(ω)/(1-2z⁻¹cos(ω)+z⁻²)";
        X->numerator.order = 1; X->numerator.coeff[0] = 0.0;
        X->numerator.coeff[1] = sin(w);
        X->denominator.order = 2; X->denominator.coeff[0] = 1.0;
        X->denominator.coeff[1] = -2.0 * cos(w);
        X->denominator.coeff[2] = 1.0;
        break;
    }
    case Z_PAIR_COSINE: {
        double w = (param > 0) ? param : 1.0;
        pair->time_name = "cos(ωn)·u[n]";
        pair->z_domain_name = "(1-z⁻¹cos(ω))/(1-2z⁻¹cos(ω)+z⁻²)";
        X->numerator.order = 1; X->numerator.coeff[0] = 1.0;
        X->numerator.coeff[1] = -cos(w);
        X->denominator.order = 2; X->denominator.coeff[0] = 1.0;
        X->denominator.coeff[1] = -2.0 * cos(w);
        X->denominator.coeff[2] = 1.0;
        break;
    }
    case Z_PAIR_DAMPED_SINE: {
        double a = param;
        double w = 1.0;
        pair->time_name = "aⁿsin(n)·u[n]";
        pair->z_domain_name = "a·z⁻¹sin(ω)/(1-2a·z⁻¹cos(ω)+a²z⁻²)";
        X->numerator.order = 1; X->numerator.coeff[0] = 0.0;
        X->numerator.coeff[1] = a * sin(w);
        X->denominator.order = 2; X->denominator.coeff[0] = 1.0;
        X->denominator.coeff[1] = -2.0 * a * cos(w);
        X->denominator.coeff[2] = a * a;
        break;
    }
    case Z_PAIR_DAMPED_COSINE: {
        double a = param;
        double w = 1.0;
        pair->time_name = "aⁿcos(n)·u[n]";
        pair->z_domain_name = "(1-a·z⁻¹cos(ω))/(1-2a·z⁻¹cos(ω)+a²z⁻²)";
        X->numerator.order = 1; X->numerator.coeff[0] = 1.0;
        X->numerator.coeff[1] = -a * cos(w);
        X->denominator.order = 2; X->denominator.coeff[0] = 1.0;
        X->denominator.coeff[1] = -2.0 * a * cos(w);
        X->denominator.coeff[2] = a * a;
        break;
    }
    default:
        errno = EINVAL;
        return -1;
    }
    return 0;
}

/*----------------------------------------------------------------------
 * L3 — Custom Pair Construction
 *----------------------------------------------------------------------*/

int laplace_pair_custom(int num_order, const double *num_coeff,
                         int den_order, const double *den_coeff,
                         LaplaceRational *result)
{
    if (!num_coeff || !den_coeff || !result ||
        num_order > LAPLACE_MAX_ORDER || den_order > LAPLACE_MAX_ORDER)
    { errno = EINVAL; return -1; }

    memset(result, 0, sizeof(*result));
    result->numerator.order = num_order;
    result->denominator.order = den_order;
    for (int i = 0; i <= num_order; i++) result->numerator.coeff[i] = num_coeff[i];
    for (int i = 0; i <= den_order; i++) result->denominator.coeff[i] = den_coeff[i];
    return 0;
}

/*----------------------------------------------------------------------
 * L3 — Pair Listing
 *----------------------------------------------------------------------*/

int laplace_pair_list(char *buffer)
{
    if (!buffer) { errno = EINVAL; return 0; }
    const char *names[] = {
        "δ(t) ↔ 1",
        "u(t) ↔ 1/s",
        "t·u(t) ↔ 1/s²",
        "t²·u(t) ↔ 2/s³",
        "e^(-at)·u(t) ↔ 1/(s+a)",
        "t·e^(-at)·u(t) ↔ 1/(s+a)²",
        "sin(ωt)·u(t) ↔ ω/(s²+ω²)",
        "cos(ωt)·u(t) ↔ s/(s²+ω²)",
        "e^(-at)sin(ωt)·u(t) ↔ ω/((s+a)²+ω²)",
        "e^(-at)cos(ωt)·u(t) ↔ (s+a)/((s+a)²+ω²)",
        "sinh(at)·u(t) ↔ a/(s²-a²)",
        "cosh(at)·u(t) ↔ s/(s²-a²)"
    };
    int n = sizeof(names) / sizeof(names[0]);
    int pos = 0;
    for (int i = 0; i < n; i++)
        pos += sprintf(buffer + pos, "%d: %s\n", i, names[i]);
    buffer[pos] = '\0';
    return n;
}

int z_pair_list(char *buffer)
{
    if (!buffer) { errno = EINVAL; return 0; }
    const char *names[] = {
        "δ[n] ↔ 1",
        "u[n] ↔ 1/(1-z⁻¹)",
        "n·u[n] ↔ z⁻¹/(1-z⁻¹)²",
        "n²·u[n] ↔ z⁻¹(1+z⁻¹)/(1-z⁻¹)³",
        "aⁿ·u[n] ↔ 1/(1-a·z⁻¹)",
        "n·aⁿ·u[n] ↔ a·z⁻¹/(1-a·z⁻¹)²",
        "sin(ωn)·u[n] ↔ z⁻¹sin(ω)/(1-2z⁻¹cos(ω)+z⁻²)",
        "cos(ωn)·u[n] ↔ (1-z⁻¹cos(ω))/(1-2z⁻¹cos(ω)+z⁻²)",
        "aⁿsin(ωn)·u[n] ↔ a·z⁻¹sin(ω)/(1-2a·z⁻¹cos(ω)+a²z⁻²)",
        "aⁿcos(ωn)·u[n] ↔ (1-a·z⁻¹cos(ω))/(1-2a·z⁻¹cos(ω)+a²z⁻²)"
    };
    int n = sizeof(names) / sizeof(names[0]);
    int pos = 0;
    for (int i = 0; i < n; i++)
        pos += sprintf(buffer + pos, "%d: %s\n", i, names[i]);
    buffer[pos] = '\0';
    return n;
}

/*----------------------------------------------------------------------
 * L3 — Time Signal Generation from Pair
 *----------------------------------------------------------------------*/

int laplace_pair_time_signal(LaplacePairID id, double param,
                              int n_points, double T_max,
                              TimeSignal *sig)
{
    if (!sig || n_points < 2 || T_max <= 0.0) { errno = EINVAL; return -1; }

    sig->n_samples = n_points;
    sig->t_start   = 0.0;
    sig->t_step    = T_max / (n_points - 1);
    sig->values    = (double*)malloc(n_points * sizeof(double));
    if (!sig->values) return -1;

    for (int i = 0; i < n_points; i++) {
        double t = i * sig->t_step;
        double val = 0.0;
        switch (id) {
        case LAPLACE_PAIR_IMPULSE:  val = (i == 0) ? 1e10 : 0.0; break; /* approximate */
        case LAPLACE_PAIR_STEP:     val = 1.0; break;
        case LAPLACE_PAIR_RAMP:     val = t; break;
        case LAPLACE_PAIR_PARABOLA: val = t * t; break;
        case LAPLACE_PAIR_EXP_DECAY: val = exp(-param * t); break;
        case LAPLACE_PAIR_EXP_RAMP:  val = t * exp(-param * t); break;
        case LAPLACE_PAIR_SINE:      val = sin(param * t); break;
        case LAPLACE_PAIR_COSINE:    val = cos(param * t); break;
        case LAPLACE_PAIR_DAMPED_SINE:  val = exp(-param * t) * sin(t); break;
        case LAPLACE_PAIR_DAMPED_COSINE: val = exp(-param * t) * cos(t); break;
        case LAPLACE_PAIR_SINH:     val = sinh(param * t); break;
        case LAPLACE_PAIR_COSH:     val = cosh(param * t); break;
        default: val = 0.0; break;
        }
        sig->values[i] = val;
    }
    return 0;
}

int z_pair_sequence(ZPairID id, double param, int length,
                     DiscreteSignal *seq)
{
    if (!seq || length < 1) { errno = EINVAL; return -1; }

    seq->length      = length;
    seq->start_index = 0;
    seq->values      = (double*)malloc(length * sizeof(double));
    if (!seq->values) return -1;

    for (int n = 0; n < length; n++) {
        double val = 0.0;
        switch (id) {
        case Z_PAIR_IMPULSE: val = (n == 0) ? 1.0 : 0.0; break;
        case Z_PAIR_STEP:    val = 1.0; break;
        case Z_PAIR_RAMP:    val = (double)n; break;
        case Z_PAIR_PARABOLA: val = (double)(n * n); break;
        case Z_PAIR_EXP_DECAY: val = pow(param, n); break;
        case Z_PAIR_EXP_RAMP:  val = n * pow(param, n); break;
        case Z_PAIR_SINE:      val = sin(param * n); break;
        case Z_PAIR_COSINE:    val = cos(param * n); break;
        case Z_PAIR_DAMPED_SINE:  val = pow(param, n) * sin((double)n); break;
        case Z_PAIR_DAMPED_COSINE: val = pow(param, n) * cos((double)n); break;
        default: val = 0.0; break;
        }
        seq->values[n] = val;
    }
    return 0;
}
