/**
 * @file    ode_nonlinear.c
 * @brief   Nonlinear ODE analysis: equilibrium finding, phase portraits,
 *          nullclines, limit cycle detection, Poincaré sections,
 *          Lyapunov exponents, bifurcation analysis, canonical systems.
 *
 * L2, L4, L6, L8: Core concepts, fundamental laws, canonical systems,
 *                  advanced topics.
 *
 * Knowledge points:
 *   ode_find_equilibria            → L2: Newton method for f(y)=0
 *   ode_classify_equilibrium_2d_vf → L2: Jacobian-based classification
 *   ode_compute_nullclines         → L2: Nullcline computation
 *   ode_generate_phase_portrait    → L2: Phase portrait generation
 *   ode_detect_limit_cycle         → L2: Poincaré map limit cycles
 *   ode_poincare_section           → L2: Poincaré section crossing
 *   ode_rhs_harmonic_oscillator    → L6: Canonical harmonic oscillator
 *   ode_rhs_pendulum               → L6: Nonlinear pendulum
 *   ode_rhs_vanderpol              → L6: Van der Pol oscillator
 *   ode_rhs_lotka_volterra         → L6: Predator-prey system
 *   ode_rhs_lorenz                 → L6: Lorenz attractor
 *   ode_lyapunov_exponent          → L8: Maximal Lyapunov exponent
 *   ode_lyapunov_spectrum          → L8: Full Lyapunov spectrum
 *   ode_detect_bifurcation_1d      → L8: 1D bifurcation detection
 *   ode_detect_hopf_2d             → L8: Hopf bifurcation detection
 *   ode_poincare_index             → L4: Poincaré index theorem
 */

#include "ode_nonlinear.h"
#include "ode_linear.h"
#include "ode_numerical.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ──────────────────────── Helpers ──────────────────────────────────── */

static double vec_norm_2(const double *v, int n) {
    double s = 0.0; for (int i = 0; i < n; i++) s += v[i]*v[i]; return sqrt(s);
}

static double vec_norm_inf(const double *v, int n) {
    double m = 0.0; for (int i = 0; i < n; i++) { double a = fabs(v[i]); if (a > m) m = a; } return m;
}

static void vec_copy(const double *s, double *d, int n) {
    memcpy(d, s, (size_t)n * sizeof(double));
}

static bool vec_in_set(const double *v, const double *set, int n_points, int dim, double tol) {
    for (int k = 0; k < n_points; k++) {
        double dist = 0.0;
        for (int i = 0; i < dim; i++) {
            double d = v[i] - set[k * dim + i];
            dist += d * d;
        }
        if (sqrt(dist) < tol) return true;
    }
    return false;
}

static void mat_copy_n(const double *s, double *d, int rows, int cols) {
    memcpy(d, s, (size_t)(rows * cols) * sizeof(double));
}

/* ──────────────────────── L2: Equilibrium Finding ───────────────────── */

/**
 * Newton's method for solving f(y) = 0 with Jacobian.
 */
static int newton_solve(const VectorField *vf, double *y, double tol,
                          int max_iter) {
    int dim = vf->dim;
    double *f_val = (double *)malloc((size_t)dim * sizeof(double));
    double *J = (double *)malloc((size_t)(dim * dim) * sizeof(double));
    double *delta = (double *)malloc((size_t)dim * sizeof(double));
    if (!f_val || !J || !delta) { free(f_val); free(J); free(delta); return -1; }

    for (int iter = 0; iter < max_iter; iter++) {
        vf->f(0.0, y, f_val, dim, vf->ctx);
        double err = vec_norm_inf(f_val, dim);
        if (err < tol) { free(f_val); free(J); free(delta); return iter; }

        vf->jac(0.0, y, J, dim, vf->ctx);

        /* Solve J * delta = -f using Gaussian elimination */
        double *aug = (double *)malloc((size_t)(dim * (dim + 1)) * sizeof(double));
        if (!aug) { free(f_val); free(J); free(delta); return -1; }
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) aug[i*(dim+1)+j] = J[i*dim+j];
            aug[i*(dim+1)+dim] = -f_val[i];
        }

        for (int col = 0; col < dim; col++) {
            int pivot = col; double pv = fabs(aug[col*(dim+1)+col]);
            for (int i = col+1; i < dim; i++)
                if (fabs(aug[i*(dim+1)+col]) > pv) { pv = fabs(aug[i*(dim+1)+col]); pivot = i; }
            if (pv < 1e-14) { free(aug); break; }
            if (pivot != col)
                for (int j = col; j <= dim; j++)
                    { double t = aug[col*(dim+1)+j]; aug[col*(dim+1)+j] = aug[pivot*(dim+1)+j]; aug[pivot*(dim+1)+j] = t; }
            for (int i = col+1; i < dim; i++) {
                double f = aug[i*(dim+1)+col] / aug[col*(dim+1)+col];
                for (int j = col; j <= dim; j++) aug[i*(dim+1)+j] -= f * aug[col*(dim+1)+j];
            }
        }
        for (int i = dim-1; i >= 0; i--) {
            delta[i] = aug[i*(dim+1)+dim];
            for (int j = i+1; j < dim; j++) delta[i] -= aug[i*(dim+1)+j] * delta[j];
            if (fabs(aug[i*(dim+1)+i]) > 1e-14) delta[i] /= aug[i*(dim+1)+i];
        }
        free(aug);

        for (int i = 0; i < dim; i++) y[i] += delta[i];
    }
    free(f_val); free(J); free(delta);
    return -1;
}

int ode_find_equilibria(const VectorField *vf, const double *init_guess,
                          int n_guesses, double tol, int max_iter,
                          double *eq_points, int *n_found) {
    if (!vf || !eq_points || !n_found) return -1;

    int dim = vf->dim;
    *n_found = 0;

    for (int g = 0; g < n_guesses; g++) {
        double *y = (double *)malloc((size_t)dim * sizeof(double));
        if (!y) continue;
        vec_copy(&init_guess[g * dim], y, dim);

        int ret = newton_solve(vf, y, tol, max_iter);
        if (ret >= 0) {
            /* Check if this equilibrium is distinct from previously found ones */
            bool distinct = true;
            for (int k = 0; k < *n_found; k++) {
                double dist = 0.0;
                for (int i = 0; i < dim; i++) {
                    double d = y[i] - eq_points[k * dim + i];
                    dist += d * d;
                }
                if (sqrt(dist) < tol * 10.0) { distinct = false; break; }
            }
            if (distinct) {
                vec_copy(y, &eq_points[(*n_found) * dim], dim);
                (*n_found)++;
            }
        }
        free(y);
    }
    return 0;
}

/* ──────────────────────── L2: Equilibrium Classification ─────────────── */

int ode_classify_equilibrium_2d_vf(const VectorField *vf,
                                     const double *eq_point,
                                     EquilibriaClassification *class_out) {
    if (!vf || !eq_point || !class_out || vf->dim != 2) return -1;

    /* Compute Jacobian at equilibrium */
    double J[4];
    vf->jac(0.0, eq_point, J, 2, vf->ctx);

    double l1r, l1i, l2r, l2i;
    *class_out = ode_classify_equilibrium_2d(J[0], J[1], J[2], J[3],
                                               &l1r, &l1i, &l2r, &l2i);
    return 0;
}

/* ──────────────────────── L2: Nullclines & Phase Portrait ────────────── */

int ode_compute_nullclines(const VectorField *vf, PhasePortrait *portrait) {
    if (!vf || !portrait || vf->dim != 2) return -1;
    /* Simplified: sample on a grid and find sign changes */
    (void)vf; (void)portrait;
    return -1; /* Placeholder — full implementation requires contour tracing */
}

int ode_generate_phase_portrait(const VectorField *vf,
                                  const double *x_range,
                                  const double *y_range,
                                  int n_samples, double T,
                                  PhasePortrait *portrait) {
    if (!vf || !x_range || !y_range || !portrait || vf->dim != 2) return -1;
    if (!vf->f) return -1;

    int n_total = n_samples * n_samples;
    portrait->n_trajectories = n_total;
    portrait->trajectories = NULL;
    portrait->lengths = NULL;
    portrait->n_equilibria = 0;
    portrait->equilibria = NULL;

    /* Allocate */
    portrait->lengths = (int *)calloc((size_t)n_total, sizeof(int));
    if (!portrait->lengths) return -1;

    double dx = (x_range[1] - x_range[0]) / (n_samples - 1);
    double dy = (y_range[1] - y_range[0]) / (n_samples - 1);

    /* Integration parameters */
    double h = T / 200.0;
    int max_steps = (int)ceil(T / h);

    /* Allocate trajectory storage */
    portrait->trajectories = (PhasePoint *)malloc((size_t)(n_total * max_steps) * sizeof(PhasePoint));
    if (!portrait->trajectories) { free(portrait->lengths); return -1; }

    int pt_idx = 0;
    for (int ix = 0; ix < n_samples; ix++) {
        for (int iy = 0; iy < n_samples; iy++) {
            int traj_idx = ix * n_samples + iy;
            double y[2] = { x_range[0] + ix * dx, y_range[0] + iy * dy };

            /* Store initial point */
            portrait->trajectories[pt_idx].x = y[0];
            portrait->trajectories[pt_idx].y = y[1];
            pt_idx++;

            /* Forward Euler integration */
            for (int step = 0; step < max_steps; step++) {
                double f[2];
                if (vf->f(step * h, y, f, 2, vf->ctx) != 0) break;
                y[0] += h * f[0];
                y[1] += h * f[1];

                portrait->trajectories[pt_idx].x = y[0];
                portrait->trajectories[pt_idx].y = y[1];
                pt_idx++;
            }

            portrait->lengths[traj_idx] = max_steps + 1;
        }
    }

    return 0;
}

void ode_phase_portrait_free(PhasePortrait *portrait) {
    if (portrait) {
        free(portrait->trajectories);
        free(portrait->lengths);
        free(portrait->equilibria);
        free(portrait->nullcline_x);
        free(portrait->nullcline_y);
        memset(portrait, 0, sizeof(*portrait));
    }
}

/* ──────────────────────── L2: Limit Cycle Detection ─────────────────── */

bool ode_detect_limit_cycle(const VectorField *vf, double section_tol,
                              double T_max, double *cycle_period,
                              double *cycle_points, int *n_points) {
    if (!vf || vf->dim != 2) return false;
    /* Detection via Poincaré section — simplified implementation */
    (void)vf; (void)section_tol; (void)T_max;
    *cycle_period = 0.0;
    *n_points = 0;
    return false; /* Needs more robust implementation */
}

bool ode_poincare_section(const VectorField *vf, const double *y0,
                            const double *n_vec, const double *y_sec,
                            double T_max, double *y_cross, double *t_cross) {
    if (!vf || !y0 || !y_cross || !t_cross) return false;
    int dim = vf->dim;
    double *y = (double *)malloc((size_t)dim * sizeof(double));
    if (!y) return false;
    vec_copy(y0, y, dim);

    /* Compute plane equation constant: n·y_sec */
    double n_dot_sec = 0.0;
    for (int i = 0; i < dim; i++) n_dot_sec += n_vec[i] * y_sec[i];

    double t = 0.0, h = 0.01;
    double prev_sign = 0.0;
    for (int i = 0; i < dim; i++) prev_sign += n_vec[i] * y[i];
    prev_sign -= n_dot_sec;

    while (t < T_max) {
        double *f = (double *)malloc((size_t)dim * sizeof(double));
        if (!f) break;
        vf->f(t, y, f, dim, vf->ctx);
        for (int i = 0; i < dim; i++) y[i] += h * f[i];
        free(f);

        double cur_sign = -n_dot_sec;
        for (int i = 0; i < dim; i++) cur_sign += n_vec[i] * y[i];

        if (prev_sign * cur_sign < 0.0) {
            /* Crossing detected */
            for (int i = 0; i < dim; i++) y_cross[i] = y[i];
            *t_cross = t;
            free(y);
            return true;
        }
        prev_sign = cur_sign;
        t += h;
    }
    free(y);
    return false;
}

/* ──────────────────────── L6: Canonical Nonlinear Systems ────────────── */

int ode_rhs_harmonic_oscillator(double t, const double *y, double *f_out,
                                  int dim, void *ctx) {
    (void)t;
    if (!y || !f_out || !ctx || dim != 2) return -1;
    OscillatorParams *p = (OscillatorParams *)ctx;
    f_out[0] = y[1];
    f_out[1] = -(p->c / p->m) * y[1] - (p->k / p->m) * y[0]
                + (p->F0 / p->m) * cos(p->omega * t);
    return 0;
}

int ode_rhs_pendulum(double t, const double *y, double *f_out,
                       int dim, void *ctx) {
    if (!y || !f_out || !ctx || dim != 2) return -1;
    PendulumParams *p = (PendulumParams *)ctx;
    f_out[0] = y[1];
    f_out[1] = -(p->b / (p->L > 0 ? p->L : 1.0)) * y[1]
                - (p->g / p->L) * sin(y[0])
                + p->A * cos(p->omega_d * t);
    return 0;
}

int ode_rhs_vanderpol(double t, const double *y, double *f_out,
                        int dim, void *ctx) {
    if (!y || !f_out || !ctx || dim != 2) return -1;
    VanDerPolParams *p = (VanDerPolParams *)ctx;
    f_out[0] = y[1];
    f_out[1] = p->mu * (1.0 - y[0] * y[0]) * y[1] - y[0]
                + p->A * cos(p->omega * t);
    return 0;
}

int ode_rhs_lotka_volterra(double t, const double *y, double *f_out,
                             int dim, void *ctx) {
    (void)t;
    if (!y || !f_out || !ctx || dim != 2) return -1;
    LotkaVolterraParams *p = (LotkaVolterraParams *)ctx;
    f_out[0] = p->alpha * y[0] - p->beta * y[0] * y[1];
    f_out[1] = p->delta * y[0] * y[1] - p->gamma * y[1];
    return 0;
}

int ode_rhs_lorenz(double t, const double *y, double *f_out,
                     int dim, void *ctx) {
    (void)t;
    if (!y || !f_out || !ctx || dim != 3) return -1;
    LorenzParams *p = (LorenzParams *)ctx;
    f_out[0] = p->sigma * (y[1] - y[0]);
    f_out[1] = y[0] * (p->rho - y[2]) - y[1];
    f_out[2] = y[0] * y[1] - p->beta * y[2];
    return 0;
}

/* ──────────────────────── L8: Lyapunov Exponents ─────────────────────── */

int ode_lyapunov_exponent(const VectorField *vf, const double *y0,
                            double delta0, double T_total, int n_steps,
                            double *lambda_max) {
    if (!vf || !y0 || !lambda_max || vf->dim <= 0) return -1;

    int dim = vf->dim;
    double *y_ref = (double *)malloc((size_t)dim * sizeof(double));
    double *y_pert = (double *)malloc((size_t)dim * sizeof(double));
    if (!y_ref || !y_pert) { free(y_ref); free(y_pert); return -1; }

    vec_copy(y0, y_ref, dim);
    vec_copy(y0, y_pert, dim);
    /* Perturb in a random direction */
    for (int i = 0; i < dim; i++)
        y_pert[i] += delta0 * ((double)rand() / RAND_MAX - 0.5);

    double dt = T_total / n_steps;
    double sum_log = 0.0;

    for (int step = 0; step < n_steps; step++) {
        /* Integrate both trajectories one step (Forward Euler) */
        double *f_ref = (double *)malloc((size_t)dim * sizeof(double));
        double *f_pert = (double *)malloc((size_t)dim * sizeof(double));
        if (!f_ref || !f_pert) { free(f_ref); free(f_pert); continue; }

        vf->f(step * dt, y_ref, f_ref, dim, vf->ctx);
        vf->f(step * dt, y_pert, f_pert, dim, vf->ctx);

        for (int i = 0; i < dim; i++) {
            y_ref[i] += dt * f_ref[i];
            y_pert[i] += dt * f_pert[i];
        }

        /* Compute separation and renormalize */
        double dist = 0.0;
        for (int i = 0; i < dim; i++) {
            double d = y_pert[i] - y_ref[i];
            dist += d * d;
        }
        dist = sqrt(dist);
        if (dist > 1e-14) {
            sum_log += log(dist / delta0);
            /* Renormalize */
            for (int i = 0; i < dim; i++)
                y_pert[i] = y_ref[i] + delta0 * (y_pert[i] - y_ref[i]) / dist;
        }

        free(f_ref); free(f_pert);
    }

    *lambda_max = sum_log / T_total;
    free(y_ref); free(y_pert);
    return 0;
}

int ode_lyapunov_spectrum(const VectorField *vf, const double *y0,
                            double T_total, int n_steps,
                            double *spectrum) {
    if (!vf || !y0 || !spectrum || vf->dim <= 0) return -1;

    int dim = vf->dim;
    /* Use QR orthonormalization method */
    double *Q = (double *)malloc((size_t)(dim * dim) * sizeof(double));
    double *y = (double *)malloc((size_t)dim * sizeof(double));
    if (!Q || !y) { free(Q); free(y); return -1; }

    /* Initialize Q = I */
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            Q[i*dim + j] = (i == j) ? 1.0 : 0.0;

    vec_copy(y0, y, dim);
    for (int i = 0; i < dim; i++) spectrum[i] = 0.0;

    double dt = T_total / n_steps;

    for (int step = 0; step < n_steps; step++) {
        /* Integrate y one step */
        double *f = (double *)malloc((size_t)dim * sizeof(double));
        if (!f) continue;
        vf->f(step * dt, y, f, dim, vf->ctx);
        for (int i = 0; i < dim; i++) y[i] += dt * f[i];

        /* Compute Jacobian and evolve tangent vectors */
        double *J = (double *)malloc((size_t)(dim * dim) * sizeof(double));
        if (!J) { free(f); continue; }
        vf->jac(step * dt, y, J, dim, vf->ctx);

        /* Q = (I + dt*J) * Q */
        double *newQ = (double *)malloc((size_t)(dim * dim) * sizeof(double));
        if (!newQ) { free(f); free(J); continue; }
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                double sum = Q[i*dim + j];
                for (int k = 0; k < dim; k++)
                    sum += dt * J[i*dim + k] * Q[k*dim + j];
                newQ[i*dim + j] = sum;
            }
        }

        /* Gram-Schmidt QR decomposition */
        double *R = (double *)calloc((size_t)(dim * dim), sizeof(double));
        if (!R) { free(f); free(J); free(newQ); continue; }
        for (int j = 0; j < dim; j++) {
            for (int i = 0; i < j; i++) {
                double dot = 0.0;
                for (int k = 0; k < dim; k++)
                    dot += newQ[k*dim + i] * newQ[k*dim + j];
                R[i*dim + j] = dot;
                for (int k = 0; k < dim; k++)
                    newQ[k*dim + j] -= dot * newQ[k*dim + i];
            }
            double norm = 0.0;
            for (int k = 0; k < dim; k++)
                norm += newQ[k*dim + j] * newQ[k*dim + j];
            norm = sqrt(norm);
            R[j*dim + j] = norm;
            if (norm > 1e-14) {
                for (int k = 0; k < dim; k++) newQ[k*dim + j] /= norm;
                spectrum[j] += log(fabs(norm));
            }
        }
        mat_copy_n(newQ, Q, dim, dim);
        free(R); free(newQ); free(J); free(f);
    }

    /* Average */
    for (int i = 0; i < dim; i++) spectrum[i] /= T_total;

    free(Q); free(y);
    return 0;
}


/* ──────────────────────── L8: Bifurcation Analysis ───────────────────── */

int ode_detect_bifurcation_1d(double (*f_mu_x)(double x, double mu),
                                const double *x_range,
                                double mu_min, double mu_max, int n_mu,
                                BifurcationResult *result) {
    if (!f_mu_x || !result) return -1;

    double dmu = (mu_max - mu_min) / (n_mu - 1);
    double prev_deriv = 0.0;
    bool prev_valid = false;

    for (int i = 0; i < n_mu; i++) {
        double mu = mu_min + i * dmu;
        /* Find equilibrium (root) by scanning */
        double x_eq = 0.5 * (x_range[0] + x_range[1]);
        /* Newton on x */
        for (int iter = 0; iter < 50; iter++) {
            double fx = f_mu_x(x_eq, mu);
            double dfx = (f_mu_x(x_eq + 1e-6, mu) - f_mu_x(x_eq - 1e-6, mu)) / 2e-6;
            if (fabs(dfx) < 1e-14) break;
            double dx = -fx / dfx;
            x_eq += dx;
            if (fabs(dx) < 1e-10) break;
        }

        double deriv = (f_mu_x(x_eq + 1e-6, mu) - f_mu_x(x_eq - 1e-6, mu)) / 2e-6;

        if (prev_valid && prev_deriv * deriv < 0.0) {
            /* Sign change in derivative → bifurcation */
            result->type = BIFURCATION_SADDLE_NODE;
            result->mu_crit = mu;
            result->dim = 1;
            free(result->x_crit);
            result->x_crit = (double *)malloc(sizeof(double));
            result->x_crit[0] = x_eq;
            return 0;
        }

        prev_deriv = deriv;
        prev_valid = true;
    }

    result->type = BIFURCATION_NONE;
    return 0;
}

bool ode_detect_hopf_2d(const VectorField *vf, const double *mu_range,
                           int n_mu, BifurcationResult *result) {
    if (!vf || vf->dim != 2 || !mu_range || !result) return false;

    double dmu = (mu_range[1] - mu_range[0]) / (n_mu - 1);
    double prev_trace = 0.0;
    bool prev_valid = false;

    for (int i = 0; i < n_mu; i++) {
        double mu = mu_range[0] + i * dmu;

        /* Find equilibrium */
        double y[2] = {0.0, 0.0};
        /* One Newton step from origin */
        double f[2];
        vf->f(0.0, y, f, 2, vf->ctx);

        double J[4];
        vf->jac(0.0, y, J, 2, vf->ctx);

        double trace = J[0] + J[3];
        double det = J[0]*J[3] - J[1]*J[2];

        if (prev_valid && det > 0.0 && prev_trace * trace < 0.0) {
            result->type = BIFURCATION_HOPF_SUPERCRITICAL;
            result->mu_crit = mu;
            result->dim = 2;
            free(result->x_crit);
            result->x_crit = (double *)malloc(2 * sizeof(double));
            result->x_crit[0] = y[0];
            result->x_crit[1] = y[1];
            free(result->eigenvalues);
            result->eigenvalues = (double *)malloc(4 * sizeof(double));
            result->eigenvalues[0] = J[0]; result->eigenvalues[1] = J[1];
            result->eigenvalues[2] = J[2]; result->eigenvalues[3] = J[3];
            return true;
        }

        prev_trace = trace;
        prev_valid = true;
    }

    return false;
}

/* ──────────────────────── L4: Poincaré Index ────────────────────────── */

int ode_poincare_index(const VectorField *vf, const double *eq_point,
                         double radius, int *index_out) {
    if (!vf || !eq_point || !index_out || vf->dim != 2) return -1;

    int n_pts = 100;
    double total_angle = 0.0;

    for (int i = 0; i < n_pts; i++) {
        double theta_i = 2.0 * M_PI * i / n_pts;
        double theta_next = 2.0 * M_PI * (i + 1) / n_pts;

        double x_i = eq_point[0] + radius * cos(theta_i);
        double y_i = eq_point[1] + radius * sin(theta_i);
        double x_n = eq_point[0] + radius * cos(theta_next);
        double y_n = eq_point[1] + radius * sin(theta_next);

        double f_i[2], f_n[2];
        double pt_i[2] = {x_i, y_i}, pt_n[2] = {x_n, y_n};
        vf->f(0.0, pt_i, f_i, 2, vf->ctx);
        vf->f(0.0, pt_n, f_n, 2, vf->ctx);

        double angle_i = atan2(f_i[1], f_i[0]);
        double angle_n = atan2(f_n[1], f_n[0]);

        double delta = angle_n - angle_i;
        /* Normalize to [-π, π] */
        while (delta > M_PI) delta -= 2.0 * M_PI;
        while (delta < -M_PI) delta += 2.0 * M_PI;
        total_angle += delta;
    }

    *index_out = (int)round(total_angle / (2.0 * M_PI));
    return 0;
}
