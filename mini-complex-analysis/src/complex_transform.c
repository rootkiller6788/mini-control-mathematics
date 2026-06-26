/**
 * @file complex_transform.c
 * @brief Integral transforms: Laplace, Fourier, Z-transforms,
 *        partial fraction expansion
 *
 * L4-L5: Transform theory bridging time/frequency domains
 *
 * Textbook: Oppenheim & Willsky (1996)
 */

#include "complex_transform.h"
#include "complex_number.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

Complex laplace_transform_numeric(const double *t, const double *f,
                                   int n, Complex s) {
    Complex sum = {0.0, 0.0};
    if (!t || !f || n < 2) return sum;
    for (int i = 1; i < n; i++) {
        double dt = t[i] - t[i-1];
        double t_mid = 0.5 * (t[i] + t[i-1]);
        double f_mid = 0.5 * (f[i] + f[i-1]);
        Complex st = {s.re * t_mid, s.im * t_mid};
        Complex e_st = complex_exp((Complex){-st.re, -st.im});
        Complex term = complex_scale(e_st, f_mid * dt);
        sum = complex_add(sum, term);
    }
    return sum;
}

/* Inverse Laplace via Bromwich integral using residue theorem
 * f(t) = (1/2pi*i) * integral_{c-i*inf}^{c+i*inf} F(s) e^{st} ds
 *      = sum of residues of F(s)*e^{st} at poles
 */
double inverse_laplace_rational(const Complex *numer_coeffs, int deg_num,
                                 const Complex *denom_coeffs, int deg_den,
                                 double t) {
    if (!numer_coeffs || !denom_coeffs || deg_den <= 0) return 0.0;
    /* Find poles */
    Complex *poles = (Complex *)calloc(deg_den, sizeof(Complex));
    int n_poles = 0;
    for (int guess = 0; guess < deg_den * 2; guess++) {
        Complex z = {cos(2*M_PI*guess/(deg_den*2)), sin(2*M_PI*guess/(deg_den*2))};
        for (int it = 0; it < 100; it++) {
            Complex fz = {0,0}, dfz = {0,0};
            for (int k = 0; k <= deg_den; k++) {
                Complex zk = complex_pow_int(z, k);
                fz = complex_add(fz, complex_mul(denom_coeffs[k], zk));
                if (k > 0) {
                    Complex kc = {(double)k, 0};
                    Complex zkm1 = complex_pow_int(z, k-1);
                    dfz = complex_add(dfz, complex_mul(complex_mul(kc, denom_coeffs[k]), zkm1));
                }
            }
            if (complex_abs(fz) < 1e-10) {
                int dup = 0;
                for (int j = 0; j < n_poles; j++)
                    if (complex_distance(z, poles[j]) < 1e-8) { dup = 1; break; }
                if (!dup) { poles[n_poles++] = z; }
                break;
            }
            if (complex_abs2(dfz) < 1e-30) break;
            z = complex_sub(z, complex_div(fz, dfz));
        }
    }
    double result = 0.0;
    for (int i = 0; i < n_poles; i++) {
        /* Residue at simple pole p: N(p)/D'(p) * e^{pt} */
        Complex num_val = {0,0}, den_deriv = {0,0};
        for (int k = 0; k <= deg_num; k++) {
            Complex zk = complex_pow_int(poles[i], k);
            num_val = complex_add(num_val, complex_mul(numer_coeffs[k], zk));
        }
        for (int k = 1; k <= deg_den; k++) {
            Complex kc = {(double)k, 0};
            Complex zkm1 = complex_pow_int(poles[i], k-1);
            den_deriv = complex_add(den_deriv,
                complex_mul(complex_mul(kc, denom_coeffs[k]), zkm1));
        }
        Complex residue = complex_div(num_val, den_deriv);
        Complex exp_term;
        exp_term.re = exp(poles[i].re * t) * cos(poles[i].im * t);
        exp_term.im = exp(poles[i].re * t) * sin(poles[i].im * t);
        Complex contrib = complex_mul(residue, exp_term);
        result += contrib.re;
    }
    free(poles);
    return result;
}

/* Numeric Fourier transform */
void fourier_transform_numeric(const double *t, const double *f, int n,
                                double omega_min, double omega_max,
                                int n_freq, Complex *out) {
    if (!t || !f || !out || n < 2 || n_freq <= 0) return;
    for (int k = 0; k < n_freq; k++) {
        double omega = omega_min + (omega_max - omega_min) * k / (n_freq > 1 ? n_freq - 1 : 1);
        Complex s = {0.0, -omega};
        out[k] = laplace_transform_numeric(t, f, n, s);
    }
}


/* Partial fraction expansion of rational function G(s) = N(s)/D(s) */
PartialFraction partial_fraction_expand(const Complex *num, int deg_num,
                                         const Complex *den, int deg_den) {
    PartialFraction pf;
    memset(&pf, 0, sizeof(PartialFraction));
    if (!num || !den || deg_den <= 0) return pf;
    /* Find poles of denominator */
    Complex *poles = (Complex *)calloc(deg_den, sizeof(Complex));
    int *mults = (int *)calloc(deg_den, sizeof(int));
    int n_poles = 0;
    /* Newton for each guess */
    for (int guess = 0; guess < deg_den * 3; guess++) {
        Complex z;
        double r = 1.0 + guess * 0.3;
        double th = 2.0*M_PI * (guess % 17) / 17.0;
        z.re = r * cos(th); z.im = r * sin(th);
        for (int it = 0; it < 100; it++) {
            Complex fz = {0,0}, dfz = {0,0};
            for (int k = 0; k <= deg_den; k++) {
                Complex zk = complex_pow_int(z, k);
                fz = complex_add(fz, complex_mul(den[k], zk));
                if (k > 0) {
                    Complex kc = {(double)k, 0};
                    dfz = complex_add(dfz,
                        complex_mul(complex_mul(kc, den[k]), complex_pow_int(z, k-1)));
                }
            }
            if (complex_abs(fz) < 1e-10) {
                int dup = 0;
                for (int j = 0; j < n_poles; j++)
                    if (complex_distance(z, poles[j]) < 1e-6) { dup = 1; mults[j]++; break; }
                if (!dup) { poles[n_poles] = z; mults[n_poles] = 1; n_poles++; }
                break;
            }
            if (complex_abs2(dfz) < 1e-30) break;
            z = complex_sub(z, complex_div(fz, dfz));
        }
    }
    /* Compute residues for each distinct pole (assuming simple poles) */
    pf.n_poles = n_poles;
    pf.poles = (Complex *)malloc(n_poles * sizeof(Complex));
    pf.mults = (int *)malloc(n_poles * sizeof(int));
    pf.residues = (Complex *)malloc(n_poles * sizeof(Complex));
    if (pf.poles && pf.mults && pf.residues) {
        memcpy(pf.poles, poles, n_poles * sizeof(Complex));
        memcpy(pf.mults, mults, n_poles * sizeof(int));
        for (int i = 0; i < n_poles; i++) {
            Complex num_val = {0,0}, den_deriv = {0,0};
            for (int k = 0; k <= deg_num; k++)
                num_val = complex_add(num_val, complex_mul(num[k], complex_pow_int(poles[i], k)));
            for (int k = 1; k <= deg_den; k++) {
                Complex kc = {(double)k, 0};
                den_deriv = complex_add(den_deriv,
                    complex_mul(complex_mul(kc, den[k]), complex_pow_int(poles[i], k-1)));
            }
            pf.residues[i] = complex_div(num_val, den_deriv);
        }
    }
    free(poles); free(mults);
    return pf;
}

Complex partial_fraction_eval(const PartialFraction *pf, Complex s) {
    Complex sum = {0.0, 0.0};
    if (!pf) return sum;
    for (int i = 0; i < pf->n_poles; i++) {
        Complex diff = complex_sub(s, pf->poles[i]);
        Complex term = complex_div(pf->residues[i], diff);
        sum = complex_add(sum, term);
    }
    /* Add direct/polynomial part */
    for (int k = 0; k < pf->direct_deg; k++) {
        Complex sk = complex_pow_int(s, k);
        Complex term = complex_mul(pf->direct[k], sk);
        sum = complex_add(sum, term);
    }
    return sum;
}

void partial_fraction_free(PartialFraction *pf) {
    if (pf) {
        free(pf->poles);
        free(pf->mults);
        free(pf->residues);
        free(pf->direct);
        memset(pf, 0, sizeof(PartialFraction));
    }
}


/* Z-transform of finite sequence: X(z) = sum_{n=0}^{N-1} x[n] * z^{-n} */
Complex z_transform(const double *x, int n, Complex z) {
    Complex sum = {0.0, 0.0};
    if (!x || n <= 0) return sum;
    Complex z_inv = complex_div((Complex){1,0}, z);
    Complex z_pow = {1.0, 0.0};
    for (int i = 0; i < n; i++) {
        Complex term = complex_scale(z_pow, x[i]);
        sum = complex_add(sum, term);
        z_pow = complex_mul(z_pow, z_inv);
    }
    return sum;
}

/* Inverse Z-transform via residue method:
 * x[k] = sum of residues of X(z)*z^{k-1} at poles inside ROC
 */
double inverse_z_transform_residue(const Complex *numer_coeffs, int deg_num,
                                    const Complex *denom_coeffs, int deg_den,
                                    int k) {
    if (!numer_coeffs || !denom_coeffs || deg_den <= 0) return 0.0;
    /* x[k] = sum Res( X(z)*z^{k-1} ) at poles */
    /* Find poles via Newton */
    Complex *poles = (Complex *)calloc(deg_den, sizeof(Complex));
    int n_poles = 0;
    for (int guess = 0; guess < deg_den * 2 && n_poles < deg_den; guess++) {
        Complex z = {cos(2*M_PI*guess/(deg_den*2)), sin(2*M_PI*guess/(deg_den*2))};
        for (int it = 0; it < 100; it++) {
            Complex fz = {0,0}, dfz = {0,0};
            for (int j = 0; j <= deg_den; j++) {
                Complex zj = complex_pow_int(z, j);
                fz = complex_add(fz, complex_mul(denom_coeffs[j], zj));
                if (j > 0) {
                    Complex jc = {(double)j, 0};
                    dfz = complex_add(dfz, complex_mul(complex_mul(jc, denom_coeffs[j]),
                                        complex_pow_int(z, j-1)));
                }
            }
            if (complex_abs(fz) < 1e-10) {
                int dup = 0;
                for (int j = 0; j < n_poles; j++)
                    if (complex_distance(z, poles[j]) < 1e-8) { dup = 1; break; }
                if (!dup) poles[n_poles++] = z;
                break;
            }
            if (complex_abs2(dfz) < 1e-30) break;
            z = complex_sub(z, complex_div(fz, dfz));
        }
    }
    double result = 0.0;
    for (int i = 0; i < n_poles; i++) {
        Complex num_val = {0,0}, den_deriv = {0,0};
        Complex z_pow = complex_pow_int(poles[i], k - 1);
        for (int j = 0; j <= deg_num; j++)
            num_val = complex_add(num_val, complex_mul(numer_coeffs[j],
                                  complex_pow_int(poles[i], j)));
        Complex weighted_num = complex_mul(num_val, z_pow);
        for (int j = 1; j <= deg_den; j++) {
            Complex jc = {(double)j, 0};
            den_deriv = complex_add(den_deriv, complex_mul(complex_mul(jc, denom_coeffs[j]),
                                    complex_pow_int(poles[i], j-1)));
        }
        Complex residue = complex_div(weighted_num, den_deriv);
        result += residue.re;
    }
    free(poles);
    return result;
}

/* Bilinear transform (Tustin): s = (2/T) * (z-1)/(z+1) */
Complex bilinear_s_to_z(Complex s, double T) {
    Complex two_over_T = {2.0/T, 0.0};
    Complex num = complex_add(s, two_over_T);
    Complex den = complex_sub(two_over_T, s);
    return complex_div(num, den);
}

Complex bilinear_z_to_s(Complex z, double T) {
    Complex two_over_T = {2.0/T, 0.0};
    Complex z_minus_1 = complex_sub(z, (Complex){1,0});
    Complex z_plus_1 = complex_add(z, (Complex){1,0});
    Complex ratio = complex_div(z_minus_1, z_plus_1);
    return complex_mul(two_over_T, ratio);
}

/* Convolution theorem verification: L{f*g} = L{f}*L{g} */
double convolution_numeric(const double *f, const double *g, int n, double t) {
    if (!f || !g || n < 2) return 0.0;
    double dt = t / (n - 1);
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        double tau = i * dt;
        int idx = (int)((t - tau) / dt);
        if (idx >= 0 && idx < n) {
            result += f[i] * g[idx] * dt;
        }
    }
    return result;
}

/* ========================================================================
 * Additional utility: compute system response using convolution
 * y(t) = h(t) * u(t) = int_0^t h(tau) u(t-tau) dtau
 * where h(t) is the impulse response.
 * ======================================================================== */

double system_response_convolution(const double *impulse_resp, int n_h,
                                    const double *input, int n_in,
                                    double t, double dt) {
    double result = 0.0;
    int max_k = (int)(t / dt);
    for (int k = 0; k < max_k && k < n_h; k++) {
        double tau = k * dt;
        int in_idx = (int)((t - tau) / dt);
        if (in_idx >= 0 && in_idx < n_in) {
            result += impulse_resp[k] * input[in_idx] * dt;
        }
    }
    return result;
}
