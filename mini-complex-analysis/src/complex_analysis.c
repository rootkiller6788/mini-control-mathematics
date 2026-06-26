/**
 * @file complex_analysis.c
 * @brief Analytic function theory: Cauchy-Riemann, harmonic functions,
 *        conformal mappings, analytic continuation, Newton method
 *
 * L1-L5: Core definitions through numerical methods
 *
 * Textbook: Ahlfors (1979), Needham (1997), Churchill & Brown (2014)
 */

#include "complex_analysis.h"
#include "complex_number.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

Complex wirtinger_dz(ComplexFunc f, Complex z0, double h) {
    Complex z_ph = {z0.re + h, z0.im};
    Complex z_mh = {z0.re - h, z0.im};
    Complex z_pih = {z0.re, z0.im + h};
    Complex z_mih = {z0.re, z0.im - h};
    Complex fp = f(z_ph); Complex fm = f(z_mh);
    Complex fip = f(z_pih); Complex fim = f(z_mih);
    Complex dx = complex_sub(fp, fm);
    Complex dy = complex_sub(fip, fim);
    Complex i_dy = {-dy.im, dy.re};
    Complex diff = complex_sub(dx, i_dy);
    Complex result = complex_scale(diff, 0.25 / h);
    return result;
}

Complex wirtinger_dzbar(ComplexFunc f, Complex z0, double h) {
    Complex z_ph = {z0.re + h, z0.im};
    Complex z_mh = {z0.re - h, z0.im};
    Complex z_pih = {z0.re, z0.im + h};
    Complex z_mih = {z0.re, z0.im - h};
    Complex fp = f(z_ph); Complex fm = f(z_mh);
    Complex fip = f(z_pih); Complex fim = f(z_mih);
    Complex dx = complex_sub(fp, fm);
    Complex dy = complex_sub(fip, fim);
    Complex i_dy = {dy.im, -dy.re};
    Complex summ = complex_add(dx, i_dy);
    Complex result = complex_scale(summ, 0.25 / h);
    return result;
}

int cauchy_riemann_check(ComplexFunc f, Complex z0, double h, double tol) {
    Complex dzbar = wirtinger_dzbar(f, z0, h);
    return (complex_abs(dzbar) < tol);
}

double laplacian_2d(double (*u)(double x, double y), double x, double y, double h) {
    double u_xx = (u(x+h, y) - 2.0*u(x, y) + u(x-h, y)) / (h*h);
    double u_yy = (u(x, y+h) - 2.0*u(x, y) + u(x, y-h)) / (h*h);
    return u_xx + u_yy;
}

int is_harmonic(double (*u)(double x, double y), double x, double y, double h, double tol) {
    return (fabs(laplacian_2d(u, x, y, h)) < tol);
}

Complex mobius_apply(const MobiusTransform *m, Complex z) {
    Complex num1 = complex_mul(m->a, z);
    Complex num = complex_add(num1, m->b);
    Complex den1 = complex_mul(m->c, z);
    Complex den = complex_add(den1, m->d);
    return complex_div(num, den);
}

Complex mobius_det(const MobiusTransform *m) {
    Complex ad = complex_mul(m->a, m->d);
    Complex bc = complex_mul(m->b, m->c);
    return complex_sub(ad, bc);
}

MobiusTransform mobius_inverse(const MobiusTransform *m) {
    MobiusTransform inv;
    inv.a = m->d;
    inv.b = complex_neg(m->b);
    inv.c = complex_neg(m->c);
    inv.d = m->a;
    return inv;
}

MobiusTransform cayley_transform(void) {
    MobiusTransform m;
    m.a.re = 1.0; m.a.im = 0.0;
    m.b.re = -1.0; m.b.im = 0.0;
    m.c.re = 0.0; m.c.im = 1.0;
    m.d.re = 0.0; m.d.im = 1.0;
    return m;
}

double convergence_radius(const Complex *coeffs, int max_deg) {
    if (!coeffs || max_deg <= 0) return 0.0;
    double limsup = 0.0;
    for (int n = 0; n < max_deg; n++) {
        double abs_coeff = complex_abs(coeffs[n]);
        if (abs_coeff > 0.0) {
            double root = pow(abs_coeff, 1.0 / (n + 1));
            if (root > limsup) limsup = root;
        }
    }
    if (limsup == 0.0) return INFINITY;
    return 1.0 / limsup;
}


PowerSeries analytic_continuation(const PowerSeries *ps, Complex new_center, int new_deg) {
    PowerSeries result;
    result.center = new_center;
    result.max_deg = new_deg;
    result.coeffs = (Complex *)calloc(new_deg, sizeof(Complex));
    result.radius = 0.0;
    if (!ps || !ps->coeffs || !result.coeffs || ps->max_deg <= 0 || new_deg <= 0)
        return result;
    Complex shift = complex_sub(new_center, ps->center);
    for (int k = 0; k < new_deg && k < ps->max_deg; k++) {
        Complex sum = {0.0, 0.0};
        for (int n = k; n < ps->max_deg; n++) {
            double binom = 1.0;
            for (int j = 1; j <= k; j++) {
                binom *= (double)(n - k + j) / j;
            }
            Complex shift_pow = complex_pow_int(shift, n - k);
            Complex term = complex_scale(complex_mul(ps->coeffs[n], shift_pow), binom);
            sum = complex_add(sum, term);
        }
        result.coeffs[k] = sum;
    }
    result.radius = convergence_radius(result.coeffs, new_deg);
    return result;
}

Complex complex_newton(ComplexFunc f, ComplexFunc df, Complex z0, int maxit, double tol) {
    Complex z = z0;
    for (int it = 0; it < maxit; it++) {
        Complex fz = f(z);
        if (complex_abs(fz) < tol) return z;
        Complex dfz = df(z);
        double d_abs2 = complex_abs2(dfz);
        if (d_abs2 < 1e-30) break;
        Complex step = complex_div(fz, dfz);
        z = complex_sub(z, step);
    }
    return z;
}

Complex dominant_zero(const ComplexPoly *poly, Complex z0, int maxit, double tol) {
    Complex z = z0;
    for (int it = 0; it < maxit; it++) {
        Complex fz = complex_poly_eval(poly, z);
        if (complex_abs(fz) < tol) return z;
        Complex dfz = {0.0, 0.0};
        for (int k = 1; k <= poly->degree; k++) {
            Complex kc = {(double)k, 0.0};
            Complex z_pow = complex_pow_int(z, k - 1);
            Complex term = complex_mul(complex_mul(kc, poly->coeffs[k]), z_pow);
            dfz = complex_add(dfz, term);
        }
        double d_abs2 = complex_abs2(dfz);
        if (d_abs2 < 1e-30) break;
        Complex step = complex_div(fz, dfz);
        z = complex_sub(z, step);
    }
    return z;
}


Complex joukowsky_transform(Complex z) {
    Complex one = {1.0, 0.0};
    Complex inv = complex_div(one, z);
    return complex_add(z, inv);
}

Complex joukowsky_inverse(Complex w) {
    Complex half = {0.5, 0.0};
    Complex w_sq = complex_mul(w, w);
    Complex four = {4.0, 0.0};
    Complex diff = complex_sub(w_sq, four);
    Complex sqrt_term = complex_sqrt(diff);
    Complex sum = complex_add(w, sqrt_term);
    return complex_mul(half, sum);
}

Complex exp_map(Complex z) { return complex_exp(z); }
Complex log_map(Complex z) { return complex_log(z); }

Complex power_map(Complex z, double alpha) {
    Complex alpha_c = {alpha, 0.0};
    return complex_pow(z, alpha_c);
}

Complex upper_half_to_disk(Complex z) {
    Complex i = {0.0, 1.0};
    Complex num = complex_sub(z, i);
    Complex den = complex_add(z, i);
    return complex_div(num, den);
}

Complex disk_to_upper_half(Complex w) {
    Complex i = {0.0, 1.0};
    Complex one_c = {1.0, 0.0};
    Complex num = complex_add(one_c, w);
    Complex den = complex_sub(one_c, w);
    return complex_mul(i, complex_div(num, den));
}


Complex schwarz_christoffel_integrand(Complex z,
                                       const double *prevertices, int n,
                                       const double *angles) {
    Complex prod = {1.0, 0.0};
    for (int k = 0; k < n; k++) {
        Complex xk = {prevertices[k], 0.0};
        Complex diff = complex_sub(z, xk);
        double alpha = angles[k] / M_PI - 1.0;
        Complex alpha_c = {alpha, 0.0};
        Complex diff_pow = complex_pow(diff, alpha_c);
        prod = complex_mul(prod, diff_pow);
    }
    return prod;
}

int check_max_modulus_principle(ComplexFunc f, const ContourData *boundary,
                                 Complex *interior_pts, int n_interior, double tol) {
    if (!f || !boundary || !interior_pts || n_interior <= 0) return 1;
    double max_boundary = 0.0;
    for (int i = 0; i < boundary->n_seg; i++) {
        double mag = complex_abs(f(boundary->points[i]));
        if (mag > max_boundary) max_boundary = mag;
    }
    for (int i = 0; i < n_interior; i++) {
        double mag = complex_abs(f(interior_pts[i]));
        if (mag > max_boundary + tol) return 0;
    }
    return 1;
}

double liouville_bound_check(ComplexFunc f, Complex z0, double h,
                              double M, double tol) {
    Complex dz = wirtinger_dz(f, z0, h);
    double deriv_mag = complex_abs(dz);
    (void)M; (void)tol;
    return deriv_mag;
}

int rouche_check(ComplexFunc f, ComplexFunc g, const ContourData *gamma, double tol) {
    if (!f || !g || !gamma) return 1;
    for (int i = 0; i <= gamma->n_seg; i++) {
        Complex z = gamma->points[i];
        double f_mag = complex_abs(f(z));
        double g_mag = complex_abs(g(z));
        if (f_mag <= g_mag + tol) return 0;
    }
    return 1;
}


int analytic_continue_along_path(ComplexFunc f, Complex *path, int n_path,
                                  TaylorSeries *ts_out, double tol) {
    if (!f || !path || n_path < 2 || !ts_out) return 0;
    Complex z_current = path[0];
    *ts_out = taylor_from_derivatives(f, z_current, 10, 1e-4);
    for (int i = 1; i < n_path; i++) {
        Complex z_next = path[i];
        if (complex_distance(z_current, z_next) > ts_out->radius * 0.5) {
            Complex z_mid;
            z_mid.re = 0.5*(z_current.re + z_next.re);
            z_mid.im = 0.5*(z_current.im + z_next.im);
            Complex f_mid = taylor_eval(ts_out, z_mid);
            Complex f_actual = f(z_mid);
            if (complex_distance(f_mid, f_actual) > tol) return i;
            taylor_free(ts_out);
            *ts_out = taylor_from_derivatives(f, z_mid, 10, 1e-4);
            z_current = z_mid;
        } else {
            Complex f_val = taylor_eval(ts_out, z_next);
            Complex f_actual = f(z_next);
            if (complex_distance(f_val, f_actual) > tol) return i;
            z_current = z_next;
        }
    }
    return n_path;
}

int verify_schwarz_reflection(ComplexFunc f, Complex z, double tol) {
    Complex conj_z = complex_conj(z);
    Complex f_conj_z = f(conj_z);
    Complex conj_f_z = complex_conj(f(z));
    return (complex_distance(f_conj_z, conj_f_z) < tol);
}

double poisson_kernel(double r, double theta) {
    double num = 1.0 - r * r;
    double den = 1.0 - 2.0 * r * cos(theta) + r * r;
    return num / den;
}

double poisson_integral_formula(double (*boundary_u)(double phi),
                                 double r, double theta, int n_pts) {
    double sum = 0.0;
    double dphi = 2.0 * M_PI / n_pts;
    for (int i = 0; i < n_pts; i++) {
        double phi = i * dphi;
        double kernel = poisson_kernel(r, theta - phi);
        sum += kernel * boundary_u(phi) * dphi;
    }
    return sum / (2.0 * M_PI);
}
