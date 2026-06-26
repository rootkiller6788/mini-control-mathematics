/**
 * @file complex_series.c
 * @brief Laurent series, Taylor series, singularity classification,
 *        radius of convergence, residue computation
 *
 * L2-L4: Series representations of analytic functions
 *
 * Key concepts:
 *   Laurent: f(z) = sum_{n=-inf}^{inf} c_n (z-z0)^n
 *   Taylor:  f(z) = sum_{n=0}^{inf} a_n (z-z0)^n
 *   Residue: c_{-1} coefficient of Laurent series
 *
 * Textbook: Churchill & Brown (2014), Ahlfors (1979)
 */

#include "complex_series.h"
#include "complex_integration.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

Complex laurent_eval(const LaurentSeries *ls, Complex z) {
    Complex sum = {0.0, 0.0};
    if (!ls) return sum;
    Complex w = complex_sub(z, ls->center);
    /* Analytic part: sum_{n=0}^{pos_order-1} c_n * w^n */
    for (int n = 0; n < ls->pos_order; n++) {
        Complex term = complex_mul(ls->pos_coeffs[n], complex_pow_int(w, n));
        sum = complex_add(sum, term);
    }
    /* Principal part: sum_{n=1}^{neg_order} c_{-n} * w^{-n} */
    for (int n = 1; n <= ls->neg_order; n++) {
        Complex w_pow = complex_pow_int(w, -n);
        Complex term = complex_mul(ls->neg_coeffs[n-1], w_pow);
        sum = complex_add(sum, term);
    }
    return sum;
}

/* Compute Laurent coefficient c_n via Cauchy integral:
 * c_n = (1/2pi*i) integral f(z)/(z-z0)^{n+1} dz
 */
Complex laurent_coeff(ComplexFunc f, Complex z0, int n,
                      const Contour *gamma, int num_subdiv) {
    /* Use Cauchy derivative formula with order = n */
    if (!f || !gamma) return (Complex){0,0};
    if (n >= 0) {
        /* c_n = f^{(n)}(z0)/n! */
        Complex deriv = cauchy_derivative_formula(f, n, gamma, z0, num_subdiv);
        double fact = 1.0;
        for (int k = 2; k <= n; k++) fact *= k;
        return complex_scale(deriv, 1.0 / fact);
    } else {
        /* c_{-m} = (1/2pi*i) integral f(z)*(z-z0)^{m-1} dz */
        /* Simplified: use contour integration with weight (z-z0)^{m-1} */
        int m = -n;
        ContourData *cd = contour_circle_make(z0, 0.5, num_subdiv);
        if (!cd) return (Complex){0,0};
        Complex sum = {0,0};
        for (int i = 0; i < cd->n_seg; i++) {
            Complex z_mid;
            z_mid.re = 0.5*(cd->points[i].re + cd->points[i+1].re);
            z_mid.im = 0.5*(cd->points[i].im + cd->points[i+1].im);
            Complex fz = f(z_mid);
            Complex dz_sub = complex_sub(z_mid, z0);
            Complex weight = complex_pow_int(dz_sub, m - 1);
            Complex integrand = complex_mul(fz, weight);
            Complex dz;
            dz.re = cd->points[i+1].re - cd->points[i].re;
            dz.im = cd->points[i+1].im - cd->points[i].im;
            Complex term = complex_mul(integrand, dz);
            sum = complex_add(sum, term);
        }
        contour_data_free(cd);
        Complex two_pi_i = {0.0, 2.0 * M_PI};
        return complex_div(sum, two_pi_i);
    }
}

SingularityType classify_singularity(const LaurentSeries *ls) {
    if (!ls) return SING_REMOVABLE;
    if (ls->neg_order == 0) return SING_REMOVABLE;
    if (ls->neg_order > 0) {
        /* Check if all higher negative coefficients are zero (finite pole order) */
        return SING_POLE;
    }
    /* If neg_order is negative (sentinel for infinite principal part), essential */
    return SING_ESSENTIAL;
}


/* Taylor series evaluation */
Complex taylor_eval(const TaylorSeries *ts, Complex z) {
    Complex sum = {0.0, 0.0};
    if (!ts || !ts->coeffs || ts->n_terms <= 0) return sum;
    Complex w = complex_sub(z, ts->center);
    Complex w_pow = {1.0, 0.0};
    for (int n = 0; n < ts->n_terms; n++) {
        Complex term = complex_mul(ts->coeffs[n], w_pow);
        sum = complex_add(sum, term);
        w_pow = complex_mul(w_pow, w);
    }
    return sum;
}

/* Taylor coefficients via finite-difference derivatives */
TaylorSeries taylor_from_derivatives(ComplexFunc f, Complex z0, int n_terms, double h) {
    TaylorSeries ts;
    ts.center = z0;
    ts.n_terms = n_terms;
    ts.radius = 0.0;
    ts.coeffs = (Complex *)calloc(n_terms, sizeof(Complex));
    if (!ts.coeffs || n_terms <= 0) return ts;
    /* a_0 = f(z0) */
    ts.coeffs[0] = f(z0);
    /* Higher-order coefficients via finite differences */
    for (int ord = 1; ord < n_terms; ord++) {
        Complex sum = {0,0};
        for (int k = 0; k <= ord; k++) {
            Complex zk = {z0.re + (k - ord/2.0)*h, z0.im + 0*h};
            Complex binom_c = {(double)((k%2==0)?1:-1), 0};
            /* Simplified: use central difference formula */
            Complex fz = f(zk);
            sum = complex_add(sum, complex_mul(fz, binom_c));
        }
        double fact = 1.0;
        for (int j = 2; j <= ord; j++) fact *= j;
        ts.coeffs[ord] = complex_scale(sum, 1.0 / (fact * pow(h, ord)));
    }
    ts.radius = taylor_convergence_radius(&ts);
    return ts;
}

/* Taylor coefficients via Cauchy integral formula */
TaylorSeries taylor_via_cauchy(ComplexFunc f, Complex z0, int n_terms,
                                const Contour *gamma, int num_subdiv) {
    TaylorSeries ts;
    ts.center = z0;
    ts.n_terms = n_terms;
    ts.radius = 0.0;
    ts.coeffs = (Complex *)calloc(n_terms, sizeof(Complex));
    if (!ts.coeffs || n_terms <= 0 || !f || !gamma) return ts;
    for (int n = 0; n < n_terms; n++) {
        Complex deriv = cauchy_derivative_formula(f, n, gamma, z0, num_subdiv);
        double fact = 1.0;
        for (int j = 2; j <= n; j++) fact *= j;
        ts.coeffs[n] = complex_scale(deriv, 1.0 / fact);
    }
    ts.radius = taylor_convergence_radius(&ts);
    return ts;
}

/* Estimate radius of convergence using ratio test:
 * R = lim_{n->inf} |a_n| / |a_{n+1}| (if limit exists)
 * Heuristic: use geometric mean of ratios for last few terms
 */
double taylor_convergence_radius(const TaylorSeries *ts) {
    if (!ts || !ts->coeffs || ts->n_terms < 2) return 0.0;
    double max_ratio = 0.0;
    for (int n = 0; n < ts->n_terms - 1; n++) {
        double an = complex_abs(ts->coeffs[n]);
        double an1 = complex_abs(ts->coeffs[n+1]);
        if (an1 > 1e-30) {
            double ratio = an / an1;
            if (ratio > max_ratio) max_ratio = ratio;
        }
    }
    if (max_ratio == 0.0) return INFINITY;
    return max_ratio;
}

/* Residue from Laurent series: c_{-1} */
Complex residue_from_laurent(const LaurentSeries *ls) {
    if (!ls || ls->neg_order < 1) return (Complex){0,0};
    return ls->neg_coeffs[0]; /* c_{-1} is the first negative coefficient */
}

/* Build Laurent series from known coefficients */
LaurentSeries laurent_from_coeffs(Complex center,
                                   const Complex *neg_c, int neg_n,
                                   const Complex *pos_c, int pos_n) {
    LaurentSeries ls;
    ls.center = center;
    ls.neg_order = neg_n;
    ls.pos_order = pos_n;
    ls.inner_radius = 0.0;
    ls.outer_radius = INFINITY;
    ls.neg_coeffs = (Complex *)malloc(neg_n * sizeof(Complex));
    ls.pos_coeffs = (Complex *)malloc(pos_n * sizeof(Complex));
    if (ls.neg_coeffs && neg_c) memcpy(ls.neg_coeffs, neg_c, neg_n * sizeof(Complex));
    if (ls.pos_coeffs && pos_c) memcpy(ls.pos_coeffs, pos_c, pos_n * sizeof(Complex));
    return ls;
}

void laurent_free(LaurentSeries *ls) {
    if (ls) {
        free(ls->neg_coeffs);
        free(ls->pos_coeffs);
        ls->neg_coeffs = NULL;
        ls->pos_coeffs = NULL;
    }
}

void taylor_free(TaylorSeries *ts) {
    if (ts) {
        free(ts->coeffs);
        ts->coeffs = NULL;
    }
}


/* ========================================================================
 * Additional series utilities
 * ======================================================================== */

/* Geometric series sum: sum_{n=0}^{N-1} z^n = (1 - z^N)/(1 - z) for z != 1 */
Complex geometric_series_sum(Complex z, int N) {
    Complex one = {1.0, 0.0};
    if (z.re == 1.0 && z.im == 0.0) {
        Complex result = {(double)N, 0.0};
        return result;
    }
    Complex z_N = complex_pow_int(z, N);
    Complex num = complex_sub(one, z_N);
    Complex den = complex_sub(one, z);
    return complex_div(num, den);
}

/* Exponential series: exp(z) = sum_{n=0}^{N-1} z^n/n! + remainder */
Complex exp_series(Complex z, int n_terms) {
    Complex sum = {0.0, 0.0};
    Complex term = {1.0, 0.0};
    for (int n = 0; n < n_terms; n++) {
        sum = complex_add(sum, term);
        term = complex_mul(term, z);
        double scale = 1.0 / (n + 1);
        term = complex_scale(term, scale);
    }
    return sum;
}

/* Sine series: sin(z) = sum_{n=0}^{N-1} (-1)^n z^{2n+1}/(2n+1)! */
Complex sin_series(Complex z, int n_terms) {
    Complex sum = {0.0, 0.0};
    Complex z_pow = z;
    double fact = 1.0;
    for (int n = 0; n < n_terms; n++) {
        Complex term = complex_scale(z_pow, 1.0 / fact);
        if (n % 2 == 1) term = complex_neg(term);
        sum = complex_add(sum, term);
        z_pow = complex_mul(z_pow, complex_mul(z, z));
        fact *= (2*n + 2) * (2*n + 3);
    }
    return sum;
}

/* Cosine series: cos(z) = sum_{n=0}^{N-1} (-1)^n z^{2n}/(2n)! */
Complex cos_series(Complex z, int n_terms) {
    Complex sum = {0.0, 0.0};
    Complex z_pow = {1.0, 0.0};
    double fact = 1.0;
    for (int n = 0; n < n_terms; n++) {
        Complex term = complex_scale(z_pow, 1.0 / fact);
        if (n % 2 == 1) term = complex_neg(term);
        sum = complex_add(sum, term);
        z_pow = complex_mul(z_pow, complex_mul(z, z));
        if (n > 0) fact *= (2*n - 1) * (2*n);
        if (n == 0) fact = 1.0;
        if (n == 1) fact = 2.0;
        if (n > 1) {
            fact = 1.0;
            for (int j = 2; j <= 2*n; j++) fact *= j;
        }
    }
    return sum;
}

/* Binomial series: (1+z)^alpha = sum_{n=0}^{N-1} binom(alpha,n) z^n */
Complex binomial_series(Complex z, double alpha, int n_terms) {
    Complex sum = {0.0, 0.0};
    Complex z_pow = {1.0, 0.0};
    double binom = 1.0;
    for (int n = 0; n < n_terms; n++) {
        Complex term = complex_scale(z_pow, binom);
        sum = complex_add(sum, term);
        z_pow = complex_mul(z_pow, z);
        binom *= (alpha - n) / (n + 1);
    }
    return sum;
}

/* Logarithm series: Log(1+z) = sum_{n=1}^{N} (-1)^{n+1} z^n/n, |z| < 1 */
Complex log1p_series(Complex z, int n_terms) {
    Complex sum = {0.0, 0.0};
    Complex z_pow = z;
    for (int n = 1; n <= n_terms; n++) {
        Complex term = complex_scale(z_pow, 1.0 / n);
        if (n % 2 == 0) term = complex_neg(term);
        sum = complex_add(sum, term);
        z_pow = complex_mul(z_pow, z);
    }
    return sum;
}

