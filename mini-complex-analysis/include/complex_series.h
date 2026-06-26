/**
 * @file complex_series.h
 * @brief Power series, Laurent series, Taylor expansion,
 *        singularity classification, radius of convergence
 *
 * L2-L4: Series representations of analytic functions
 *
 * Course alignment:
 *   MIT 18.04 -- Power series, Laurent expansion, singularities
 *   Stanford MATH 106 -- Series methods
 *   Berkeley MATH 185 -- Laurent series and classification
 *
 * Textbook: Churchill and Brown; Ahlfors Complex Analysis
 */

#ifndef COMPLEX_SERIES_H
#define COMPLEX_SERIES_H

#include "complex_number.h"
#include "complex_integration.h"
#include <stddef.h>

/* Contour and ComplexFunc types from complex_integration.h */

/* --------------------------------------------------------------------------
 * L1: Singularity classification -- removable, pole, essential, branch
 * -------------------------------------------------------------------------- */

typedef enum {
    SING_REMOVABLE = 0,
    SING_POLE      = 1,
    SING_ESSENTIAL = 2,
    SING_BRANCH    = 3
} SingularityType;

/* --------------------------------------------------------------------------
 * L2: Laurent series -- f(z) = sum_{n=-inf}^{inf} c_n (z-z0)^n
 * -------------------------------------------------------------------------- */

typedef struct {
    Complex  center;
    Complex *neg_coeffs;
    int      neg_order;
    Complex *pos_coeffs;
    int      pos_order;
    double   inner_radius;
    double   outer_radius;
} LaurentSeries;

Complex laurent_eval(const LaurentSeries *ls, Complex z);
Complex laurent_coeff(ComplexFunc f, Complex z0, int n,
                      const Contour *gamma, int num_subdiv);
SingularityType classify_singularity(const LaurentSeries *ls);

/* --------------------------------------------------------------------------
 * L3: Taylor series -- f(z) = sum_{n=0}^{inf} a_n (z-z0)^n
 * -------------------------------------------------------------------------- */

typedef struct {
    Complex  center;
    Complex *coeffs;
    int      n_terms;
    double   radius;
} TaylorSeries;

Complex taylor_eval(const TaylorSeries *ts, Complex z);
TaylorSeries taylor_from_derivatives(ComplexFunc f, Complex z0,
                                      int n_terms, double h);
TaylorSeries taylor_via_cauchy(ComplexFunc f, Complex z0, int n_terms,
                                const Contour *gamma, int num_subdiv);
double taylor_convergence_radius(const TaylorSeries *ts);

/* --------------------------------------------------------------------------
 * L4: Residue computation from Laurent series
 * -------------------------------------------------------------------------- */

Complex residue_from_laurent(const LaurentSeries *ls);


/* Additional series utilities */
LaurentSeries laurent_from_coeffs(Complex center, const Complex *neg_c, int neg_n, const Complex *pos_c, int pos_n);
void laurent_free(LaurentSeries *ls);
void taylor_free(TaylorSeries *ts);
Complex geometric_series_sum(Complex z, int N);
Complex exp_series(Complex z, int n_terms);
Complex sin_series(Complex z, int n_terms);
Complex cos_series(Complex z, int n_terms);
Complex binomial_series(Complex z, double alpha, int n_terms);
Complex log1p_series(Complex z, int n_terms);

#endif /* COMPLEX_SERIES_H */
