/**
 * @file complex_integration.h
 * @brief Complex contour integration -- Cauchy theorems, residue calculus,
 *        winding numbers, argument principle
 *
 * L2-L5: Contour integrals, residue theorem, real integral evaluation
 *
 * Course alignment:
 *   MIT 18.04 -- Residue theorem, real integrals via complex methods
 *   Princeton MAT 330 -- Contour integration
 *   ETH 401-0353 -- Residuensatz
 *
 * Textbook: Churchill & Brown; Ahlfors; Lang "Complex Analysis"
 */

#ifndef COMPLEX_INTEGRATION_H
#define COMPLEX_INTEGRATION_H

#include "complex_number.h"
#include <stddef.h>

/* Contour parametric path */
typedef Complex (*ContourParam)(double t);

/* Numeric description of a contour segment */
typedef struct {
    ContourParam param;
    double       t0;
    double       t1;
    int          closed;
} Contour;

typedef Complex (*ComplexFunc)(Complex z);


/* Pre-sampled contour data (internal representation) */
typedef struct {
    Complex *points;
    int      n_seg;
    int      closed;
} ContourData;

ContourData *contour_line_make(Complex z0, Complex z1, int n);
ContourData *contour_circle_make(Complex center, double radius, int n);
ContourData *contour_semicircle_upper(Complex center, double radius, int n);
ContourData *contour_rectangle_make(Complex v0, Complex v1, Complex v2, Complex v3, int n);
void contour_data_free(ContourData *cd);

Contour contour_line(Complex z0, Complex z1);
Contour contour_circle_arc(Complex center, double radius,
                            double theta0, double theta1);
Contour contour_circle(Complex center, double radius);
Contour contour_rectangle(Complex v0, Complex v1, Complex v2, Complex v3);
Complex contour_integral_simpson(ComplexFunc f, const Contour *gamma, int n);
Complex contour_integral_trapezoidal(ComplexFunc f, const Contour *gamma, int n);
double winding_number(const Contour *gamma, Complex z0, int n);
Complex cauchy_integral_formula(ComplexFunc f, const Contour *gamma,
                                 Complex z0, int n);
Complex cauchy_derivative_formula(ComplexFunc f, int order,
                                   const Contour *gamma, Complex z0, int n);
Complex residue_at_pole(ComplexFunc f, Complex z0, int pole_order, double h);
Complex residue_theorem_sum(ComplexFunc f, const Complex *sing_points,
                             const int *pole_orders, int n_sing, double h);
double real_integral_via_residues(const Complex *numer_coeffs, int deg_num,
                                   const Complex *denom_coeffs, int deg_den,
                                   double tol);
double argument_principle(ComplexFunc f, ComplexFunc df,
                           const Contour *gamma, int n);

#endif /* COMPLEX_INTEGRATION_H */
