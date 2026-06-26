/**
 * @file complex_transform.h
 * @brief Integral transforms in the complex plane -- Laplace, Fourier, Z
 *
 * L4-L5: Transform theory bridges time/frequency domains via complex analysis.
 *
 * Course alignment:
 *   MIT 6.003 -- Signals and Systems
 *   Stanford EE 102 -- Systems and Signals
 *   Michigan EECS 216 -- Signals and Systems
 *
 * Textbook: Oppenheim and Willsky, Signals and Systems (1996)
 */

#ifndef COMPLEX_TRANSFORM_H
#define COMPLEX_TRANSFORM_H

#include "complex_number.h"
#include <stddef.h>

/* Numeric Laplace transform of sampled data (composite Simpson) */
Complex laplace_transform_numeric(const double *t, const double *f,
                                   int n, Complex s);

/* Inverse Laplace via Bromwich integral (residue-based for rational F(s)) */
double inverse_laplace_rational(const Complex *numer_coeffs, int deg_num,
                                 const Complex *denom_coeffs, int deg_den,
                                 double t);

/* Numeric Fourier transform (DFT-based) */
void fourier_transform_numeric(const double *t, const double *f, int n,
                                double omega_min, double omega_max,
                                int n_freq, Complex *out);

/* Partial fraction expansion results */
typedef struct {
    Complex *poles;
    int     *mults;
    Complex *residues;
    int      n_poles;
    Complex *direct;
    int      direct_deg;
} PartialFraction;

PartialFraction partial_fraction_expand(const Complex *num, int deg_num,
                                         const Complex *den, int deg_den);
Complex partial_fraction_eval(const PartialFraction *pf, Complex s);
void partial_fraction_free(PartialFraction *pf);

/* Z-transform of a finite sequence */
Complex z_transform(const double *x, int n, Complex z);

/* Inverse Z-transform via residue method */
double inverse_z_transform_residue(const Complex *numer_coeffs, int deg_num,
                                    const Complex *denom_coeffs, int deg_den,
                                    int k);

#endif /* COMPLEX_TRANSFORM_H */
