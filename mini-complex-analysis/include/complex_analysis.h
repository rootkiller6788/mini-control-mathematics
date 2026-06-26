/**
 * @file complex_analysis.h
 * @brief Analytic function theory — Cauchy-Riemann, harmonic functions,
 *        conformal mapping, singularities, analytic continuation
 *
 * L1-L4: Core definitions through fundamental theorems of complex analysis
 *
 * Course alignment:
 *   MIT 18.04 — Cauchy-Riemann, harmonic conjugates, analytic functions
 *   Stanford MATH 106 — analyticity criteria
 *   Berkeley MATH 185 — holomorphic function theory
 *   Cambridge Part II — Complex Methods
 *   ETH 401-0353 — Funktionentheorie
 *
 * Textbook: Ahlfors, "Complex Analysis" (1979); Needham, "Visual Complex Analysis" (1997)
 */

#ifndef COMPLEX_ANALYSIS_H
#define COMPLEX_ANALYSIS_H

#include "complex_number.h"
#include "complex_series.h"
#include "complex_integration.h"
#include <stddef.h>

/* --------------------------------------------------------------------------
 * L1: Derivative in complex sense — Wirtinger derivatives ∂/∂z and ∂/∂z̄
 *
 * ∂f/∂z = ½(∂f/∂x - i·∂f/∂y)
 * ∂f/∂z̄ = ½(∂f/∂x + i·∂f/∂y)
 *
 * f is analytic ⟺ ∂f/∂z̄ = 0 (Cauchy-Riemann condition)
 * -------------------------------------------------------------------------- */

/** Function pointer type for complex-valued functions of a complex variable */
typedef Complex (*ComplexFunc)(Complex z);

/**
 * @brief  Approximate ∂f/∂z at z0 using central difference (h real + i·h)
 * @param  f   function to differentiate
 * @param  z0  point of evaluation
 * @param  h   step size
 * @return  ∂f/∂z approximation
 */
Complex wirtinger_dz(ComplexFunc f, Complex z0, double h);

/**
 * @brief  Approximate ∂f/∂z̄ at z0 (should be ~0 for analytic f)
 * @return  ∂f/∂z̄ approximation
 */
Complex wirtinger_dzbar(ComplexFunc f, Complex z0, double h);

/* --------------------------------------------------------------------------
 * L2: Cauchy-Riemann equations — analyticity test
 *
 * u_x = v_y,  u_y = -v_x   where f(z) = u(x,y) + i·v(x,y)
 *
 * CR check reads real and imaginary parts of Wirtinger dzbar.
 * -------------------------------------------------------------------------- */

/**
 * @brief  Check Cauchy-Riemann at z0 numerically
 * @return 1 if CR holds within tol, 0 otherwise
 */
int cauchy_riemann_check(ComplexFunc f, Complex z0, double h, double tol);

/* --------------------------------------------------------------------------
 * L2: Harmonic functions and harmonic conjugates
 *
 * If f = u+iv is analytic, then u and v are harmonic: ∇²u = ∇²v = 0.
 * Given u, find harmonic conjugate v (up to constant).
 *
 * Laplace equation: u_xx + u_yy = 0 in 2D.
 * -------------------------------------------------------------------------- */

/**
 * @brief  Compute Laplacian ∇²u(x,y) numerically (5-point stencil)
 * @param  u   real function of two real variables
 * @param  x   x-coordinate
 * @param  y   y-coordinate
 * @param  h   grid spacing
 */
double laplacian_2d(double (*u)(double x, double y), double x, double y, double h);

/**
 * @brief  Check if u is harmonic at (x,y): |∇²u| < tol
 */
int is_harmonic(double (*u)(double x, double y), double x, double y, double h, double tol);

/* --------------------------------------------------------------------------
 * L2: Conformal mapping — angle-preserving transformations
 *
 * Möbius transform: T(z) = (a·z + b) / (c·z + d), ad - bc ≠ 0
 *
 * Conformal maps preserve angles locally. They are used to transform
 * complex boundary value problems into simpler domains (e.g., unit disk).
 * -------------------------------------------------------------------------- */

/** Parameters for Möbius transformation T(z) = (a·z + b) / (c·z + d) */
typedef struct {
    Complex a, b, c, d;
} MobiusTransform;

/** Apply Möbius transform: T(z) = (a·z + b)/(c·z + d) */
Complex mobius_apply(const MobiusTransform *m, Complex z);

/** Compute determinant ad - bc (non-zero ⟺ invertible) */
Complex mobius_det(const MobiusTransform *m);

/** Compute inverse Möbius transform parameters */
MobiusTransform mobius_inverse(const MobiusTransform *m);

/** Möbius mapping the upper half-plane to unit disk (Cayley transform) */
MobiusTransform cayley_transform(void);

/* --------------------------------------------------------------------------
 * L3: Analytic continuation — extending domain of analytic functions
 *
 * Given power series coefficients at z0 with radius r,
 * re-expand at z1 inside convergence disk.
 * -------------------------------------------------------------------------- */

/** Power series: Σ a_n (z - center)^n */
typedef struct {
    Complex  center;     /**< expansion center */
    int      max_deg;    /**< number of terms */
    Complex *coeffs;     /**< coeffs[0..max_deg-1] */
    double   radius;     /**< radius of convergence (estimated) */
} PowerSeries;

/**
 * @brief  Estimate radius of convergence via Cauchy-Hadamard:
 *         R = 1 / limsup |a_n|^{1/n}
 */
double convergence_radius(const Complex *coeffs, int max_deg);

/**
 * @brief  Re-expand power series from center0 to center1
 *         (binomial expansion of (z-c1) terms)
 */
PowerSeries analytic_continuation(const PowerSeries *ps, Complex new_center, int new_deg);

/* --------------------------------------------------------------------------
 * L5: Numerical methods for analytic functions
 * -------------------------------------------------------------------------- */

/**
 * @brief  Find zeros of f(z) using Newton-Raphson in complex plane
 * @param  f     analytic function
 * @param  df    derivative of f
 * @param  z0    initial guess
 * @param  maxit max iterations
 * @param  tol   convergence tolerance |f(z)| < tol
 * @return root approximation
 */
Complex complex_newton(ComplexFunc f, ComplexFunc df, Complex z0, int maxit, double tol);

/**
 * @brief  Duration-Power method: given f(z) polynomial, find dominant zero
 *         z_{k+1} = z_k - f(z_k)/f'(z_k) with deflation
 */
Complex dominant_zero(const ComplexPoly *poly, Complex z0, int maxit, double tol);


/* Additional conformal mappings */
Complex joukowsky_transform(Complex z);
Complex joukowsky_inverse(Complex w);
Complex exp_map(Complex z);
Complex log_map(Complex z);
Complex power_map(Complex z, double alpha);
Complex upper_half_to_disk(Complex z);
Complex disk_to_upper_half(Complex w);

/* Schwarz-Christoffel mapping */
Complex schwarz_christoffel_integrand(Complex z, const double *prevertices, int n, const double *angles);

/* Theorem verifications */
int check_max_modulus_principle(ComplexFunc f, const ContourData *boundary, Complex *interior_pts, int n, double tol);
double liouville_bound_check(ComplexFunc f, Complex z0, double h, double M, double tol);
int rouche_check(ComplexFunc f, ComplexFunc g, const ContourData *gamma, double tol);
int analytic_continue_along_path(ComplexFunc f, Complex *path, int n, TaylorSeries *ts, double tol);
int verify_schwarz_reflection(ComplexFunc f, Complex z, double tol);

/* Poisson integral formula */
double poisson_kernel(double r, double theta);
double poisson_integral_formula(double (*u)(double phi), double r, double theta, int n);

#endif /* COMPLEX_ANALYSIS_H */
