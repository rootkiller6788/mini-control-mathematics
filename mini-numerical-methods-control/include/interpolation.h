/**
 * interpolation.h — Interpolation and Approximation Methods
 *
 * Covers: L5 Methods (interpolation), L6 Problems (data fitting)
 *
 * Interpolation is essential in control for: signal reconstruction from
 * sampled data, lookup tables in gain-scheduled control, trajectory
 * generation for robotics, and numerical differentiation.
 *
 * Reference: de Boor, A Practical Guide to Splines (1978)
 *            Powell, Approximation Theory and Methods (1981)
 *            Quarteroni, Sacco, Saleri, Numerical Mathematics (2007)
 */

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "numerical_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Interpolation Data Types
 * ========================================================================== */

/**
 * InterpData — 1D interpolation data points
 *
 * L1 Definition: A set of ordered pairs (x_i, y_i) for i = 0,...,n-1
 * with x_0 < x_1 < ... < x_{n-1} (strictly increasing).
 */
typedef struct {
    double *x;      /**< Abscissae (independent variable), length n */
    double *y;      /**< Ordinates (dependent variable), length n */
    size_t  n;      /**< Number of data points */
} InterpData;

/**
 * Spline — cubic spline representation
 *
 * L1 Definition: A piecewise cubic polynomial s(x) that is C² continuous.
 * Natural spline: s''(x₀) = s''(x_n) = 0.
 * Clamped spline: s'(x₀) = f'(x₀), s'(x_n) = f'(x_n).
 */
typedef struct {
    InterpData data;      /**< Knot points */
    double    *a;         /**< Constant coefficients a_i = y_i */
    double    *b;         /**< Linear coefficients */
    double    *c;         /**< Quadratic coefficients (spline curvature) */
    double    *d;         /**< Cubic coefficients */
    size_t     n_segments;/**< Number of polynomial segments = n-1 */
} Spline;

/**
 * InterpData2D — 2D gridded interpolation data
 *
 * L1 Definition: Z = f(X, Y) on a rectangular grid.
 * Used for 2D lookup tables in gain-scheduled controllers.
 */
typedef struct {
    double *x;      /**< X-grid points, length nx */
    double *y;      /**< Y-grid points, length ny */
    double *z;      /**< Function values, row-major: z[i*ny + j] = f(x_i, y_j) */
    size_t  nx;     /**< Number of X-grid points */
    size_t  ny;     /**< Number of Y-grid points */
} InterpData2D;

/* ==========================================================================
 * L5: Polynomial Interpolation
 * ========================================================================== */

/**
 * interp_linear — piecewise linear interpolation
 *
 * L5 Method: O(log n) via binary search. Continuous but C⁰ only.
 * Simplest and most robust interpolation. Used in real-time control
 * for efficiency.
 */
double interp_linear(const InterpData *data, double x);

/**
 * interp_lagrange — Lagrange polynomial interpolation
 *
 * L5 Method: p(x) = Σ y_i · ℓ_i(x) where ℓ_i(x) = Π_{j≠i} (x - x_j)/(x_i - x_j).
 * O(n²) for n points. Prone to Runge phenomenon for high-degree polynomials.
 * Best for ≤ 5-7 points unless using Chebyshev nodes.
 */
double interp_lagrange(const InterpData *data, double x);

/**
 * interp_newton_divided_diff — Newton's divided difference interpolation
 *
 * L5 Method: Builds interpolant via Newton form:
 * p(x) = f[x₀] + f[x₀,x₁](x-x₀) + f[x₀,x₁,x₂](x-x₀)(x-x₁) + ...
 * O(n²) to build divided difference table, O(n) to evaluate.
 * Convenient for adding new points incrementally.
 */
double interp_newton_divided_diff(const InterpData *data, double x);

/**
 * interp_barycentric_lagrange — numerically stable barycentric Lagrange
 *
 * L5 Method: Barycentric formula avoids overflow issues of standard Lagrange.
 * p(x) = Σ (w_i/(x-x_i))·y_i / Σ (w_i/(x-x_i)).
 * O(n) evaluation. Stable even for large n with Chebyshev nodes.
 */
double interp_barycentric_lagrange(const InterpData *data, double x);

/* ==========================================================================
 * L5/L8: Spline Interpolation
 * ========================================================================== */

/**
 * spline_cubic_natural — construct natural cubic spline
 *
 * L5 Method: Solves tridiagonal system for second derivatives.
 * Natural BC: s''(x₀) = s''(x_n) = 0. O(n).
 * Provides C² smooth interpolation. Standard for trajectory generation
 * in robotics and CNC control.
 */
Spline* spline_cubic_natural(const InterpData *data);

/**
 * spline_cubic_clamped — construct clamped cubic spline
 *
 * L5 Method: Clamped BC: s'(x₀) = deriv_left, s'(x_n) = deriv_right.
 * More accurate near boundaries when derivatives are known.
 */
Spline* spline_cubic_clamped(const InterpData *data,
                              double deriv_left, double deriv_right);

/**
 * spline_evaluate — evaluate spline at point x
 *
 * L5 Method: O(log n_segments) via binary search for containing interval,
 * then O(1) Horner evaluation of cubic polynomial on that segment.
 */
double spline_evaluate(const Spline *s, double x);

/**
 * spline_derivative — evaluate spline first derivative at x
 */
double spline_derivative(const Spline *s, double x);

/**
 * spline_second_derivative — evaluate spline second derivative at x
 */
double spline_second_derivative(const Spline *s, double x);

/**
 * spline_integral — ∫_{x0}^{xn} s(x) dx
 *
 * L5 Method: Exact integration of piecewise cubic polynomial.
 * O(n). Used for path length and area computations in control.
 */
double spline_integral(const Spline *s);

/**
 * spline_free — deallocate spline
 */
void spline_free(Spline *s);

/* ==========================================================================
 * L6: 2D Interpolation (Lookup Tables for Gain Scheduling)
 * ========================================================================== */

/**
 * interp_bilinear — bilinear interpolation on rectangular grid
 *
 * L6 Method: O(log nx + log ny). Interpolates z = f(x,y) at (x_q, y_q):
 * First interpolate in x at y_j and y_{j+1}, then interpolate in y.
 * Used for 2D lookup tables in gain-scheduled controllers.
 */
double interp_bilinear(const InterpData2D *data, double x, double y);

/**
 * interpdata_free — free 1D interpolation data
 */
void interpdata_free(InterpData *data);

/**
 * interpdata2d_free — free 2D interpolation data
 */
void interpdata2d_free(InterpData2D *data);

#ifdef __cplusplus
}
#endif

#endif /* INTERPOLATION_H */
