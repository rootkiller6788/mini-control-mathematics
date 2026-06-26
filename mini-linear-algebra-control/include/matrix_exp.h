/**
 * @file matrix_exp.h
 * @brief Matrix exponential, logarithm, and related functions for control
 * @module mini-linear-algebra-control
 *
 * Knowledge coverage:
 *   L1: Matrix exponential definition: e^A = sum_{k=0}^{inf} A^k/k!
 *   L2: State transition matrix Phi(t) = e^{At}
 *   L4: Properties: e^{A+B} = e^A e^B iff AB = BA; d/dt e^{At} = A e^{At}
 *   L5: Computation methods: scaling+squaring, Pade approximation, Krylov
 *   L6: Discretization: A_d = e^{A T_s}, B_d = integral_0^{T_s} e^{A tau} B d_tau
 *
 * Reference: Moler & Van Loan, Nineteen Dubious Ways to Compute the Matrix
 *            Exponential (1978, updated 2003)
 *            Higham, Functions of Matrices: Theory and Computation (2008)
 */

#ifndef MINI_LINALG_MATRIX_EXP_H
#define MINI_LINALG_MATRIX_EXP_H

#include "matrix.h"
#include "vector.h"

/* =========================================================================
 * L5: Matrix Exponential Computation
 * ========================================================================= */

/**
 * @brief Compute matrix exponential e^A using scaling and squaring
 * with Pade approximation.
 *
 * Algorithm (L5): Scale A by 2^{-s}, compute Pade approximant R_{p,q}(A/2^s),
 * then square result s times: e^A = (R_{p,q}(A/2^s))^{2^s}.
 *
 * This is the algorithm used by MATLAB's expm (Higham, 2005).
 * Complexity: O(n^3).
 *
 * @param A      Input nxn matrix
 * @param result Output nxn matrix exponential
 * @return       0 on success
 */
int matrix_exp(const Matrix *A, Matrix *result);

/**
 * @brief Compute matrix exponential via Taylor series (naive method).
 * e^A = I + A + A^2/2! + ... + A^k/k!
 * Simple but can have cancellation issues for large ||A||.
 * (L5: Pedagogical method - illustrates the definition)
 */
int matrix_exp_taylor(const Matrix *A, Matrix *result, int terms);

/**
 * @brief Compute matrix exponential via eigendecomposition.
 * If A = V diag(lambda_i) V^{-1}, then e^A = V diag(e^{lambda_i}) V^{-1}.
 * Only works for diagonalizable matrices.
 * (L5: Spectral method - exact when eigendecomposition is available)
 */
int matrix_exp_eigen(const Matrix *A, Matrix *result);

/* =========================================================================
 * L6: State Transition and Discretization
 * ========================================================================= */

/**
 * @brief Compute state transition matrix Phi(t) = e^{A t}.
 * (L6: Solution of dx/dt = Ax is x(t) = Phi(t) * x(0))
 */
int matrix_exp_t(const Matrix *A, double t, Matrix *Phi);

/**
 * @brief Discretize continuous-time state-space system.
 * Given dx/dt = Ax + Bu, with sampling period T_s:
 *   x_{k+1} = A_d x_k + B_d u_k
 * where A_d = e^{A T_s} and B_d = (A_d - I) A^{-1} B.
 *
 * Uses matrix exponential for A_d and numerical integration for B_d.
 * (L6: Zero-order-hold discretization, fundamental for digital control)
 */
int discretize_zoh(const Matrix *A, const Matrix *B, double Ts,
                    Matrix *Ad, Matrix *Bd);

/**
 * @brief Compute matrix logarithm: inverse of matrix exponential.
 * A = logm(X) satisfies e^A = X.
 * Uses inverse scaling-and-squaring + Pade approximation.
 * (L8: Matrix logarithm for continuous-time system identification)
 */
int matrix_log(const Matrix *X, Matrix *A);

/* =========================================================================
 * L8: Matrix Functions
 * ========================================================================= */

/**
 * @brief Compute matrix cosine: cos(A) = (e^{iA} + e^{-iA}) / 2.
 * (L8: Arises in second-order systems, wave equations)
 */
int matrix_cos(const Matrix *A, Matrix *result);

/**
 * @brief Compute matrix sine: sin(A) = (e^{iA} - e^{-iA}) / (2i).
 */
int matrix_sin(const Matrix *A, Matrix *result);

/**
 * @brief Compute A^t for real t (matrix power with real exponent).
 * A^t = e^{t * log(A)} for A having no eigenvalues on negative real axis.
 * (L8: Fractional powers in fractional-order control systems, L9)
 */
int matrix_pow_real(const Matrix *A, double t, Matrix *result);

/**
 * @brief Frechet derivative of matrix exponential: L(A, E) = d/dt e^{A+tE}|_{t=0}
 * Used for sensitivity analysis of state transition.
 * (L8: Perturbation analysis of matrix functions)
 */
int matrix_exp_frechet(const Matrix *A, const Matrix *E, Matrix *L);

#endif /* MINI_LINALG_MATRIX_EXP_H */