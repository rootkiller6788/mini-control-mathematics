/**
 * @file control_linalg.h
 * @brief Control-theoretic linear algebra: controllability, observability,
 *        Lyapunov/Riccati equations, pole placement, Kalman decomposition
 * @module mini-linear-algebra-control
 *
 * Knowledge coverage:
 *   L1: Controllability matrix, observability matrix, Gramians
 *   L2: Controllability/observability concepts, duality, Kalman decomposition
 *   L3: Controllability/Observability indices, condition numbers
 *   L4: PBH test, Lyapunov stability theorem, separation principle
 *   L5: Lyapunov equation solvers, Riccati solvers, pole placement algorithms
 *   L6: Stabilizability, detectability, minimal realization
 *
 * Reference: Kailath, Linear Systems (1980)
 *            Zhou, Doyle & Glover, Robust and Optimal Control (1996)
 *            Ogata, Modern Control Engineering (5th Ed, 2010)
 *            Antsaklis & Michel, Linear Systems (2006)
 */

#ifndef MINI_LINALG_CONTROL_LINALG_H
#define MINI_LINALG_CONTROL_LINALG_H

#include "matrix.h"
#include "vector.h"

/* =========================================================================
 * L1/L6: Controllability and Observability Matrices
 * ========================================================================= */

/**
 * @brief Build controllability matrix C = [B AB A^2B ... A^{n-1}B].
 * For system dx/dt = Ax + Bu, the pair (A,B) is controllable iff
 * rank(C) = n (L1 definition, L6 canonical problem).
 *
 * @param A  nxn system matrix
 * @param B  nxp input matrix
 * @param C  Output nx(np) controllability matrix
 * @return   0 on success
 */
int control_controllability_matrix(const Matrix *A, const Matrix *B, Matrix *C);

/**
 * @brief Build observability matrix O = [C^T A^TC^T ... (A^{n-1})^TC^T]^T.
 * For system dx/dt = Ax + Bu, y = Cx, the pair (A,C) is observable iff
 * rank(O) = n.
 */
int control_observability_matrix(const Matrix *A, const Matrix *C, Matrix *O);

/**
 * @brief Test controllability of (A,B) via rank check.
 * @return 1 if controllable, 0 if not, -1 on error.
 */
int control_is_controllable(const Matrix *A, const Matrix *B, double tol);

/**
 * @brief Test observability of (A,C) via rank check.
 * @return 1 if observable, 0 if not, -1 on error.
 */
int control_is_observable(const Matrix *A, const Matrix *C, double tol);

/**
 * @brief Test stabilizability: all uncontrollable eigenvalues are stable.
 * (L6: Decomposes system via Kalman form, checks unstabilizable modes)
 */
int control_is_stabilizable(const Matrix *A, const Matrix *B, double tol);

/**
 * @brief Test detectability: all unobservable eigenvalues are stable.
 */
int control_is_detectable(const Matrix *A, const Matrix *C, double tol);

/* =========================================================================
 * L4/L6: PBH Test (Popov-Belevitch-Hautus)
 * ========================================================================= */

/**
 * @brief PBH controllability test: (A,B) is controllable iff
 * rank([lambda*I - A,  B]) = n for all eigenvalues lambda of A.
 *
 * Theorem (L4): The PBH rank test is necessary and sufficient for
 * controllability. This provides an eigenvalue-based characterization.
 *
 * @param A      nxn system matrix
 * @param B      nxp input matrix
 * @param lambda Eigenvalue to test
 * @return       1 if rank condition holds at lambda, 0 otherwise
 */
int control_pbh_controllability_test(const Matrix *A, const Matrix *B,
                                      double lambda_real, double lambda_imag,
                                      double tol);

/**
 * @brief PBH observability test at eigenvalue lambda.
 */
int control_pbh_observability_test(const Matrix *A, const Matrix *C,
                                    double lambda_real, double lambda_imag,
                                    double tol);

/* =========================================================================
 * L5/L6: Kalman Decomposition
 * ========================================================================= */

/**
 * @brief Kalman canonical decomposition.
 * Transforms system to reveal controllable/observable subsystems.
 * Uses orthogonal transformation to produce:
 *   [Aco A12  A13  A14 ]    [Bco]
 *   [0   Aco_ 0    A24 ]    [Bco_]
 *   [0   0    Aco~ A34 ]  B=[0   ]
 *   [0   0    0    Aco~~]    [0   ]
 *
 * where co=controllable+observable, co_=controllable+unobservable,
 * co~=uncontrollable+observable, co~~=uncontrollable+unobservable.
 *
 * (L6: Canonical problem - structural decomposition)
 */
int control_kalman_decomposition(const Matrix *A, const Matrix *B,
                                  const Matrix *C,
                                  Matrix *A_k, Matrix *B_k, Matrix *C_k,
                                  int *n_co, int *n_co_, int *n_co_tilde, int *n_co_tilde_tilde);

/* =========================================================================
 * L1/L5/L6: Gramians
 * ========================================================================= */

/**
 * @brief Compute controllability Gramian Wc = integral_0^inf e^{At} B B^T e^{A^T t} dt.
 * Satisfies Lyapunov equation: A Wc + Wc A^T + B B^T = 0.
 */
int control_controllability_gramian(const Matrix *A, const Matrix *B, Matrix *Wc);

/**
 * @brief Compute observability Gramian Wo from A^T Wo + Wo A + C^T C = 0.
 */
int control_observability_gramian(const Matrix *A, const Matrix *C, Matrix *Wo);

/**
 * @brief Hankel singular values: sqrt(eigenvalues(Wc * Wo)).
 * Used for balanced realization and model reduction (L8).
 */
int control_hankel_singular_values(const Matrix *Wc, const Matrix *Wo, Vector *hsv);

/* =========================================================================
 * L5/L6: Lyapunov Equations
 * ========================================================================= */

/**
 * @brief Solve continuous-time Lyapunov equation: A X + X A^T + Q = 0.
 *
 * Algorithm (L5): Bartels-Stewart method using Schur decomposition.
 * Complexity: O(n^3).
 *
 * Theorem (L4): If A is Hurwitz (all eigenvalues have negative real parts),
 * the Lyapunov equation has a unique symmetric solution X.
 * If Q > 0, then X > 0 iff A is Hurwitz.
 *
 * @param A  nxn system matrix (must be Hurwitz for unique solution)
 * @param Q  nxn symmetric RHS matrix
 * @param X  Output nxn symmetric solution matrix
 * @return   0 on success, -1 on failure
 */
int lyapunov_solve_continuous(const Matrix *A, const Matrix *Q, Matrix *X);

/**
 * @brief Solve discrete-time Lyapunov equation: A X A^T - X + Q = 0.
 *
 * Theorem (L4): If A is Schur stable (all eigenvalues inside unit circle),
 * there exists a unique symmetric solution X. If Q > 0, X > 0 iff A is Schur.
 */
int lyapunov_solve_discrete(const Matrix *A, const Matrix *Q, Matrix *X);

/* =========================================================================
 * L5/L6: Algebraic Riccati Equations
 * ========================================================================= */

/**
 * @brief Solve continuous-time algebraic Riccati equation (CARE):
 * A^T P + P A - P B R^{-1} B^T P + Q = 0.
 *
 * Algorithm (L5): Hamiltonian matrix / Schur vector method (Laub, 1979).
 * Complexity: O(n^3).
 *
 * Theorem (L4): If (A,B) is stabilizable and (A,Q^{1/2}) is detectable,
 * there exists a unique stabilizing solution P >= 0.
 *
 * (L6: Core of LQR optimal control design)
 */
int riccati_solve_care(const Matrix *A, const Matrix *B,
                        const Matrix *Q, const Matrix *R, Matrix *P);

/**
 * @brief Solve discrete-time algebraic Riccati equation (DARE):
 * A^T P A - P - A^T P B (R + B^T P B)^{-1} B^T P A + Q = 0.
 */
int riccati_solve_dare(const Matrix *A, const Matrix *B,
                        const Matrix *Q, const Matrix *R, Matrix *P);

/* =========================================================================
 * L6: Pole Placement
 * ========================================================================= */

/**
 * @brief Pole placement via Ackermann's formula for SISO systems.
 *
 * Problem (L6): Given (A,b) controllable, find K such that
 * eigenvalues(A - b*K) = {lambda_1, ..., lambda_n}.
 *
 * Algorithm: K = e_n^T * C^{-1} * p(A)
 * where C = controllability matrix, p(s) = prod(s - lambda_i).
 *
 * @param A          nxn system matrix
 * @param b          nx1 input vector (SISO)
 * @param poles      Desired pole locations (n elements)
 * @param K          Output 1xn state feedback gain
 * @return           0 on success, -1 if uncontrollable
 */
int pole_placement_ackermann(const Matrix *A, const Matrix *b,
                              const Vector *poles, Matrix *K);

/**
 * @brief Pole placement via Bass-Gura formula.
 * Alternative to Ackermann, often better conditioned.
 */
int pole_placement_bass_gura(const Matrix *A, const Matrix *b,
                              const Vector *poles, Matrix *K);

/* =========================================================================
 * L7: State Feedback and Observer Design
 * ========================================================================= */

/**
 * @brief Compute LQR gain K = R^{-1} B^T P where P solves CARE.
 * (L7: Application - optimal state feedback control)
 */
int lqr_gain(const Matrix *A, const Matrix *B, const Matrix *Q,
             const Matrix *R, Matrix *K);

/**
 * @brief Design Luenberger observer gain L for system (A,C).
 * Uses duality: place eigenvalues of (A^T - C^T L^T) = A - L C.
 * (L7: Application - state estimation)
 */
int observer_gain_luenberger(const Matrix *A, const Matrix *C,
                              const Vector *poles, Matrix *L);

/**
 * @brief Design Kalman filter gain: L = P C^T R^{-1} where P solves CARE.
 * (L7: Application - optimal state estimation under noise)
 */
int kalman_gain(const Matrix *A, const Matrix *C, const Matrix *Q,
                const Matrix *R, Matrix *L);

/* =========================================================================
 * L7/L8: Balanced Realization and Model Reduction
 * ========================================================================= */

/**
 * @brief Compute balanced realization transformation.
 * Finds T such that controllability and observability Gramians are
 * equal and diagonal: T Wc T^T = T^{-T} Wo T^{-1} = Sigma.
 * (L8: Advanced method - balanced truncation for model reduction)
 */
int balanced_realization(const Matrix *A, const Matrix *B, const Matrix *C,
                          Matrix *T, Matrix *A_bal, Matrix *B_bal,
                          Matrix *C_bal, Vector *hsv);

/**
 * @brief Balanced truncation model reduction.
 * Reduces system order from n to r by keeping states with largest
 * Hankel singular values. (L8: Model reduction)
 */
int balanced_truncation(const Matrix *A, const Matrix *B, const Matrix *C,
                         int r, Matrix *A_r, Matrix *B_r, Matrix *C_r);

#endif /* MINI_LINALG_CONTROL_LINALG_H */