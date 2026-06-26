/**
 * control_solvers.h — Control-Specific Numerical Solvers
 *
 * Covers: L6 Engineering Problems, L7 Applications
 *
 * Implements numerical algorithms for canonical control problems:
 * LQR design, Kalman filtering, pole placement, controllability/observability
 * analysis, and H2/H∞ norm computation.
 *
 * Reference: Laub, "A Schur method for solving algebraic Riccati equations" (1979)
 *            Van Dooren, "A generalized eigenvalue approach..." (1981)
 *            Zhou, Doyle, Glover, Robust and Optimal Control (1996)
 */

#ifndef CONTROL_SOLVERS_H
#define CONTROL_SOLVERS_H

#include "numerical_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Control System Representations
 * ========================================================================== */

/**
 * LTISystem — continuous-time linear time-invariant (LTI) system
 *
 * L1 Definition: ẋ = Ax + Bu, y = Cx + Du
 * A ∈ ℝ^{n×n}: state matrix
 * B ∈ ℝ^{n×m}: input matrix
 * C ∈ ℝ^{p×n}: output matrix
 * D ∈ ℝ^{p×m}: feedthrough matrix
 */
typedef struct {
    Matrix *A;
    Matrix *B;
    Matrix *C;
    Matrix *D;
    size_t  n;   /**< Number of states */
    size_t  m;   /**< Number of inputs */
    size_t  p;   /**< Number of outputs */
} LTISystem;

/**
 * DiscreteLTISystem — discrete-time LTI system
 *
 * L1 Definition: x_{k+1} = A_d x_k + B_d u_k, y_k = C_d x_k + D_d u_k
 */
typedef struct {
    Matrix *A;
    Matrix *B;
    Matrix *C;
    Matrix *D;
    size_t  n, m, p;
} DiscreteLTISystem;

/**
 * LQRResult — Linear Quadratic Regulator solution
 *
 * L6 Problem: min ∫₀^∞ (xᵀQx + uᵀRu) dt  subject to ẋ = Ax + Bu
 * Solution: u = -Kx with K = R⁻¹BᵀP where P solves CARE.
 */
typedef struct {
    Matrix *K;       /**< Optimal state feedback gain (m×n) */
    Matrix *P;       /**< CARE solution (n×n) */
    double  J_opt;   /**< Optimal cost J* = x₀ᵀPx₀ */
    int     solved;   /**< 1 if CARE solved successfully */
} LQRResult;

/**
 * KalmanFilter — Kalman filter data structure
 *
 * L6 Problem: Estimate state x̂ of ẋ = Ax + Bu + w, y = Cx + v
 * where w ~ N(0, Q) is process noise, v ~ N(0, R) is measurement noise.
 */
typedef struct {
    Matrix *x_hat;     /**< Current state estimate (n×1 as Matrix) */
    Matrix *P;         /**< Error covariance (n×n) */
    Matrix *A;         /**< System dynamics */
    Matrix *B;         /**< Input matrix */
    Matrix *C;         /**< Output matrix */
    Matrix *Q;         /**< Process noise covariance */
    Matrix *R;         /**< Measurement noise covariance */
    size_t  n, m, p;
} KalmanFilter;

/* ==========================================================================
 * L6: LQR Design
 * ========================================================================== */

/**
 * lqr_design — compute LQR gain K for continuous-time system
 *
 * L6 Problem: Solve CARE: AᵀP + PA - PBR⁻¹BᵀP + Q = 0
 * Then K = R⁻¹BᵀP.
 *
 * L7 Application: Satellite attitude control, quadrotor stabilization,
 * Apollo Lunar Module digital autopilot.
 */
LQRResult* lqr_design(const Matrix *A, const Matrix *B,
                      const Matrix *Q, const Matrix *R);

/**
 * lqr_design_discrete — compute DLQR gain for discrete-time system
 *
 * L6 Problem: Solve DARE for discrete LQR.
 */
LQRResult* lqr_design_discrete(const Matrix *A, const Matrix *B,
                               const Matrix *Q, const Matrix *R);

/**
 * lqr_result_free — deallocate LQR result
 */
void lqr_result_free(LQRResult *r);

/* ==========================================================================
 * L6: Kalman Filtering
 * ========================================================================== */

/**
 * kalman_filter_create — initialize Kalman filter
 */
KalmanFilter* kalman_filter_create(const Matrix *A, const Matrix *B,
                                   const Matrix *C, const Matrix *Q,
                                   const Matrix *R, const Matrix *x0,
                                   const Matrix *P0);

/**
 * kalman_predict — prediction step: x̂⁻ = A x̂ + B u, P⁻ = A P Aᵀ + Q
 *
 * L6 Problem: Propagate state estimate and covariance through dynamics.
 * Called at each time step before measurement arrives.
 */
void kalman_predict(KalmanFilter *kf, const Matrix *u, double dt);

/**
 * kalman_update — update step: K = P⁻ Cᵀ (C P⁻ Cᵀ + R)⁻¹
 *
 * L6 Problem: Incorporate measurement y to refine estimate.
 * x̂ = x̂⁻ + K(y - C x̂⁻), P = (I - KC)P⁻
 */
void kalman_update(KalmanFilter *kf, const Matrix *y);

/**
 * kalman_filter_free — deallocate Kalman filter
 */
void kalman_filter_free(KalmanFilter *kf);

/* ==========================================================================
 * L6: Pole Placement
 * ========================================================================== */

/**
 * pole_placement_acker — Ackermann's formula for SISO state feedback
 *
 * L6 Method: K = eₙᵀ · C⁻¹ · Δ(A) where Δ(A) = p(A) using Cayley-Hamilton.
 * Places closed-loop eigenvalues at desired locations for single-input systems.
 * Complexity: O(n³).
 *
 * L7 Application: Automotive engine control, aircraft flight control
 * (Boeing 747 control law design uses eigenvalue placement).
 */
Matrix* pole_placement_acker(const Matrix *A, const Matrix *B,
                              const double *desired_poles, size_t n_poles);

/**
 * pole_placement_bass_gura — Bass-Gura formula for SISO pole placement
 *
 * L6 Method: Transform to controllable canonical form, assign poles,
 * transform back. Alternative to Ackermann with different numerical properties.
 */
Matrix* pole_placement_bass_gura(const Matrix *A, const Matrix *B,
                                  const double *desired_poles, size_t n_poles);

/* ==========================================================================
 * L6: Controllability and Observability
 * ========================================================================== */

/**
 * controllability_matrix — compute C = [B AB A²B ... A^{n-1}B]
 *
 * L4 Conservation/L6 Problem: System is controllable iff rank(C) = n.
 * Used to check if all states can be driven to arbitrary values.
 */
Matrix* controllability_matrix(const Matrix *A, const Matrix *B);

/**
 * observability_matrix — compute O = [C; CA; CA²; ...; CA^{n-1}]
 *
 * L4 Conservation/L6 Problem: System is observable iff rank(O) = n.
 * Required for state estimator design.
 */
Matrix* observability_matrix(const Matrix *A, const Matrix *C);

/**
 * controllability_gramian — compute W_c = ∫₀^∞ e^{Aτ} B Bᵀ e^{Aᵀτ} dτ
 *
 * L6 Problem: Controllability Gramian. W_c ≻ 0 ⇔ (A, B) controllable.
 * Also bounds the minimum control energy: E_min = x₀ᵀ W_c⁻¹ x₀.
 */
Matrix* controllability_gramian(const Matrix *A, const Matrix *B);

/**
 * observability_gramian — compute W_o = ∫₀^∞ e^{Aᵀτ} Cᵀ C e^{Aτ} dτ
 *
 * L6 Problem: Observability Gramian. W_o ≻ 0 ⇔ (A, C) observable.
 * Used in balanced truncation model reduction.
 */
Matrix* observability_gramian(const Matrix *A, const Matrix *C);

/**
 * is_controllable — check controllability via rank test
 *
 * L6 Problem: rank(ctrb(A,B)) == n
 */
int is_controllable(const Matrix *A, const Matrix *B, double tol);

/**
 * is_observable — check observability via rank test
 */
int is_observable(const Matrix *A, const Matrix *C, double tol);

/* ==========================================================================
 * L7/L8: System Norms and Robust Control Metrics
 * ========================================================================== */

/**
 * h2_norm — compute H₂ norm of LTI system
 *
 * L7 Application: ||H||₂² = tr(BᵀW_o B) = tr(C W_c Cᵀ).
 * Measures RMS response to white noise. Used in LQG controller evaluation.
 */
double h2_norm(const LTISystem *sys);

/**
 * hinf_norm_bisection — compute H∞ norm via bisection on Hamiltonian eigenvalues
 *
 * L8 Method: ||H||_∞ = sup_ω σ_max(H(jω)). Computed via bisection:
 * ||H||_∞ < γ ⇔ Hamiltonian H_γ has no eigenvalues on imaginary axis.
 *
 * L7 Application: Robust stability margin, loop shaping, H∞ controller synthesis.
 * Examples: F-35 flight control, hard disk drive servo control.
 */
double hinf_norm_bisection(const LTISystem *sys, double tol, size_t max_iter);

/**
 * lti_system_free — deallocate LTI system
 */
void lti_system_free(LTISystem *sys);

/**
 * discrete_lti_system_free — deallocate discrete LTI system
 */
void discrete_lti_system_free(DiscreteLTISystem *sys);

/**
 * lti_controllability_staircase — staircase form for controllability analysis
 *
 * L8 Method: Orthogonal transformation revealing controllable/uncontrollable
 * subspaces. Tᵀ A T in block triangular form. More numerically reliable
 * than rank test on controllability matrix.
 */
int lti_controllability_staircase(Matrix *A, Matrix *B, size_t *n_controllable);

#ifdef __cplusplus
}
#endif

#endif /* CONTROL_SOLVERS_H */
