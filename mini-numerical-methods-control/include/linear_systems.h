/**
 * linear_systems.h — Numerical Linear Algebra for Control
 *
 * Covers: L3 Math Structures, L4 Conservation/Linear Constraints, L5 Methods
 *
 * Solving linear systems Ax = b is fundamental to control: state estimation,
 * MPC optimization, and linearized dynamics. Matrix decompositions are the
 * computational backbone of modern control algorithms.
 *
 * Reference: Golub & Van Loan, Matrix Computations, 4th ed. (2013)
 *            Demmel, Applied Numerical Linear Algebra (1997)
 *            Stewart, Matrix Algorithms (1998)
 */

#ifndef LINEAR_SYSTEMS_H
#define LINEAR_SYSTEMS_H

#include "numerical_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L3: Decomposition Result Types
 * ========================================================================== */

/**
 * LUResult — Stores LU decomposition A = P·L·U
 *
 * L3 Structure: For a square matrix A, Gaussian elimination with partial
 * pivoting produces A = P·L·U where P is permutation, L is unit lower
 * triangular, U is upper triangular.
 */
typedef struct {
    Matrix *L;          /**< Unit lower triangular (n×n) */
    Matrix *U;          /**< Upper triangular (n×n) */
    size_t *pivots;     /**< Row permutation indices (0-based), length n */
    size_t  n;          /**< Matrix dimension */
    int     singular;   /**< 1 if matrix is numerically singular */
} LUResult;

/**
 * QRResult — Stores QR decomposition A = Q·R
 *
 * L3 Structure: Q is orthogonal (QᵀQ = I), R is upper triangular.
 * Computed via Householder reflections.
 */
typedef struct {
    Matrix *Q;          /**< Orthogonal matrix (m×m) */
    Matrix *R;          /**< Upper triangular (m×n) */
    size_t  m, n;       /**< Dimensions */
} QRResult;

/**
 * CholeskyResult — Stores Cholesky decomposition A = L·Lᵀ (A ≻ 0)
 *
 * L3 Structure: For symmetric positive definite A, the Cholesky decomposition
 * is A = L·Lᵀ where L is lower triangular. Key to efficient solution of
 * Lyapunov equations and Kalman filtering.
 */
typedef struct {
    Matrix *L;          /**< Lower triangular (n×n) */
    size_t  n;          /**< Dimension */
} CholeskyResult;

/* ==========================================================================
 * L5: Direct Linear System Solvers
 * ========================================================================== */

/**
 * gaussian_elimination_solve — solve Ax = b by Gaussian elimination with
 * partial pivoting.
 *
 * L5 Method: O(n³). Augments [A|b], forward eliminates to upper triangular
 * form, back-substitutes. The foundational direct method.
 *
 * Returns: 0 on success, -1 if singular.
 */
int gaussian_elimination_solve(const Matrix *A, const Vector *b, Vector *x);

/**
 * lu_decompose — compute LU decomposition A = P·L·U
 *
 * L5 Method: O(2n³/3) for n×n matrix. Partial pivoting for numerical stability.
 * L and U stored compactly in the output matrices.
 */
LUResult* lu_decompose(const Matrix *A);

/**
 * lu_solve — solve Ax = b using precomputed LU decomposition
 *
 * L5 Method: O(n²). Forward substitution (Ly = Pᵀb) + backward substitution (Ux = y).
 * Much faster than re-decomposing for multiple right-hand sides.
 */
Vector* lu_solve(const LUResult *lu, const Vector *b);

/**
 * lu_determinant — compute det(A) from LU decomposition
 *
 * L5 Method: det(A) = det(P)·det(L)·det(U) = (-1)^(#swaps) · Π u_{ii}
 */
double lu_determinant(const LUResult *lu);

/**
 * lu_free — free LU decomposition
 */
void lu_free(LUResult *lu);

/**
 * cholesky_decompose — Cholesky decomposition A = L·Lᵀ for A ≻ 0
 *
 * L5 Method: O(n³/3). Twice as fast as LU for SPD matrices.
 * Essential for solving Lyapunov equations and Riccati equations in control.
 *
 * Returns NULL if A is not positive definite.
 */
CholeskyResult* cholesky_decompose(const Matrix *A);

/**
 * cholesky_solve — solve Ax = b using Cholesky decomposition
 *
 * L5 Method: O(n²). Forward (Ly = b) + backward (Lᵀx = y) substitution.
 */
Vector* cholesky_solve(const CholeskyResult *chol, const Vector *b);

/**
 * cholesky_logdet — log|det(A)| = 2 Σ log(L_{ii}) for numerical stability
 *
 * L3 Engineering Quantity: Used in maximum likelihood estimation for
 * Kalman filters and system identification.
 */
double cholesky_logdet(const CholeskyResult *chol);

/**
 * cholesky_free — free Cholesky decomposition
 */
void cholesky_free(CholeskyResult *chol);

/* ==========================================================================
 * L5: Orthogonal Decompositions
 * ========================================================================== */

/**
 * qr_decompose_householder — QR decomposition via Householder reflections
 *
 * L5 Method: O(2mn²) for m×n matrix (m ≥ n). Numerically stable.
 * Householder reflection: H = I - 2vvᵀ/(vᵀv) maps a vector to a multiple of e₁.
 */
QRResult* qr_decompose_householder(const Matrix *A);

/**
 * qr_solve — solve least squares min ||Ax - b||₂ via QR
 *
 * L5 Method: O(mn²). For overdetermined systems (m > n), solves
 * R₁x = Qᵀb where R₁ is the upper n×n block of R.
 * This is the standard method for least squares problems in control.
 */
Vector* qr_solve(const QRResult *qr, const Vector *b);

/**
 * qr_free — free QR decomposition
 */
void qr_free(QRResult *qr);

/* ==========================================================================
 * L5: Special Linear System Solvers for Control
 * ========================================================================== */

/**
 * solve_triangular_upper — solve Ux = b where U is upper triangular
 *
 * L5 Method: Back substitution, O(n²). Used in LQR Riccati solvers.
 */
void solve_triangular_upper(const Matrix *U, const Vector *b, Vector *x);

/**
 * solve_triangular_lower — solve Lx = b where L is lower triangular
 *
 * L5 Method: Forward substitution, O(n²).
 */
void solve_triangular_lower(const Matrix *L, const Vector *b, Vector *x);

/**
 * solve_sylvester_equation — solve AX + XB = C (Sylvester equation)
 *
 * L5 Method: Bartels-Stewart algorithm using real Schur decomposition.
 * O(n³). Key to solving Lyapunov equations (B = Aᵀ for continuous case).
 *
 * L6 Problem: Stability of linear systems ẋ = Ax depends on solving
 * AᵀP + PA = -Q (Lyapunov equation) — a special Sylvester equation.
 */
Matrix* solve_sylvester_equation(const Matrix *A, const Matrix *B, const Matrix *C);

/**
 * solve_lyapunov_continuous — solve AᵀP + PA + Q = 0 (continuous Lyapunov)
 *
 * L6 Problem: The matrix P ≻ 0 proves A is Hurwitz (stable). Used in
 * controllability/observability Gramian computation and LQR design.
 */
Matrix* solve_lyapunov_continuous(const Matrix *A, const Matrix *Q);

/**
 * solve_lyapunov_discrete — solve AᵀPA - P + Q = 0 (discrete Lyapunov)
 *
 * L6 Problem: The matrix P ≻ 0 proves A is Schur (discrete-time stable).
 * Used in discrete Kalman filter and DLQR design.
 */
Matrix* solve_lyapunov_discrete(const Matrix *A, const Matrix *Q);

/**
 * solve_riccati_care — solve AᵀP + PA - PBR⁻¹BᵀP + Q = 0 (CARE)
 *
 * L6 Problem: Continuous Algebraic Riccati Equation. The solution P
 * gives the optimal LQR gain K = R⁻¹BᵀP. Solved via Hamiltonian/Schur method.
 *
 * Complexity: O(n³) using real Schur decomposition.
 */
Matrix* solve_riccati_care(const Matrix *A, const Matrix *B,
                           const Matrix *Q, const Matrix *R);

/**
 * solve_riccati_dare — solve AᵀPA - P - AᵀPB(R + BᵀPB)⁻¹BᵀPA + Q = 0 (DARE)
 *
 * L6 Problem: Discrete Algebraic Riccati Equation for discrete-time LQR.
 */
Matrix* solve_riccati_dare(const Matrix *A, const Matrix *B,
                           const Matrix *Q, const Matrix *R);

/* ==========================================================================
 * L8: Advanced Factorizations
 * ========================================================================== */

/**
 * svd_decompose — compute thin SVD A = U Σ Vᵀ
 *
 * L8 Method: Golub-Kahan-Reinsch algorithm (bidiagonalization + QR iteration).
 * O(mn²) for thin SVD. The most informative matrix decomposition — reveals
 * rank, range, nullspace, and optimal low-rank approximation.
 *
 * L8 Application: Model reduction via balanced truncation, PCA for
 * process monitoring, and total least squares for EIV identification.
 */
SVDResult* svd_decompose(const Matrix *A);

/**
 * svd_pseudo_inverse — compute Moore-Penrose pseudoinverse A⁺ = V Σ⁺ Uᵀ
 *
 * L8 Method: For rank-deficient or rectangular systems, A⁺ gives the
 * minimum-norm least-squares solution. Used in control allocation.
 */
Matrix* svd_pseudo_inverse(const SVDResult *svd);

/**
 * svd_rank — compute effective numerical rank of A
 *
 * L3 Engineering Quantity: rank(A) = number of singular values > ε·σ₁.
 * Critical for controllability/observability analysis.
 */
size_t svd_rank(const SVDResult *svd, double tol);

/**
 * svd_low_rank_approx — best rank-k approximation in Frobenius norm
 *
 * L8 Method: Eckart-Young theorem — Σ_{i=1}^{k} σᵢ uᵢ vᵢᵀ is optimal.
 * Foundation of POD (Proper Orthogonal Decomposition) for model reduction.
 */
Matrix* svd_low_rank_approx(const SVDResult *svd, size_t k);

/**
 * svd_free — free SVD result
 */
void svd_free(SVDResult *svd);

#ifdef __cplusplus
}
#endif

#endif /* LINEAR_SYSTEMS_H */
