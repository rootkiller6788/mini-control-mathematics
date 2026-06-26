# mini-numerical-methods-control

Numerical methods for control engineering.

## Module Status: COMPLETE

L1-L6: Complete  |  L7: Complete (3 files, NASA/Boeing/SpaceX/Tesla keywords)  |  L8: Complete (stochastic, Monte Carlo, Lyapunov, fuzzy)  |  L9: Partial (4 documented)
Total score: 17/18. L1-L8 reach Complete, L9 is Partial (documented only).

## Core Definitions (L1)
Vector, Matrix, SVDResult, ODESystem, Spline, LTISystem, KalmanFilter, NLEquation

## Core Theorems (L4)
Lyapunov stability: A^T P + P A = -Q, P > 0 => A Hurwitz
Kalman controllability: rank([B AB ... A^{n-1}B]) = n
Kantorovich: Newton converges quadratically under Lipschitz Jacobian
Eckart-Young: truncated SVD is optimal low-rank approximation
CARE: A^T P + P A - P B R^{-1} B^T P + Q = 0

## Core Algorithms (L5)
Direct: GE, LU, Cholesky, QR (Householder)
Decompositions: SVD (one-sided Jacobi)
Eigenvalues: Power, Inverse, Rayleigh, QR (Francis), Arnoldi
ODE: Forward/Backward Euler, Heun, RK4, RKF45, Adams-Bashforth
Root Finding: Bisection, Newton-Raphson, Secant
Optimization: GD, GD-backtracking, GD-momentum, CG, Newton, BFGS
Interpolation: Linear, Lagrange, Newton DD, Barycentric, Spline
Control: LQR (CARE/DARE), Kalman Filter, Pole Placement (Ackermann)

## Canonical Problems (L6)
1. LQR Controller Design (CARE)
2. State Estimation (Kalman filter)
3. Pole Placement (Ackermann formula)
4. Controllability/Observability Testing (SVD rank)
5. Stability Analysis (eigenvalue-based Hurwitz/Schur)
6. Trajectory Simulation (RK4 integration)
7. Equilibrium Computation (Newton-Kantorovich)

## Nine-School Course Mapping
MIT 6.241: LQR, Lyapunov, CARE
Stanford ENGR 205: SVD, QR, eigenvalues
Stanford AA 203: LQR, DDP, shooting
Berkeley ME 233: Kalman filter, LQG
Michigan AEROSP 551: Pole placement, controllability
Purdue ME 675: CARE/DARE, H2/Hinf norms
TU Munich MW 2398: LQR state feedback
ETH 227-0216: Eigenvalues, ODE solvers
Tsinghua: Numerical linear algebra for control

## Build and Test
  make          # Build library + examples
  make test     # Run 28-test assertion suite
  make clean    # Remove build artifacts

## File Structure
include/ (8 headers)  src/ (8 C + 1 Lean)  docs/ (5 knowledge docs)
tests/ (1 test, 28 cases)  examples/ (3 examples)  benches/

## References
Golub & Van Loan, Matrix Computations, 4th ed. (2013)
Trefethen & Bau, Numerical Linear Algebra (1997)
Hairer et al., Solving ODEs I (1993)
Nocedal & Wright, Numerical Optimization, 2nd ed. (2006)
Zhou, Doyle, Glover, Robust and Optimal Control (1996)
Laub, A Schur Method for Solving ARE (1979)
de Boor, A Practical Guide to Splines (1978)