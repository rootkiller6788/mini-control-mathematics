# Course Dependency Tree

## Prerequisites
- Calculus (multivariate, Taylor series)
- Linear Algebra (vector spaces, eigenvalues)
- Differential Equations (ODE theory)
- Probability Theory (stochastic processes)
- Control Theory (state-space fundamentals)

## Internal Dependencies
numerical_core -> linear_systems -> eigen_methods -> control_solvers
numerical_core -> integration, optimization, interpolation, nonlinear_methods

## Upward Dependencies
- mini-model-predictive-control (QP, ODE integration)
- mini-robust-control (H-inf, LMI)
- mini-system-identification (least squares, SVD)
- mini-nonlinear-control (Newton-Kantorovich, shooting)
- mini-stochastic-control (Kalman, Riccati)

## Key Theorems
- Kantorovich theorem (Newton convergence)
- Eckart-Young theorem (optimal low-rank)
- Francis QR algorithm (eigenvalue computation)
- Bartels-Stewart algorithm (Sylvester equation)
- Laub CARE method (Riccati solution)
- Kalman controllability criterion
