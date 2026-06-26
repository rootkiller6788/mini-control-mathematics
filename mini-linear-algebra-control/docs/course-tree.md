# Course Tree -- mini-linear-algebra-control

## Prerequisite Dependency Graph

```
mini-linear-algebra-control
  |
  +-- Requires: mini-calculus-variations (matrix calculus)
  |   |
  |   +-- Matrix derivatives for gradient-based methods
  |
  +-- Requires: mini-general-system-theory (state-space models)
  |   |
  |   +-- State-space representation (A,B,C,D)
  |   +-- System interconnections
  |
  +-- Provides foundation for:
      |
      +-- mini-optimal-control (LQR, LQG, H-infinity)
      +-- mini-state-estimation (Kalman filter, observers)
      +-- mini-system-identification (subspace methods)
      +-- mini-robust-control (Lyapunov-based analysis)
      +-- mini-model-predictive-control (quadratic programming)
      +-- mini-nonlinear-control (linearization, Lyapunov methods)
```

## Internal Dependency Map
```
matrix.h           vector.h
  |                   |
  v                   v
solvers.h         decompositions.h
  |                   |
  v                   v
norms.h            eigen.h
  |                   |
  +---------+---------+
            |
            v
     control_linalg.h
            |
            v
      matrix_exp.h
```

## Knowledge Progression
1. Matrix/Vector basics (matrix.h, vector.h)
2. Linear system solving (solvers.h)
3. Matrix decompositions (decompositions.h)
4. Norms and conditioning (norms.h)
5. Eigenvalue computation (eigen.h)
6. Control-specific operations (control_linalg.h)
7. Matrix functions (matrix_exp.h)