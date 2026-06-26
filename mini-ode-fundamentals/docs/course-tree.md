# Course Tree (Prerequisite Dependencies) — mini-ode-fundamentals

## Dependency Graph

```
Calculus (limits, derivatives, integrals)
  │
  ├─► L1: ODE Definitions
  │     └─► L2: Core Concepts (Lipschitz, superposition, equilibria)
  │           ├─► L3: Mathematical Structures (matrices, QR, eigenvalues)
  │           │     └─► L4: Fundamental Laws (matrix exp, Lyapunov eq, Floquet)
  │           │           ├─► L5: Computational Methods (RK, AB, BDF, BVP)
  │           │           │     ├─► L6: Canonical Systems (oscillator, Lorenz, LV)
  │           │           │     │     └─► L7: Applications (DC motor, RLC, population)
  │           │           │     │           └─► L8: Advanced (Lyapunov exp, bifurcation)
  │           │           │     │                 └─► L9: Frontiers (DDE, SDE, neural ODE)
  │           │           │     └─► L7: Stiffness detection, method selection
  │           │           └─► L8: ROA estimation, ISS
  │           └─► L8: Poincaré index, limit cycles
  └─► Linear Algebra (eigenvalues, matrix decompositions)
        └─► L3, L4, L5 (all numerical linear algebra)
```

## Prerequisites for This Module

1. **Single-variable calculus**: Limits, derivatives, integrals, Taylor series
2. **Linear algebra**: Matrix multiplication, determinants, eigenvalues/eigenvectors, Gaussian elimination
3. **Basic programming**: C99, dynamic memory allocation, function pointers

## What This Module Enables (Downstream)

1. **mini-complex-analysis**: Eigenvalues for stability of complex-valued systems
2. **mini-laplace-z-transform**: ODE-to-transfer-function conversion
3. **mini-linear-algebra-control**: Matrix exponential, Lyapunov equation (shared infrastructure)
4. **mini-state-space-theory**: State transition matrix, controllability/observability
5. **mini-nonlinear-control** (future): Lyapunov functions, bifurcation analysis
6. **mini-optimal-control** (future): Shooting methods for two-point BVPs

## L9 Research Frontiers Detail

### Delay Differential Equations (DDEs)
- y'(t) = f(t, y(t), y(t-τ))
- Require history functions on [t₀-τ, t₀]
- Method of steps, continuous Runge-Kutta extensions
- Reference: Bellen & Zennaro "Numerical Methods for Delay Differential Equations" (2003)

### Stochastic Differential Equations (SDEs)
- dX_t = a(X_t)dt + b(X_t)dW_t
- Euler-Maruyama: X_{n+1} = X_n + a(X_n)·Δt + b(X_n)·√Δt·ξ_n
- Milstein method for better strong convergence
- Reference: Kloeden & Platen "Numerical Solution of SDEs" (1992)

### Fractional-Order ODEs
- D^α y(t) = f(t, y), α ∈ (0,2)
- Grünwald-Letnikov discretization
- Memory effects in viscoelasticity and control

### Neural ODEs
- dz/dt = f_θ(z(t), t) where f_θ is a neural network
- Adjoint sensitivity method for backpropagation
- Reference: Chen et al. "Neural Ordinary Differential Equations" (NeurIPS 2018)
