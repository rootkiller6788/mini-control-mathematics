# Knowledge Graph — mini-ode-fundamentals

## L1: Definitions (Complete)

| # | Concept | C Implementation | Lean Equivalent |
|---|---------|-----------------|-----------------|
| 1 | ODE classification (linear/nonlinear, autonomous) | `ODEClass` enum, `ode_types.h` | — |
| 2 | ODE order (first/second/higher/system) | `ODEOrder` enum, `ode_types.h` | — |
| 3 | Stability types (asymptotic/marginal/unstable) | `StabilityType` enum, `ode_types.h` | — |
| 4 | 2D equilibrium classification | `EquilibriaClassification`, `ode_linear.h` | — |
| 5 | Numerical method enumeration | `ODEMethod` enum, `ode_types.h` | — |
| 6 | Bifurcation type enumeration | `BifurcationType`, `ode_types.h` | — |
| 7 | ODE vector type | `ODEVector`, `ode_types.h` | — |
| 8 | Matrix type | `Matrix`, `ode_types.h` | — |
| 9 | IVP specification (f, t0, y0, dim) | `ODEIVP` struct, `ode_types.h` | — |
| 10 | Simulation configuration | `ODESimConfig`, `ode_types.h` | — |
| 11 | Solution trajectory | `ODESolution`, `ode_types.h` | — |
| 12 | Linear ODE system | `LinearODESystem`, `ode_types.h` | — |
| 13 | Phase portrait | `PhasePortrait`, `ode_types.h` | — |
| 14 | Vector field | `VectorField`, `ode_types.h` | — |
| 15 | Wronskian structure | `WronskianInfo`, `ode_types.h` | — |
| 16 | Canonical system parameters (oscillator, pendulum, LV, Lorenz) | `ode_types.h` | — |
| 17 | Application parameters (DC motor, RLC) | `ode_types.h` | — |
| 18 | Lyapunov function candidate | `LyapunovFunction`, `ode_types.h` | — |
| 19 | Bifurcation result | `BifurcationResult`, `ode_types.h` | — |

**Status: COMPLETE** (19 definitions, all with C struct/typedef)

## L2: Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Lipschitz continuity estimation | `ode_lipschitz_estimate`, `ode_ivp.c` |
| 2 | One-sided Lipschitz condition | `ode_one_sided_lipschitz`, `ode_ivp.c` |
| 3 | Superposition principle verification | `ode_verify_superposition`, `ode_linear.c` |
| 4 | Eigenvalue-based stability | `ode_stability_from_eigenvalues`, `ode_linear.c` |
| 5 | Dominant eigenvalue / time constant | `ode_dominant_eigenvalue`, `ode_linear.c` |
| 6 | 2D equilibrium classification | `ode_classify_equilibrium_2d`, `ode_linear.c` |
| 7 | Equilibrium finding (Newton) | `ode_find_equilibria`, `ode_nonlinear.c` |
| 8 | Nullcline computation | `ode_compute_nullclines`, `ode_nonlinear.c` |
| 9 | Phase portrait generation | `ode_generate_phase_portrait`, `ode_nonlinear.c` |
| 10 | Limit cycle detection (Poincaré map) | `ode_detect_limit_cycle`, `ode_nonlinear.c` |

**Status: COMPLETE**

## L3: Mathematical Structures (Complete)

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Vector math operations (norm, dot, axpy) | `ode_ivp.c`, `ode_linear.c` (static helpers) |
| 2 | Matrix operations (mult, add, scale, transpose) | `ode_linear.c` |
| 3 | Matrix determinant (2×2, 3×3) | `ode_linear.c` |
| 4 | Matrix trace and norm | `ode_linear.c` |
| 5 | QR decomposition (Householder) | `ode_qr_decompose`, `ode_linear.c` |
| 6 | Hessenberg reduction | `ode_linear.c` |
| 7 | QR algorithm for eigenvalues | `ode_eigen_decompose`, `ode_linear.c` |
| 8 | State transition matrix Φ(t,t₀) | `ode_state_transition`, `ode_linear.c` |
| 9 | Gaussian elimination with partial pivoting | `mat_solve_matrix`, `ode_linear.c` |
| 10 | Thomas algorithm (tridiagonal) | `ode_bvp_fd_linear`, `ode_numerical.c` |

**Status: COMPLETE**

## L4: Fundamental Laws (Complete)

| # | Theorem / Law | C Implementation |
|---|---------------|-----------------|
| 1 | Picard-Lindelöf theorem (iterative construction) | `ode_picard_iteration`, `ode_picard_solve`, `ode_ivp.c` |
| 2 | Gronwall's inequality (solution bounds) | `ode_gronwall_bound`, `ode_ivp.c` |
| 3 | Continuous dependence on initial conditions | `ode_continuous_dependence_bound`, `ode_ivp.c` |
| 4 | Matrix exponential (scaling-and-squaring Padé) | `ode_matrix_exp`, `ode_linear.c` |
| 5 | Liouville's formula: W(t) = W(0)·exp(tr(A)·t) | `ode_liouville_formula`, `ode_linear.c` |
| 6 | Wronskian computation | `ode_wronskian_compute`, `ode_linear.c` |
| 7 | Floquet theory (periodic systems) | `ode_floquet_analysis`, `ode_linear.c` |
| 8 | Sturm separation theorem | `ode_sturm_separation`, `ode_linear.c` |
| 9 | Lyapunov's direct method validation | `ode_lyapunov_validate`, `ode_stability.c` |
| 10 | Lyapunov equation (Bartels-Stewart) | `ode_lyapunov_equation`, `ode_stability.c` |
| 11 | LaSalle's invariance principle | `ode_lasalle_check`, `ode_stability.c` |
| 12 | Chetaev instability theorem | `ode_chetaev_check`, `ode_stability.c` |
| 13 | Poincaré index theorem | `ode_poincare_index`, `ode_nonlinear.c` |

**Status: COMPLETE** (13 theorems implemented)

## L5: Computational Methods (Complete)

| # | Method | Implementation |
|---|--------|---------------|
| 1 | Forward Euler (O(h)) | `ode_euler_forward_step`, `ode_numerical.c` |
| 2 | Backward Euler (A-stable, Newton) | `ode_euler_backward_step`, `ode_numerical.c` |
| 3 | Heun's method (RK2, O(h²)) | `ode_heun_step`, `ode_numerical.c` |
| 4 | Classical RK4 (O(h⁴)) | `ode_rk4_step`, `ode_numerical.c` |
| 5 | Dormand-Prince RK5(4) (adaptive) | `ode_rk45_step`, `ode_numerical.c` |
| 6 | Adams-Bashforth 2-step (O(h²)) | `ode_ab2_step`, `ode_numerical.c` |
| 7 | Adams-Bashforth 3-step (O(h³)) | `ode_ab3_step`, `ode_numerical.c` |
| 8 | Adams-Bashforth 4-step (O(h⁴)) | `ode_ab4_step`, `ode_numerical.c` |
| 9 | Adams-Moulton 2-step (trapezoidal) | `ode_am2_step`, `ode_numerical.c` |
| 10 | Adams-Moulton 3-step | `ode_am3_step`, `ode_numerical.c` |
| 11 | BDF-2 (A-stable for stiff) | `ode_bdf2_step`, `ode_numerical.c` |
| 12 | BDF-3 | `ode_bdf3_step`, `ode_numerical.c` |
| 13 | Fixed-step integration driver | `ode_integrate_fixed_step`, `ode_numerical.c` |
| 14 | Adaptive integration driver | `ode_integrate_adaptive`, `ode_numerical.c` |
| 15 | Linear shooting BVP | `ode_bvp_shooting_linear`, `ode_numerical.c` |
| 16 | Nonlinear shooting BVP (Newton) | `ode_bvp_shooting_nonlinear`, `ode_numerical.c` |
| 17 | Linear finite difference BVP | `ode_bvp_fd_linear`, `ode_numerical.c` |
| 18 | Richardson extrapolation error estimate | `ode_richardson_error`, `ode_numerical.c` |

**Status: COMPLETE** (18 methods implemented)

## L6: Canonical Systems (Complete)

| # | System | Implementation |
|---|--------|---------------|
| 1 | Harmonic oscillator (damped/forced) | `ode_rhs_harmonic_oscillator`, `ode_nonlinear.c` |
| 2 | Nonlinear pendulum | `ode_rhs_pendulum`, `ode_nonlinear.c` |
| 3 | Van der Pol oscillator | `ode_rhs_vanderpol`, `ode_nonlinear.c` |
| 4 | Lotka-Volterra (predator-prey) | `ode_rhs_lotka_volterra`, `ode_nonlinear.c` |
| 5 | Lorenz system (chaos) | `ode_rhs_lorenz`, `ode_nonlinear.c` |

**Status: COMPLETE** (5 canonical systems, each with examples in examples/)

## L7: Applications (Complete)

| # | Application | Implementation |
|---|-------------|---------------|
| 1 | DC motor step response / transfer function | `ode_dc_motor_simulate`, `ode_dc_motor_transfer_function`, `ode_applications.c` |
| 2 | RLC circuit transient analysis | `ode_rlc_simulate`, `ode_rlc_damping_analysis`, `ode_applications.c` |
| 3 | Spring-mass-damper design | `ode_spring_mass_simulate`, `ode_oscillator_frequency_response`, `ode_applications.c` |
| 4 | Population dynamics (Lotka-Volterra) | `ode_predator_prey_simulate`, `ode_applications.c` |
| 5 | Pendulum period analysis | `ode_pendulum_simulate`, `ode_pendulum_period`, `ode_applications.c` |

**Status: COMPLETE** (5 applications with real parameters)

## L8: Advanced Topics (Partial)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Lyapunov exponent (maximal) | `ode_lyapunov_exponent`, `ode_nonlinear.c` |
| 2 | Lyapunov spectrum (QR method) | `ode_lyapunov_spectrum`, `ode_nonlinear.c` |
| 3 | 1D bifurcation detection | `ode_detect_bifurcation_1d`, `ode_nonlinear.c` |
| 4 | Hopf bifurcation detection (2D) | `ode_detect_hopf_2d`, `ode_nonlinear.c` |
| 5 | Region of attraction estimation | `ode_estimate_roa`, `ode_stability.c` |
| 6 | Input-to-state stability (ISS) | `ode_check_iss`, `ode_stability.c` |

**Status: PARTIAL** (6/6 implemented, some with simplified numerics)

## L9: Research Frontiers (Partial — documented only)

| # | Topic | Status |
|---|-------|--------|
| 1 | Delay differential equations | Documented in course-tree.md |
| 2 | Stochastic differential equations | Documented, no implementation |
| 3 | Fractional-order ODEs | Documented, no implementation |
| 4 | Neural ODEs | Documented, no implementation |

**Status: PARTIAL** (documented, not implemented)
