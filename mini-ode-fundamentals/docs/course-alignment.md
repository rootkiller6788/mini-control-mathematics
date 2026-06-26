# Course Alignment — mini-ode-fundamentals

## Nine-School Curriculum Mapping

### MIT
- **6.241 Dynamic Systems & Control**: State-space IVP formulation → `ODEIVP`, `ode_ivp_create`
- **6.241**: Matrix exponential, state transition matrix → `ode_matrix_exp`, `ode_state_transition`
- **6.241**: Eigenvalue-based stability → `ode_stability_from_eigenvalues`
- **6.243 Nonlinear Systems**: Lyapunov functions, LaSalle → `ode_lyapunov_validate`, `ode_lasalle_check`

### Stanford
- **ENGR 207B Linear Dynamical Systems**: Modal decomposition → `ode_eigen_decompose`
- **ENGR 209A Nonlinear Dynamics**: Phase portraits, bifurcations → `ode_generate_phase_portrait`, `ode_detect_hopf_2d`
- **ME 346A**: Heat equation ODE reduction → `ode_bvp_fd_linear`

### Berkeley
- **ME 132 Dynamic Systems**: Canonical systems (oscillator, pendulum) → `ode_rhs_harmonic_oscillator`, `ode_rhs_pendulum`
- **ME 232 Advanced Control**: Lyapunov equation → `ode_lyapunov_equation`
- **ME 234 Nonlinear**: Limit cycles → `ode_detect_limit_cycle`

### Caltech
- **CDS 101/110**: Numerical integration fundamentals → `ode_euler_forward_step`, `ode_rk4_step`
- **CDS 212 Robust**: ISS → `ode_check_iss`
- **CDS 240 Optimal**: Shooting methods → `ode_bvp_shooting_linear`

### ETH
- **151-0591 Control I**: RLC circuit, DC motor → `ode_rhs_rlc_circuit`, `ode_rhs_dc_motor`
- **151-0555 Linear Systems**: Fundamental matrix, Floquet → `ode_floquet_analysis`
- **151-0563 Robust**: Absolute stability → `ode_absolute_stability_region`

### Cambridge
- **3F2 Systems & Control**: ODE classification, Lipschitz → `ode_lipschitz_estimate`
- **4F2 Robust**: Lyapunov stability → `ode_check_lyapunov_stability`
- **4F3 Nonlinear**: Bifurcation, chaos → `ode_rhs_lorenz`, `ode_lyapunov_exponent`

### Georgia Tech
- **ECE 6550 Nonlinear**: Van der Pol, pendulum → `ode_rhs_vanderpol`, `ode_rhs_pendulum`
- **AE 6530 Optimal**: BVP shooting → `ode_bvp_shooting_nonlinear`
- **ME 6401 Linear**: Modal analysis → `ode_dominant_eigenvalue`

### Purdue
- **ECE 602 Lumped Systems**: RLC circuit analysis → `ode_rlc_damping_analysis`
- **ME 575 Industrial**: DC motor model → `ode_dc_motor_transfer_function`
- **ME 675 Multivariable**: Matrix exponential → `ode_matrix_exp`

### 清华 (Tsinghua)
- **自动控制原理**: Transfer function, step response → `ode_dc_motor_transfer_function`
- **现代控制理论**: State space, Lyapunov equation → `ode_lyapunov_equation`
- **非线性控制**: Phase plane analysis → `ode_classify_equilibrium_2d`

## Coverage Summary

| School | Courses Mapped | Coverage |
|--------|---------------|----------|
| MIT | 6.241, 6.243 | ✓ Full |
| Stanford | ENGR 207B, 209A | ✓ Full |
| Berkeley | ME 132, 232, 234 | ✓ Full |
| Caltech | CDS 101, 212, 240 | ✓ Full |
| ETH | 151-0591, 0555, 0563 | ✓ Full |
| Cambridge | 3F2, 4F2, 4F3 | ✓ Full |
| Georgia Tech | ECE 6550, AE 6530, ME 6401 | ✓ Full |
| Purdue | ECE 602, ME 575, 675 | ✓ Full |
| 清华 | 自动控制, 现代控制, 非线性 | ✓ Full |
