# Gap Report — mini-ode-fundamentals

## Identified Gaps (by priority)

### Priority 1 — Low (Minor Enhancements)

1. **Implicit method robustness**: BDF and backward Euler Newton solvers use simplified damped steps rather than full Newton. Acceptable for most problems but can fail on highly nonlinear stiff systems.

2. **Vector field Jacobian**: `ode_nonlinear.c` Newton solver requires user-provided Jacobian. No finite-difference fallback for equilibrium finding.

3. **Adaptive multistep**: AB/AM methods lack adaptive step size control with error estimation.

### Priority 2 — Medium (Would Strengthen Completeness)

4. **Delay differential equations (DDEs)**: Currently only documented. A basic constant-delay solver would bridge to L9.

5. **Differential-algebraic equations (DAEs)**: Not implemented. Important for constrained mechanical systems.

6. **Boundary value problem robustness**: Shooting method uses simplified Newton; FD method is limited to linear problems.

### Priority 3 — High (Research Frontiers)

7. **Stochastic ODE integrators**: Euler-Maruyama or Milstein methods not implemented.

8. **Neural ODE integration**: Backpropagation through ODE solvers not implemented (adjoint method).

9. **Symplectic integrators**: Not implemented (important for Hamiltonian systems).

### No Critical Gaps

All required levels (L1-L6) are complete. L7 is complete with 5 real applications. L8 has 6 advanced topics. L9 has adequate documentation.
