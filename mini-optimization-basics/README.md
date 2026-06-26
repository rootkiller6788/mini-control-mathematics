# mini-optimization-basics

**Numerical Optimization Fundamentals** — Convex analysis, gradient methods, Newton-type methods, constrained optimization, linear/quadratic programming.

## Module Status: COMPLETE

- L1-L6: Complete
- L7: Complete (3 applications)
- L8: Partial+ (7 advanced methods)
- L9: Partial (documented)

## Nine-Layer Knowledge Coverage

| Level | Name | Status |
|-------|------|--------|
| L1 | Definitions | Complete — 8+ struct/typedef definitions |
| L2 | Core Concepts | Complete — convexity, duality, optimality |
| L3 | Engineering Quantities | Complete — convergence rates, condition numbers |
| L4 | Conservation Laws | Complete — Fermat, KKT, duality theorems |
| L5 | Engineering Methods | Complete — 20+ algorithms implemented |
| L6 | Engineering Problems | Complete — LP, QP, NLS, composite optimization |
| L7 | Applications | Complete — 3 examples (LP, Rosenbrock, Newton) |
| L8 | Advanced Methods | Partial+ — L-BFGS, trust-region, Fenchel, Moreau |
| L9 | Research Frontiers | Partial — Documented (SGD, ADMM) |

## Core Definitions

- `vector_t`, `matrix_t` — Fundamental linear algebra types
- `convex_set_t` — Convex set with membership, projection, support oracles
- `convex_function_t` — Convex function with eval, gradient, Hessian, conjugate, proximal oracles
- `linear_program_t` — Standard-form LP (c^T x s.t. Ax=b, x>=0)
- `quadratic_program_t` — Convex QP ((1/2)x^T P x + q^T x s.t. Ax=b, x>=0)
- `general_nlp_t` — General NLP with KKT fields

## Core Theorems

- **Fermat Theorem**: ∇f(x*) = 0 at local minimum (optimality_cert_t)
- **KKT Theorem**: Necessary + sufficient for convex NLP with Slater condition
- **Weak Duality**: p* >= d* always; duality gap >= 0
- **Strong Duality**: p* = d* under Slater + convexity
- **Convergence: Gradient Descent**: O(1/k) convex, O((1-μ/L)^k) strongly convex
- **Convergence: Nesterov**: O(1/k²) — optimal for first-order methods
- **Convergence: Newton**: Quadratic locally, globalized by line search
- **Separating Hyperplane Theorem**: Disjoint convex sets separable
- **Moreau Decomposition**: v = prox_f(v) + prox_{f*}(v)

## Core Algorithms

| Method | Type | Convergence |
|--------|------|-------------|
| Gradient Descent + Backtracking | First-order | O(1/k) |
| Conjugate Gradient (Linear + Nonlinear) | First-order | O(√κ) |
| Nesterov Accelerated Gradient | First-order (optimal) | O(1/k²) |
| Heavy Ball (Polyak Momentum) | First-order | O(√κ) |
| Subgradient Method | Nonsmooth | O(1/√k) |
| Coordinate Descent | First-order | O(1/k) |
| Proximal Gradient (ISTA) | Composite | O(1/k) |
| Newton Method + Damped Newton | Second-order | Quadratic local |
| BFGS / DFP / SR1 | Quasi-Newton | Superlinear |
| L-BFGS | Limited-memory | R-linear |
| Trust-Region Newton | Second-order | Quadratic local |
| Gauss-Newton / Levenberg-Marquardt | NLS | Linear/Quadratic |
| Simplex Method (Two-Phase) | LP | Exponential worst, polynomial avg |
| Active-Set QP | QP | Finite termination |
| Quadratic Penalty / Log-Barrier | Constrained | Sublinear |

## Nine-School Curriculum Mapping

| School | Course | Key Topics |
|--------|--------|-----------|
| MIT | 6.255 Optimization | Gradient, Newton, duality, LP |
| Stanford | EE 364a Convex Opt. | Convexity, KKT, interior point |
| Berkeley | EE 227BT Convex Opt. | Fenchel duality, proximal |
| Michigan | IOE 510 Linear Prog. | Simplex, LP duality |
| Purdue | IE 535 Nonlinear Opt. | Newton, BFGS, trust-region |
| TU Munich | MA 3503 Nonlinear Opt. | Gradient, CG, quasi-Newton |
| ETH | 227-0427 Optimization | First/second-order, constrained |
| Tsinghua | Operations Research | LP, QP, simplex, duality |

## Build & Test

```
make all        # Build library + tests + examples
make test       # Build and run tests (9/9 pass)
make examples   # Build examples
make clean      # Clean build artifacts
```

## File Structure

```
mini-optimization-basics/
  Makefile
  README.md
  include/
    optimization_core.h       (L1: types, structures)
    gradient_methods.h        (L5: first-order methods)
    constrained_optimization.h (L4-L6: KKT, LP, QP)
    newton_methods.h          (L5/L8: Newton, quasi-Newton)
    convex_analysis.h         (L1-L4/L8: convexity, duality)
  src/
    optimization_core.c       (L1/L3/L5: core ops, Cholesky)
    gradient_methods.c        (L5/L8: GD, CG, Nesterov, prox)
    constrained_optimization.c (L4-L7: KKT, simplex, QP, penalty)
    newton_methods.c          (L5/L8: Newton, BFGS, L-BFGS, LM)
    convex_analysis.c         (L2/L4/L8: convexity, duality)
  tests/
    test_optimization.c       (9 tests, all pass)
  examples/
    example_gradient_descent.c (Rosenbrock function)
    example_linear_program.c   (Resource allocation)
    example_newton_method.c    (Newton/BFGS comparison)
  docs/
    knowledge-graph.md
    coverage-report.md
    gap-report.md
    course-alignment.md
    course-tree.md
```

## Module Status: COMPLETE

- L1-L6: Complete (all core definitions, theorems, methods, problems)
- L7: Complete (3 application examples with real-world context)
- L8: Partial (7/10 advanced topics implemented)
- L9: Partial (documented, not implemented)

Include/ + src/ total: 5249 lines
Tests: 9/9 passing
Examples: 3/3 building and running
