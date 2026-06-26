# Mini Control Mathematics

A collection of **from-scratch, zero-dependency C implementations** of engineering mathematics foundational to control theory. Each sub-module maps to MIT and Stanford (and other top-tier university) courses, covering complex analysis, ODEs, difference equations, transforms, linear algebra, optimization, numerical methods, and stochastic filtering — all built for the control engineer.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|------------|--------|-------------|
| [mini-complex-analysis](mini-complex-analysis/) | Analytic functions, Cauchy-Riemann, contour integration, residue calculus, series, Laplace/Fourier/Z transforms | MIT 18.04, Stanford MATH 116 |
| [mini-difference-equations](mini-difference-equations/) | Discrete-time dynamics, difference operators, Lyapunov stability, state-space, Z-transform | MIT 6.241, Stanford ENGR207B |
| [mini-laplace-z-transform](mini-laplace-z-transform/) | Laplace transform, Tustin bilinear mapping, Z-domain digital filter design (IIR/FIR) | MIT 6.003, Stanford EE 264 |
| [mini-linear-algebra-control](mini-linear-algebra-control/) | Controllability, observability, matrix factorizations, eigenvalue methods, state-space geometry | MIT 6.241, Stanford ENGR 205 |
| [mini-numerical-methods-control](mini-numerical-methods-control/) | Numerical linear algebra, ODE integration, eigenvalue computation, Riccati solvers, nonlinear methods, interpolation | MIT 18.335, Stanford CME 200 |
| [mini-ode-fundamentals](mini-ode-fundamentals/) | IVP theory, linear ODE systems, fundamental matrix, state transition matrix, control applications | MIT 18.03, Stanford CME 102 |
| [mini-optimization-basics](mini-optimization-basics/) | Convex analysis, gradient methods, Newton-type (BFGS/DFP), constrained optimization, simplex, interior point | MIT 6.7201, Stanford EE 364a |
| [mini-probability-stochastic](mini-probability-stochastic/) | Kalman filtering (linear, extended, unscented), Monte Carlo methods, Markov chains, Poisson processes | MIT 6.241, Stanford AA 273 |

## Design Philosophy

- **Zero external dependencies** — pure C99/C11, only standard library headers
- **Self-contained sub-modules** — each has its own `include/`, `src/`, `CMakeLists.txt`, and smoke tests
- **Theory-to-code mapping** — every module includes inline references to textbook sections and lecture notes
- **Control-engineering focus** — all mathematics is presented with control-theoretic motivation (Nyquist criterion, pole placement, LQR, state estimation)

## Building

Each sub-module is standalone. Build with CMake:

```bash
cd mini-complex-analysis
mkdir build && cd build
cmake ..
make
./smoke_test
```

Requires a **C99-compliant compiler** and **CMake ≥ 3.14**.

## Project Structure

```
0. mini-control-mathematics/
├── mini-complex-analysis/            # Analytic functions, contour integration, residue, transforms
├── mini-difference-equations/        # Discrete-time systems, stability, state-space, Z-transform
├── mini-laplace-z-transform/         # Laplace transform, Tustin mapping, digital filters
├── mini-linear-algebra-control/      # Controllability, observability, matrix decompositions
├── mini-numerical-methods-control/   # ODE integration, eigenvalue solvers, Riccati, nonlinear methods
├── mini-ode-fundamentals/            # IVP theory, linear ODE systems, state transition matrix
├── mini-optimization-basics/         # Convex optimization, gradient/Newton methods, constrained optimization
├── mini-probability-stochastic/      # Kalman filtering, Monte Carlo, Markov chains, stochastic processes
├── .gitignore
├── README.md
└── README-CN.md
```

## License

MIT
