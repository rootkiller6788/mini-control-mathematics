# Course Tree — mini-optimization-basics

## Prerequisites
```
Linear Algebra (vector, matrix, eigenvalues)
  -> Calculus (gradients, Hessians, Taylor expansion)
       -> Real Analysis (convergence, continuity, Lipschitz)
            -> mini-optimization-basics
```

## Internal Dependencies
```
optimization_core.h (vector_t, matrix_t, convex_set_t, convex_function_t)
  +-- gradient_methods.h (gradient descent, CG, line search)
  |     +-- newton_methods.h (Newton, BFGS, L-BFGS, trust-region)
  +-- constrained_optimization.h (KKT, Lagrange, LP, QP, penalty)
  +-- convex_analysis.h (convexity, duality, projections, Moreau)
```

## Postrequisites
```
mini-optimization-basics
  +-- mini-optimal-control (LQR, MPC optimization)
  +-- mini-estimation-theory (MLE, MAP via optimization)
  +-- mini-system-identification (parameter estimation)
  +-- mini-machine-learning-basics (empirical risk minimization)
```
