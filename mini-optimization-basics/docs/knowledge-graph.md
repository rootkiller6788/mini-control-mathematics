# Knowledge Graph — mini-optimization-basics

## L1: Definitions
- convex_set_t: Convex set structure with membership/projection/support oracles
- convex_function_t: Convex function with eval/grad/hess/conjugate/prox oracles  
- vector_t, matrix_t: Core linear algebra types
- gradient_t, hessian_t: Gradient/Hessian with derived quantities
- optimality_cert_t: Optimality certificate (first-order + second-order)
- linear_program_t: Standard-form LP structure
- quadratic_program_t: Convex QP structure
- general_nlp_t: General nonlinear program with KKT fields

## L2: Core Concepts
- Convexity: Set convexity (Jensen), function convexity (first-order condition)
- Duality: Weak/strong duality, Lagrangian duality, Fenchel duality
- Optimality conditions: Fermat theorem, first-order necessary, second-order sufficient
- Constraint qualifications: Slater condition
- Operations preserving convexity: Nonnegative sum, pointwise max, affine composition
- Convex conjugate (Fenchel transform), perspective function, infimal convolution
- Moreau envelope: Smooth approximation preserving minima

## L3: Engineering Quantities
- Convergence rates: O(1/k) sublinear, O((1-μ/L)^k) linear, quadratic
- Condition number κ = λ_max/λ_min
- Lipschitz constant L of gradient
- Newton decrement λ(x): affine-invariant optimality measure
- Duality gap: optimality certificate for primal-dual methods

## L4: Conservation Laws / Theorems
- Fermat theorem: ∇f(x*) = 0 at local minimum
- KKT Theorem: Necessary + sufficient for convex NLP with Slater
- Weak/Strong duality theorems
- Separating hyperplane theorem
- Moreau decomposition: v = prox_f(v) + prox_{f*}(v)
- Fenchel duality theorem
- Fundamental theorem of LP: optimum at vertex

## L5: Engineering Methods
- Gradient descent with backtracking line search (Armijo)
- Linear/nonlinear conjugate gradient (Hestenes-Stiefel, Fletcher-Reeves)
- Newton method, damped Newton (self-concordant)
- BFGS, DFP, SR1 quasi-Newton methods
- L-BFGS limited-memory BFGS
- Trust-region Newton (Levenberg-Marquardt)
- Gauss-Newton and Levenberg-Marquardt for NLS
- Proximal gradient (ISTA), subgradient method
- Coordinate descent, heavy ball (Polyak momentum)
- Nesterov accelerated gradient (optimal O(1/k²))
- Simplex method (two-phase), active-set QP
- Penalty and barrier methods
- Projected gradient, alternating projections (POCS)

## L6: Engineering Problems
- Linear programming: Resource allocation, transportation
- Quadratic programming: Portfolio optimization
- Unconstrained minimization: Parameter estimation
- Composite optimization: f(x) + g(x) structure
- Nonlinear least squares: Gauss-Newton, LM

## L7: Applications
- Factory resource allocation (simplex LP)
- Rosenbrock function minimization (benchmark)
- Quadratic model fitting (Newton/BFGS comparison)
- Transportation problem (resource allocation LP)

## L8: Advanced Methods
- L-BFGS for large-scale problems
- Trust-region methods for robust globalization
- Fenchel duality framework
- Moreau envelope/decomposition
- Proximal gradient (ISTA) for composite problems
- Nesterov acceleration (optimal first-order rates)
- Subdifferential calculus

## L9: Research Frontiers
- Stochastic gradient methods
- Distributed optimization / ADMM
- Non-convex optimization landscapes
- Self-concordant barrier theory
