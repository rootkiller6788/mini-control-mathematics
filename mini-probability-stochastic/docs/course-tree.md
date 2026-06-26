# Course Tree — mini-probability-stochastic

## Prerequisite Dependency Graph

```
                    Calculus I-II
                         |
              +----------+----------+
              |          |          |
         Linear Algebra  |    Differential Eqns
              |          |          |
              +-------+  |  +-------+
                      |  |  |
                   Probability Theory
                   (core module)
                         |
         +-------+-------+--------+--------+
         |       |        |        |        |
    Statistics  Stochastic  Monte  Estimation
    (L1-L6)    Processes   Carlo   Filtering
               (L1-L2,L6)  (L5,L8) (L5-L7)
         |       |        |        |
         +-------+--------+--------+
                         |
              Advanced Topics (L8-L9)
              - MCMC diagnostics
              - QMC methods
              - UKF / Particle filters
              - HMC, SMC (documented)
```

## Topic Dependency Tree

### Level 1: Foundations (no prerequisites beyond calculus)
- Set theory, sample spaces, events
- Random variables, distribution functions
- Probability mass/density functions

### Level 2: Core Probability (requires L1)
- Expectation, moments
- Conditional probability, independence
- Bayes theorem
- Common distributions (Gaussian, exponential, uniform, etc.)

### Level 3: Statistical Summaries (requires L2)
- Mean, variance, standard deviation
- Skewness, kurtosis
- Covariance, correlation
- Covariance matrix (requires linear algebra)

### Level 4: Limit Theorems (requires L2-L3)
- Law of large numbers
- Central limit theorem
- Chapman-Kolmogorov equation
- Gibbs inequality

### Level 5: Estimation & Methods (requires L3-L4)
- Maximum likelihood estimation
- KS goodness-of-fit test
- Kalman filter (requires linear algebra + L2)
- Monte Carlo integration
- MCMC (Metropolis-Hastings, Gibbs)
- Entropy, information divergence
- Online statistics (Welford algorithm)

### Level 6: Engineering Problems (requires L5)
- Distribution fitting
- State estimation (Kalman filter application)
- Absorption analysis for Markov chains
- Time-series analysis (ACF, PSD)
- Bootstrap confidence intervals

### Level 7: Applications (requires L6)
- Radar/sonar target tracking
- GPS/IMU sensor fusion
- Weather risk modeling
- Quality control
- Financial volatility modeling

### Level 8: Advanced Methods (requires L5-L7)
- MCMC convergence diagnostics (Gelman-Rubin R-hat)
- Effective sample size computation
- Quasi-Monte Carlo (Halton, Sobol sequences)
- Unscented Kalman filter
- Particle filters (SIR)

### Level 9: Research Frontiers (requires L8)
- Hamiltonian Monte Carlo
- Sequential Monte Carlo (adaptive)
- Gaussian process regression
- Fractional Brownian motion
- Causal inference

## Cross-Module Dependencies

| This Module Depends On | Dependency Type |
|------------------------|-----------------|
| Linear algebra (matrix operations) | `mat_mul_n()`, `mat_cholesky_solve()`, `mat_inverse()` |
| Numerical methods (special functions) | `prob_erf_approx()`, `prob_lgamma()`, incomplete gamma/beta |
| C standard library | `math.h` (sin, cos, exp, log, sqrt, pow, fabs) |

| Modules That Depend On This | Use Case |
|------------------------------|----------|
| mini-control-theory | Kalman filter for state estimation, stochastic control |
| mini-machine-learning | MLE for parameter estimation, Bayesian inference |
| mini-signal-processing | Particle filter, Kalman smoother |
| mini-financial-engineering | GBM, Monte Carlo for option pricing |
