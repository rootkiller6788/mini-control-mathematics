# Knowledge Graph — mini-difference-equations

## L1: Definitions — Complete
- Forward difference operator Δ
- Backward difference operator ∇
- Shift operator E
- Linear difference equation (constant/variable coeffs)
- Homogeneous/nonhomogeneous classification
- Characteristic polynomial
- Impulse response h[n]
- Step response s[n]
- Kronecker delta δ[n], unit step u[n]
- Fixed point / equilibrium
- z-transform definition X(z) = Σ x[n] z^{-n}
- Region of convergence (ROC)
- Transfer function H(z)
- State-space quadruple (A,B,C,D)
- Discrete signal types (impulse, step, ramp, exponential, sinusoid)

## L2: Core Concepts — Complete
- Superposition principle
- Linearity (additivity + homogeneity)
- Convolution sum
- Impulse response characterization
- BIBO stability
- Asymptotic stability
- Controllability
- Observability
- Causality
- Time-invariance
- Minimum phase

## L3: Mathematical Structures — Complete
- z-transform pairs
- Pole-zero representation
- Frequency response H(e^{jω})
- Companion matrix
- Controllable canonical form
- Observable canonical form
- Diagonal/Jordan form
- State transition matrix A^k
- Toeplitz matrix structure
- Discrete time constant τ_d
- Settling time, rise time, overshoot
- Damping ratio ζ, natural frequency ω_n
- Gain margin, phase margin (discrete)
- Signal energy, power, RMS, peak

## L4: Fundamental Laws — Complete
- Jury stability criterion (necessary & sufficient)
- Schur-Cohn stability test
- Final value theorem (z-domain)
- Initial value theorem (z-domain)
- Convolution theorem (z-domain)
- Time-shift property
- Linearity of z-transform
- Cayley-Hamilton theorem (discrete)
- Existence and uniqueness theorem
- Bilinear transform stability preservation
- Superposition principle (verified numerically)
- Convolution properties (commutative, associative, distributive)

## L5: Computational Methods — Complete
- Iterative (recursive) solution
- Characteristic root method
- Method of undetermined coefficients
- Direct z-transform computation
- Inverse z-transform (partial fractions)
- Inverse z-transform (long division)
- Partial fraction expansion
- Jury table construction
- Bilinear transform s↔z
- Forward/Backward Euler discretization
- Tustin (bilinear) discretization
- ZOH discretization
- Matched pole-zero discretization
- State transition via binary exponentiation
- Ackermann pole placement
- Deadbeat controller/observer
- FFT-based convolution
- Overlap-add/save convolution

## L6: Canonical Problems — Complete
- Fibonacci sequence (Binet formula)
- Logistic map / population dynamics
- Digital FIR/IIR filter design
- ARMA time series simulation
- Economic cobweb model
- Moving average filter
- Compound interest / loan amortization
- Gambler's ruin probability
- Malthusian / Verhulst population growth
- Forced harmonic oscillator
- Sampled-data ZOH system
- Lotka-Volterra predator-prey
- SIR/SEIR epidemic models

## L7: Applications — Complete (12 applications)
- Discrete PID controller (Tustin + anti-windup)
- Kalman filter prediction step
- Exponential smoothing
- Digital resonator
- DC motor discrete model
- Quadcopter altitude model
- Mars rover slip detection
- Smart grid LFC (load frequency control)
- Climate carbon cycle model
- Beer fermentation temperature control
- GPS constant-velocity model
- Battery SoC estimation (Coulomb counting)
- RC low-pass filter discretized
- Mass-spring-damper discrete model
- Newton's law of cooling
- Audio echo/reverb (comb filter, Schroeder)

## L8: Advanced Topics — Complete (12 methods)
- Discrete Lyapunov equation solver
- Discrete-time LQR (Linear Quadratic Regulator)
- LMS adaptive filter
- Discrete Kalman filter (prediction+update)
- Model Predictive Control (MPC)
- Discrete sliding mode control
- Receding horizon control
- Iterative learning control (ILC)
- Fuzzy gain scheduling
- Adaptive gain via Lyapunov method
- Game of Life (cellular automaton)
- Markov blanket identification
- Bayesian state estimation
- Monte Carlo random walk

## L9: Research Frontiers — Partial
- Fractional difference equations (documented)
- Neural network layers as difference equations (documented)
- Quantum difference equations / walks (documented)
- Reservoir computing / echo state networks (documented)
