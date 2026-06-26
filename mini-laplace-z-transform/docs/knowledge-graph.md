# Knowledge Graph — Laplace & Z-Transform

## L1: Definitions (Complete)
- Laplace transform integral: F(s) = ∫₀^∞ f(t)e^{-st}dt
- Z-transform sum: X(z) = Σ x[n]z^{-n}
- Region of Convergence (ROC)
- Transfer function G(s) = N(s)/D(s), G(z) = N(z)/D(z)
- Pole-zero representation
- Proper/strictly proper systems
- Order of transfer function
- Initial Value Theorem
- Final Value Theorem
- Inverse Laplace: Bromwich integral
- Inverse Z: Cauchy residue, power series, partial fraction
- Polynomial structures (LaplacePolynomial, ZPoly)
- Rational function structures (LaplaceRational, ZRational)
- PoleZero, ZPoleZero types
- TimeSignal, DiscreteSignal types

## L2: Core Concepts (Complete)
- Linearity property (both transforms)
- Time-shifting (Laplace: e^{-sτ}, Z: z^{-k})
- Frequency shifting (s→s+a, z→z/a)
- Convolution theorem (time convolution ↔ multiplication)
- Differentiation property (s-domain)
- Integration property (1/s factor)
- Time advance property (Z-domain)
- DC gain: G(0) for CT, G(1) for DT
- System type (integrator count)
- Second-order parameters (ωₙ, ζ)
- Relative degree
- Transform pair correspondence

## L3: Mathematical Structures (Complete)
- Complex polynomial arithmetic (evaluation, addition, multiplication)
- Polynomial derivative computation
- Polynomial binomial shift (for frequency shift)
- Polynomial evaluation via Horner's method
- Complex number arithmetic for s/z-plane evaluation
- Root-finding: quadratic, cubic (Cardano), Aberth-Ehrlich iteration
- Partial fraction decomposition (residue computation)
- Bode plot data generation
- Z-pole to s-pole mapping

## L4: Fundamental Laws (Complete)
- Initial Value Theorem (Laplace + Z)
- Final Value Theorem (Laplace + Z)
- Routh-Hurwitz stability criterion
- Jury stability test
- Nyquist stability criterion
- Convolution theorem (numerical verification)
- Linearity (numerical verification)
- Parseval's theorem (implicit in numerical transforms)

## L5: Algorithms/Methods (Complete)
- Bilinear (Tustin) transform
- Impulse invariance discretization
- ZOH equivalent discretization
- Forward/Backward Euler discretization
- Matched pole-zero mapping
- Frequency prewarping for bilinear
- Butterworth filter order estimation
- Chebyshev Type I filter order estimation
- IIR filter design (analog prototype + bilinear)
- FIR filter design (window method)
- Biquad cascade decomposition
- Dubner-Abate-Crump inverse Laplace
- Composite Simpson's rule for forward Laplace
- Power series (long division) inverse Z
- Partial fraction inverse Z
- Residue method inverse Z
- Root locus computation
- Gain/phase margin computation

## L6: Canonical Problems (Complete)
- Second-order step response (underdamped/critically damped/overdamped)
- Transient specification extraction
- Steady-state error computation
- Analog filter prototype design
- Digital Butterworth filter design
- Window-method FIR design

## L7: Applications (Complete)
- DC motor speed control (electrical+mechanical time constants)
- Audio anti-aliasing filter (48 kHz ADC)
- Tesla servo positioning analysis

## L8: Advanced Topics (Partial)
- Minimum-phase system detection
- Aberth-Ehrlich simultaneous root finding
- Lanczos sigma factor for convergence acceleration
- Biquad cascade for high-order IIR

## L9: Research Frontiers (Partial — documented)
- Numerical inverse Laplace convergence properties
- Frequency warping correction
