# mini-difference-equations
Discrete-time dynamical systems, recurrence relations, and digital control fundamentals.

## Module Status: COMPLETE
- L1-L6: Complete
- L7: Complete (16 applications)
- L8: Complete (12 advanced methods)
- L9: Partial (documented, not implemented)
- include/ + src/ lines: 3000+

## Knowledge Coverage Summary
| Level | Name | Status |
|-------|------|--------|
| L1 | Definitions | Complete |
| L2 | Core Concepts | Complete |
| L3 | Math Structures | Complete |
| L4 | Fundamental Laws | Complete |
| L5 | Algorithms/Methods | Complete |
| L6 | Canonical Problems | Complete |
| L7 | Applications | Complete |
| L8 | Advanced Topics | Complete |
| L9 | Research Frontiers | Partial |

## Core Definitions (L1)
- Forward difference: Delta y[n] = y[n+1] - y[n]
- Backward difference: Nabla y[n] = y[n] - y[n-1]
- Shift operator: E y[n] = y[n+1]
- Linear difference equation: sum a_i y[n-i] = sum b_j x[n-j]
- Characteristic polynomial: P(lambda) = lambda^k + a_1 lambda^{k-1} + ... + a_k
- Impulse response h[n], step response s[n]
- z-transform: X(z) = sum x[n] z^{-n}
- Transfer function: H(z) = B(z)/A(z)
- State-space: x[k+1] = A x[k] + B u[k]
- Fixed point: x* where f(x*) = x*

## Core Theorems (L4)
- Jury Stability Criterion: |a_0| > |a_n|, P(1)>0, (-1)^n P(-1)>0, table constraints
- Final Value Theorem: lim x[n] = lim (z-1)X(z) as z approaches 1
- Initial Value Theorem: x[0] = lim X(z) as z approaches infinity
- Convolution Theorem: Z{x1*x2} = X1(z) X2(z)
- Superposition Principle: T{alpha x1 + beta x2} = alpha T{x1} + beta T{x2}
- Cayley-Hamilton (1x1): A^1 - A I = 0
- Bilinear Stability Preservation: Re(s)<0 iff |z|<1
- Existence and Uniqueness: k-th order + k IC implies unique solution
- BIBO Stability: sum |h[n]| < infinity implies BIBO stable

## Core Algorithms (L5)
1. Iterative solution - forward recursion O(k*N)
2. Characteristic root method - closed-form via Vandermonde O(k^3)
3. Undetermined coefficients - constant/exponential/sinusoid RHS
4. Direct z-transform computation
5. Inverse z-transform - partial fractions + long division
6. Jury stability table construction O(n^2)
7. Bilinear transform s<->z (Tustin)
8. Forward/Backward Euler discretization
9. ZOH discretization - exact step-invariant
10. State transition A^k via binary exponentiation O(log k * n^3)
11. Ackermann pole placement - discrete domain
12. Deadbeat controller/observer - all poles at z=0
13. FFT convolution - radix-2 Cooley-Tukey
14. Overlap-add/save convolution - block processing
15. LMS adaptive filter - stochastic gradient descent

## Engineering Problems (L6)
1. Fibonacci sequence (Binet formula)
2. Logistic map (population dynamics, chaos)
3. Digital FIR/IIR filter design
4. ARMA time series simulation
5. Economic cobweb model
6. Moving average filter
7. Loan amortization (compound interest)
8. Gambler ruin probability
9. Malthusian/Verhulst population growth
10. Forced harmonic oscillator
11. Sampled-data ZOH system
12. Lotka-Volterra predator-prey
13. SIR/SEIR epidemic models

## Applications (L7)
Discrete PID controller (Tustin + anti-windup), Kalman filter prediction,
exponential smoothing, digital resonator, DC motor discrete model,
quadcopter altitude, Mars rover slip detection, smart grid LFC,
climate carbon cycle, beer fermentation control, GPS CV model,
battery SoC estimation (Tesla/SpaceX), RC low-pass, mass-spring-damper,
Newton cooling, audio echo/reverb.

## Advanced Methods (L8)
Discrete Lyapunov equation, LQR via Riccati, LMS adaptive filter,
Kalman filter (predict+update), MPC, sliding mode control,
receding horizon, ILC, fuzzy gain scheduling, adaptive gain (Lyapunov),
Game of Life, Markov blanket, Bayesian state estimation, Monte Carlo.

## Course Mapping
| School | Course | Topics |
|--------|--------|--------|
| MIT | 6.241 Dynamic Systems | Discrete SS, stability |
| MIT | 6.003 Signals and Systems | z-transform, convolution |
| Stanford | ENGR207B Linear Control | Controllability |
| Stanford | EE 264 DSP | FIR/IIR, FFT |
| Berkeley | ME 132 Dynamic Systems | Diff operators |
| Caltech | CDS 110 Linear Systems | Jury criterion |
| ETH | 151-0591 Control I | Discretization |
| Cambridge | 3F2 Systems and Control | Discrete stability |
| Tsinghua | Auto Control Theory | Diff equations |

## Build and Test
```
make          # Build library + examples + test
make test     # Run 21 assert-based tests
make examples # Build examples only
make count    # Line count
```

## References
- Elaydi, "An Introduction to Difference Equations" (2005)
- Jury, "Theory and Application of the z-Transform Method" (1964)
- Ogata, "Discrete-Time Control Systems" (1995)
- Astrom and Wittenmark, "Computer-Controlled Systems" (1997)
- Oppenheim and Schafer, "Discrete-Time Signal Processing" (2010)
