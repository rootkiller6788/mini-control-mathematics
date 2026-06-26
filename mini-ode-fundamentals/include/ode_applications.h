/**
 * @file    ode_applications.h
 * @brief   Real-world ODE applications for control engineering.
 *
 * L7: Applications – DC motor, RLC circuit, spring-mass-damper,
 *     population models, vehicle dynamics.
 *
 * Reference: Ogata "Modern Control Engineering" (2010)
 *            Franklin, Powell, Emami-Naeini "Feedback Control" (2019)
 *            Dorf & Bishop "Modern Control Systems" (2016)
 *
 * Course Mapping:
 *   MIT 6.302 – DC motor control
 *   Stanford ENGR 105 – Feedback design applications
 *   Berkeley ME 132 – Electromechanical systems
 *   Caltech CDS 110 – Application examples
 *   Tsinghua 自动控制原理 – Engineering case studies
 */

#ifndef ODE_APPLICATIONS_H
#define ODE_APPLICATIONS_H

#include "ode_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ──────────────────────── L7: DC Motor ──────────────────────────────── */

/**
 * @brief RHS for DC motor: armature-controlled model.
 *
 * States: x₁ = θ (angle), x₂ = ω (angular velocity), x₃ = i (current)
 *   ẋ₁ = x₂
 *   ẋ₂ = -(b/J)·x₂ + (K/J)·x₃
 *   ẋ₃ = -(K/L)·x₂ - (R/L)·x₃ + V/L
 *
 * This is the standard model used in servo control design.
 * Reference: Ogata (2010), Chap. 4; Franklin (2019), Chap. 2
 *
 * @param t      Time
 * @param y      [θ, ω, i] (length 3)
 * @param f_out  [θ̇, ω̇, i̇] (length 3)
 * @param dim    3
 * @param ctx    Cast to DCMotorParams*
 * @return       0
 */
int ode_rhs_dc_motor(double t, const double *y, double *f_out,
                       int dim, void *ctx);

/**
 * @brief Simulate DC motor step response.
 *
 * @param params    Motor parameters
 * @param t_end     Simulation duration
 * @param n_steps   Number of time points
 * @param theta     Output: angle trajectory (length n_steps)
 * @param omega     Output: angular velocity trajectory (length n_steps)
 * @param current   Output: current trajectory (length n_steps)
 * @param t_vals    Output: time points (length n_steps)
 * @return          0 on success
 */
int ode_dc_motor_simulate(const DCMotorParams *params, double t_end,
                            int n_steps, double *theta, double *omega,
                            double *current, double *t_vals);

/**
 * @brief Find DC motor transfer function G(s) = ω(s)/V(s) in standard form:
 *        G(s) = K / ((Js+b)(Ls+R) + K²)
 *
 * @param params  Motor parameters
 * @param num     Output: numerator coefficients [b0] (length 1)
 * @param den     Output: denominator coefficients [a2, a1, a0] (length 3)
 */
void ode_dc_motor_transfer_function(const DCMotorParams *params,
                                      double *num, double *den);

/* ──────────────────────── L7: RLC Circuit ──────────────────────────── */

/**
 * @brief RHS for series RLC circuit with voltage source v(t).
 *
 * States: x₁ = v_c (capacitor voltage), x₂ = i_L (inductor current)
 *   C·v̇_c = i_L
 *   L·i̇_L = v(t) - R·i_L - v_c
 *
 * @param t      Time
 * @param y      [v_c, i_L] (length 2)
 * @param f_out  [v̇_c, i̇_L] (length 2)
 * @param dim    2
 * @param ctx    Cast to RLCCircuitParams*
 * @return       0
 */
int ode_rhs_rlc_circuit(double t, const double *y, double *f_out,
                          int dim, void *ctx);

/**
 * @brief Simulate RLC circuit transient response to a sinusoidal source.
 *
 * @param params    Circuit parameters
 * @param t_end     Simulation duration
 * @param n_steps   Number of time points
 * @param v_c       Output: capacitor voltage (length n_steps)
 * @param i_l       Output: inductor current (length n_steps)
 * @param t_vals    Output: time points (length n_steps)
 * @return          0 on success
 */
int ode_rlc_simulate(const RLCCircuitParams *params, double t_end,
                       int n_steps, double *v_c, double *i_l,
                       double *t_vals);

/**
 * @brief Compute RLC resonant frequency ω₀ = 1/√(LC) and damping factor
 *        ζ = (R/2)·√(C/L).
 *
 * Classification:
 *   ζ = 0  → undamped (purely oscillatory)
 *   0 < ζ < 1 → underdamped
 *   ζ = 1  → critically damped
 *   ζ > 1  → overdamped
 *
 * @param params  Circuit parameters
 * @param omega0  Output: natural frequency (rad/s)
 * @param zeta    Output: damping ratio
 * @return        ODE classification based on ζ
 */
int ode_rlc_damping_analysis(const RLCCircuitParams *params,
                               double *omega0, double *zeta);

/* ──────────────────────── L7: Spring-Mass-Damper ────────────────────── */

/**
 * @brief Simulate spring-mass-damper system: m·ẍ + c·ẋ + k·x = F·cos(ωt).
 *
 * @param params  Oscillator parameters
 * @param t_end   Simulation duration
 * @param n_steps Number of time points
 * @param x       Output: displacement (length n_steps)
 * @param v       Output: velocity (length n_steps)
 * @param t_vals  Output: time points
 * @return        0 on success
 */
int ode_spring_mass_simulate(const OscillatorParams *params, double t_end,
                               int n_steps, double *x, double *v,
                               double *t_vals);

/**
 * @brief Compute the frequency response (magnification factor) for
 *        a forced oscillator.
 *
 * M(ω) = 1 / √((1 - r²)² + (2ζr)²)  where r = ω/ωₙ, ωₙ = √(k/m)
 *
 * @param params  Oscillator parameters
 * @param omega   Array of driving frequencies
 * @param n_freqs Number of frequencies
 * @param M       Output: magnification factors
 * @param phase   Output: phase shifts (radians)
 */
void ode_oscillator_frequency_response(const OscillatorParams *params,
                                         const double *omega, int n_freqs,
                                         double *M, double *phase);

/**
 * @brief Find the resonant frequency (where magnification is maximized).
 *
 * ω_res = ωₙ·√(1 - 2ζ²) for ζ < 1/√2.
 * For ζ ≥ 1/√2, the response has no peak (monotonically decreasing).
 *
 * @param params  Oscillator parameters
 * @param omega_res Output: resonant frequency (0 if no resonance peak)
 * @param M_max   Output: max magnification factor
 * @return        true if resonance exists
 */
bool ode_find_resonance(const OscillatorParams *params,
                          double *omega_res, double *M_max);

/* ──────────────────────── L7: Population Dynamics ───────────────────── */

/**
 * @brief Simulate Lotka-Volterra system with given parameters.
 *
 * @param params  Lotka-Volterra parameters
 * @param x0      Initial prey population
 * @param y0      Initial predator population
 * @param t_end   Simulation time
 * @param n_steps Number of time points
 * @param prey    Output: prey trajectory (length n_steps)
 * @param pred    Output: predator trajectory (length n_steps)
 * @param t_vals  Output: time points
 * @return        0 on success
 */
int ode_predator_prey_simulate(const LotkaVolterraParams *params,
                                 double x0, double y0, double t_end,
                                 int n_steps, double *prey, double *pred,
                                 double *t_vals);

/**
 * @brief Find interior equilibrium of Lotka-Volterra:
 *        (x*, y*) = (γ/δ, α/β).
 *
 * This equilibrium is a center (neutrally stable), surrounded by
 * periodic orbits.
 *
 * @param params  Lotka-Volterra parameters
 * @param x_eq    Output: equilibrium prey population
 * @param y_eq    Output: equilibrium predator population
 */
void ode_lotka_volterra_equilibrium(const LotkaVolterraParams *params,
                                      double *x_eq, double *y_eq);

/**
 * @brief First integral (conserved quantity) of Lotka-Volterra:
 *        H(x, y) = δ·x - γ·ln(x) + β·y - α·ln(y)
 *
 * Level sets of H give periodic orbits.
 *
 * @param params  Parameters
 * @param x       Prey population
 * @param y       Predator population
 * @return        H(x, y) value
 */
double ode_lotka_volterra_integral(const LotkaVolterraParams *params,
                                     double x, double y);

/* ──────────────────────── L7: Lorenz System Analysis ────────────────── */

/**
 * @brief Simulate the Lorenz system.
 *
 * @param params  Lorenz parameters
 * @param x0, y0, z0 Initial conditions
 * @param t_end   Simulation duration
 * @param n_steps Number of time points
 * @param x, y, z Output: trajectory components
 * @param t_vals  Output: time points
 * @return        0 on success
 */
int ode_lorenz_simulate(const LorenzParams *params,
                          double x0, double y0, double z0,
                          double t_end, int n_steps,
                          double *x, double *y, double *z, double *t_vals);

/**
 * @brief Find equilibrium points of the Lorenz system:
 *        C₀ = (0, 0, 0)
 *        C± = (±√(β(ρ-1)), ±√(β(ρ-1)), ρ-1) for ρ > 1
 *
 * @param params      Lorenz parameters
 * @param eq_points   Output: equilibrium coordinates (9 doubles: 3 points × 3 dims)
 * @param n_eq        Output: number of equilibria (1 or 3)
 */
void ode_lorenz_equilibria(const LorenzParams *params,
                             double *eq_points, int *n_eq);

/* ──────────────────────── L7: Damped Pendulum ───────────────────────── */

/**
 * @brief Simulate the pendulum with damping and forcing.
 *
 * @param params  Pendulum parameters
 * @param theta0  Initial angle (rad)
 * @param omega0  Initial angular velocity (rad/s)
 * @param t_end   Simulation duration
 * @param n_steps Number of time points
 * @param theta   Output: angle trajectory
 * @param omega   Output: angular velocity trajectory
 * @param t_vals  Output: time points
 * @return        0 on success
 */
int ode_pendulum_simulate(const PendulumParams *params,
                            double theta0, double omega0,
                            double t_end, int n_steps,
                            double *theta, double *omega, double *t_vals);

/**
 * @brief Compute the small-angle period:
 *        T₀ = 2π·√(L/g)
 *
 * Large-angle period (elliptic integral):
 *        T(θ₀) = T₀·(2/π)·K(sin²(θ₀/2))
 * where K is the complete elliptic integral of the first kind.
 *
 * @param params    Pendulum parameters
 * @param theta0    Amplitude (rad)
 * @param T_small   Output: small-angle period
 * @param T_exact   Output: exact period (via elliptic integral approximation)
 */
void ode_pendulum_period(const PendulumParams *params, double theta0,
                           double *T_small, double *T_exact);

#ifdef __cplusplus
}
#endif

#endif /* ODE_APPLICATIONS_H */
