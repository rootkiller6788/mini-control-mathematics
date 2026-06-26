/**
 * diffeq.h — Core Difference Equation Types and Definitions
 *
 * Domain: Discrete-time dynamical systems, recurrence relations,
 *         digital control, signal processing fundamentals.
 *
 * Knowledge Level Coverage:
 *   L1: Core definitions — difference operators, equation classes, equilibrium
 *   L2: Core concepts — superposition, impulse/step response, linearity
 *   L3: Engineering quantities — time constant, settling time, damping ratio
 *   L4: Fundamental laws — existence/uniqueness, superposition principle
 *
 * Reference: Elaydi, S. "An Introduction to Difference Equations" (3rd ed, 2005)
 * Reference: Ogata, K. "Discrete-Time Control Systems" (2nd ed, 1995)
 * Reference: Jury, E.I. "Theory and Application of the z-Transform Method" (1964)
 *
 * Course Mapping:
 *   MIT 6.241  — Dynamic Systems and Control (discrete-time systems)
 *   Stanford ENGR207B — Linear Control Systems (discrete domain)
 *   ETH 151-0591 — Control Systems I (difference equations)
 */

#ifndef DIFFEQ_H
#define DIFFEQ_H

#include <stddef.h> /* size_t */
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── L1: Core Definitions ────────────────────────────────────────────── */

/**
 * Difference operator types.
 * Δ (forward):  Δy[n] = y[n+1] - y[n]
 * ∇ (backward): ∇y[n] = y[n] - y[n-1]
 * E (shift):    Ey[n] = y[n+1]
 * E^{-1} (lag): E^{-1}y[n] = y[n-1]
 */
typedef enum {
    DIFF_OP_FORWARD,   /* Δ operator */
    DIFF_OP_BACKWARD,  /* ∇ operator */
    DIFF_OP_SHIFT,     /* E operator (advance) */
    DIFF_OP_LAG        /* E^{-1} operator (delay) */
} DiffOperator;

/**
 * Difference equation classification.
 */
typedef enum {
    DEQ_TYPE_LINEAR_CONSTANT,     /* Linear with constant coefficients */
    DEQ_TYPE_LINEAR_VARIABLE,     /* Linear with variable coefficients */
    DEQ_TYPE_NONLINEAR,           /* Nonlinear difference equation */
    DEQ_TYPE_HOMOGENEOUS,         /* Homogeneous (RHS = 0) */
    DEQ_TYPE_NONHOMOGENEOUS       /* Nonhomogeneous (RHS ≠ 0) */
} DiffEqType;

/**
 * Linear difference equation representation:
 *   a_0 y[n] + a_1 y[n-1] + ... + a_k y[n-k] =
 *   b_0 x[n] + b_1 x[n-1] + ... + b_m x[n-m]
 *
 * With leading coefficient a_0 normalized to 1:
 *   y[n] + a_1 y[n-1] + ... + a_k y[n-k] =
 *   b_0 x[n] + b_1 x[n-1] + ... + b_m x[n-m]
 */
typedef struct {
    int      order;         /* k = max order of y terms */
    double  *a_coeffs;      /* [a_1, ..., a_k], length k (a_0 = 1) */
    int      input_order;   /* m = max order of x terms */
    double  *b_coeffs;      /* [b_0, ..., b_m], length m+1 */
    DiffEqType type;        /* Equation classification */
} LinearDiffEq;

/**
 * Discrete-time signal: y[n] defined on integer indices.
 */
typedef struct {
    double  *values;        /* Signal values */
    int      length;        /* Number of samples */
    int      start_index;   /* Index of first element (typically 0) */
} DiscreteSignal;

/**
 * Impulse response h[n] of a linear time-invariant (LTI) discrete system.
 * Fully characterizes the system via convolution.
 */
typedef struct {
    double  *h;             /* Impulse response values h[0], h[1], ..., h[N-1] */
    int      length;        /* Number of impulse response samples */
    int      is_iir;        /* 1 if infinite impulse response, 0 if FIR */
} ImpulseResponse;

/**
 * Step response s[n] = sum_{k=0}^n h[k].
 */
typedef struct {
    double  *s;             /* Step response values */
    int      length;        /* Number of samples */
    double   steady_state;  /* Steady-state gain (if stable) */
} StepResponse;

/**
 * Characteristic polynomial roots for linear difference equation:
 *   P(λ) = λ^k + a_1 λ^{k-1} + ... + a_k = 0
 *
 * The general homogeneous solution is:
 *   y_h[n] = c_1 λ_1^n + c_2 λ_2^n + ... + c_k λ_k^n
 * (for distinct roots; modifications for repeated/complex roots)
 */
typedef struct {
    int             num_roots;      /* Number of characteristic roots (= order) */
    double complex *roots;          /* Complex roots array */
    int            *multiplicity;   /* Multiplicity of each root */
} CharPoly;

/**
 * Equilibrium / fixed point classification.
 */
typedef enum {
    FIXED_STABLE,       /* |f'(x*)| < 1 — attracting */
    FIXED_UNSTABLE,     /* |f'(x*)| > 1 — repelling */
    FIXED_NEUTRAL,      /* |f'(x*)| = 1 — neutral/indifferent */
    FIXED_UNKNOWN       /* Cannot classify */
} FixedPointType;

/**
 * Fixed point: x* such that f(x*) = x*.
 */
typedef struct {
    double         value;    /* The fixed point value x* */
    FixedPointType type;     /* Stability classification */
    double         multiplier; /* |f'(x*)| — the multiplier */
} FixedPoint;

/* ── L3: Engineering Quantities ──────────────────────────────────────── */

/**
 * Discrete-time engineering response quantities extracted from
 * step/impulse response data.
 */
typedef struct {
    double   time_constant;       /* Time to reach ~63.2% of steady state (samples) */
    double   rise_time;           /* 10%-90% rise time (samples) */
    double   settling_time;       /* Time to stay within ±2% (samples) */
    double   overshoot_percent;   /* Percentage overshoot (%) */
    double   peak_time;           /* Time to first peak (samples) */
    double   steady_state_error;  /* Steady-state error */
    double   damping_ratio;       /* Equivalent continuous ζ */
    double   natural_frequency;   /* Equivalent continuous ω_n (rad/sample) */
} DiscreteResponseMetrics;

/* ── L4: Existence & Uniqueness ──────────────────────────────────────── */

/**
 * Result of verifying existence and uniqueness conditions for
 * a k-th order difference equation with k initial conditions.
 */
typedef struct {
    int     exists;          /* 1 if a solution exists */
    int     is_unique;       /* 1 if the solution is unique */
    int     order;           /* Order of the equation */
    int     num_conditions;  /* Number of initial conditions provided */
    char    message[256];    /* Human-readable explanation */
} ExistenceResult;

/* ── API: L1 — Difference Operators ──────────────────────────────────── */

/**
 * Apply forward difference operator: Δy[n] = y[n+1] - y[n]
 * Complexity: O(1) per call, O(N) for signal
 */
double  forward_difference(const DiscreteSignal *sig, int n);
void    forward_difference_signal(const DiscreteSignal *in, DiscreteSignal *out);

/**
 * Apply backward difference operator: ∇y[n] = y[n] - y[n-1]
 * Complexity: O(1) per call, O(N) for signal
 */
double  backward_difference(const DiscreteSignal *sig, int n);
void    backward_difference_signal(const DiscreteSignal *in, DiscreteSignal *out);

/**
 * Apply shift operator: Ey[n] = y[n+1], E^{-1}y[n] = y[n-1]
 * Complexity: O(1) per call, O(N) for signal
 */
double  shift_operator(const DiscreteSignal *sig, int n, int k);
void    shift_operator_signal(const DiscreteSignal *in, DiscreteSignal *out, int k);

/* ── API: L1/L2 — Equation Classification & Analysis ─────────────────── */

/**
 * Classify a difference equation based on its coefficients.
 * Complexity: O(k) where k = order
 */
DiffEqType classify_difference_eq(int order, const double *a, const double *b);

/**
 * Construct the characteristic polynomial P(λ) = λ^k + a_1 λ^{k-1} + ... + a_k
 * Complexity: O(k)
 */
void    build_characteristic_polynomial(const LinearDiffEq *eq, CharPoly *cp);

/**
 * Find fixed points x* of 1D map f(x) satisfying f(x*) = x*.
 * Uses bisection method to find roots of g(x) = f(x) - x.
 * Complexity: O(num_guesses * log(bracketing_range / tol))
 */
int     find_fixed_points(double (*f)(double),
                          double x_min, double x_max, double tol,
                          FixedPoint *fixed_points, int max_points);

/**
 * Classify a fixed point's stability via |f'(x*)|.
 * Complexity: O(1)
 */
FixedPointType classify_fixed_point(double multiplier);

/* ── API: L5 — Solution Methods ─────────────────────────────────────── */

int     solve_iterative_lindiffeq(const LinearDiffEq *eq,
                                   const double *x,
                                   const double *initial_y,
                                   int num_samples,
                                   double *y);
int     solve_characteristic_roots(const LinearDiffEq *eq,
                                    const double *initial_y,
                                    int num_samples,
                                    double *y);
int     solve_undetermined_coeffs(const LinearDiffEq *eq,
                                   int rhs_type,
                                   double rhs_param1,
                                   double rhs_param2,
                                   const double *initial_y,
                                   int num_samples,
                                   double *y);

/* ── API: L2 — Convolution ──────────────────────────────────────────── */

int     convolve_discrete(const double *x, int Nx,
                           const double *h, int Nh,
                           double *y);
int     convolve_fft(const double *x, int Nx,
                      const double *h, int Nh,
                      double *y);
int     deconvolve_discrete(const double *y, int Ny,
                             const double *h, int Nh,
                             double *x);
void    auto_correlation_discrete(const double *x, int N,
                                   double *Rxx, int max_lag);
void    cross_correlation_discrete(const double *x, const double *y,
                                    int N, double *Rxy, int max_lag);
int     circular_convolve_discrete(const double *x, const double *h,
                                    int N, double *y);
int     convolution_properties_verify(const double *x,
                                       const double *h1, const double *h2,
                                       int Nx, int Nh, double tol);
void    toeplitz_from_impulse_response(const double *h, int Nh,
                                        int N, double **T);
int     overlap_add_convolution(const double *x, int Nx,
                                 const double *h, int Nh,
                                 int block_size, double *y);
int     overlap_save_convolution(const double *x, int Nx,
                                  const double *h, int Nh,
                                  int block_size, double *y);
/* ── API: L2 — Impulse & Step Response ───────────────────────────────── */

int     compute_impulse_response(const LinearDiffEq *eq,
                                 int num_samples,
                                 ImpulseResponse *ir);
int     compute_step_response(const ImpulseResponse *ir,
                              StepResponse *sr);
int     compute_step_response_direct(const LinearDiffEq *eq,
                                     int num_samples,
                                     StepResponse *sr);

/* ── API: L4 — Superposition & Linearity ─────────────────────────────── */

int     verify_superposition(const LinearDiffEq *eq,
                             const DiscreteSignal *x1,
                             const DiscreteSignal *x2,
                             double alpha, double beta,
                             int num_samples, double tol);
int     verify_linearity(const LinearDiffEq *eq, int num_samples, double tol);

/* ── API: L4 — Existence & Uniqueness ────────────────────────────────── */

ExistenceResult check_existence_uniqueness(const LinearDiffEq *eq,
                                           const double *initial_conditions,
                                           int num_conditions);

/* ── API: L3 — Engineering Response Metrics ──────────────────────────── */

void    compute_response_metrics(const StepResponse *sr,
                                 double sampling_period,
                                 DiscreteResponseMetrics *metrics);
double  time_constant_discrete(const CharPoly *cp);
double  settling_time_discrete(const CharPoly *cp);
double  damping_ratio_discrete(const CharPoly *cp);
double  natural_frequency_discrete(const CharPoly *cp);
double  overshoot_percent_discrete(const CharPoly *cp);

/* ── Utility ─────────────────────────────────────────────────────────── */

void    linear_diffeq_free(LinearDiffEq *eq);
void    discrete_signal_free(DiscreteSignal *sig);
void    char_poly_free(CharPoly *cp);
void    impulse_response_free(ImpulseResponse *ir);
void    step_response_free(StepResponse *sr);
DiscreteSignal* discrete_signal_alloc(int length, int start_index);
LinearDiffEq*   linear_diffeq_alloc(int order, int input_order);
int linear_diffeq_set_coeffs(LinearDiffEq *eq,
                             const double *a, const double *b);
void linear_diffeq_print(const LinearDiffEq *eq);
double nonlinear_map_eval(double (*f)(double), double y_n);
int nonlinear_map_iterate(double (*f)(double),
                          double y0, int N,
                          double *trajectory);

#ifdef __cplusplus
}
#endif

#endif /* DIFFEQ_H */

/* ── Applications & Advanced (L6-L8) ─────────────────────────────────── */

/* Fibonacci */
void fibonacci_closed_form(int n, double *Fn);
int fibonacci_iterative(int n, double *values, int count);

/* Population dynamics */
int logistic_map_iterate(double r, double x0, int N, double *trajectory);
void logistic_bifurcation_data(double r_min, double r_max, int r_steps,
    int transient, int steady, double *results);
void malthusian_growth(double P0, double r, int N, double *population);
void verhulst_logistic_discrete(double P0, double r, double K, int N, double *population);

/* ARMA / Time series */
int arma_model_simulate(const double *ar_coeffs, int p,
    const double *ma_coeffs, int q, double noise_std, int N, double *y);
int moving_average_filter(const double *x, int N, int window, double *y);

/* Cobweb / Economic */
int cobweb_model(double a, double b, double initial_price, int N,
    double *prices, double *quantities);

/* Gambler's ruin */
double gamblers_ruin_probability(double p, int initial, int N);

/* Engineering applications */
typedef struct {double Kp,Ki,Kd,Ts;double integral,prev_error,prev_output;
 double out_min,out_max;int has_aw;} DiscretePID;
void discrete_pid_init(DiscretePID *pid, double Kp, double Ki, double Kd,
    double Ts, double out_min, double out_max);
double discrete_pid_update(DiscretePID *pid, double setpoint, double measurement);

double exponential_smoothing(double alpha, double x, double *prev_y);
double digital_resonator(double r, double omega0, double x, double *y1, double *y2);

typedef struct {double J,B,Kt,Ts;double prev_w,prev_i;} DCMotorDisc;
void dc_motor_discrete_init(DCMotorDisc *motor, double J, double B, double Kt, double Ts);
double dc_motor_discrete_update(DCMotorDisc *motor, double voltage, double R);

/* RLC / Mechanical */
typedef struct {double R,C,Ts;double Vc_prev;} RCLowPass;
void rc_lowpass_init(RCLowPass *rc, double R, double C, double Ts);
double rc_lowpass_update(RCLowPass *rc, double Vin);

typedef struct {double m,c,k,Ts;double x_prev,v_prev;} MassSpringDamper;
void msd_discrete_init(MassSpringDamper *msd, double m, double c, double k, double Ts);
void msd_discrete_update(MassSpringDamper *msd, double F, double *x, double *v);

typedef struct {double k,T_amb,Ts;double T_prev;} NewtonCooling;
void newton_cooling_init(NewtonCooling *nc, double k, double T_amb, double Ts);
double newton_cooling_update(NewtonCooling *nc, double T_current);

/* Quadcopter */
typedef struct {double mass,g,k_thrust,Ts;double z_prev,v_prev;} QuadAltDisc;
void quad_alt_disc_init(QuadAltDisc *q, double mass, double k_thrust, double Ts);
void quad_alt_disc_update(QuadAltDisc *q, double throttle, double *z, double *v);

/* Sampled-data ZOH */
typedef struct {double A,B,C;double Ts;double x_prev;} SampledDataZOH;
void sampled_data_zoh_init(SampledDataZOH *sys, double A, double B, double C, double Ts);
double sampled_data_zoh_update(SampledDataZOH *sys, double u);

/* Lotka-Volterra */
void lotka_volterra_discrete(double alpha, double beta, double gamma, double delta,
    double *prey, double *predator, int N);

/* SEIR */
typedef struct {double S,E,Inf,R;double beta,sigma,gamma;int N;} SEIRModel;
void seir_model_init(SEIRModel *m, double S0, double E0, double I0, double R0,
    double beta, double sigma, double gamma);
void seir_model_step(SEIRModel *m);

/* Forced oscillator */
void forced_oscillator_discrete(double m, double c, double k, const double *F,
    int N, double Ts, double *displacement, double *velocity);

/* Advanced methods (L8) */
typedef struct {int N;double mu;double *w;double *xbuf;int bi;} LMSFilter;
void lms_filter_init(LMSFilter *lms, int order, double mu);
double lms_filter_update(LMSFilter *lms, double x, double desired);
void lms_filter_free(LMSFilter *lms);

int discrete_lyapunov_solve(const double *A, int n, const double *const *Q, double **P);

typedef struct {int rows,cols,**grid,**next;} GameOfLife;
void gol_init(GameOfLife *gol, int rows, int cols);
void gol_step(GameOfLife *gol);
void gol_free(GameOfLife *gol);
