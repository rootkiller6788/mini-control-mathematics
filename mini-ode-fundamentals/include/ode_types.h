/**
 * @file    ode_types.h
 * @brief   Core type definitions for ODE systems in control mathematics.
 *
 * L1 Definitions: All fundamental ODE types, structs, and enumerations.
 *
 * Reference: Boyce & DiPrima "Elementary Differential Equations" (2017)
 *            Hirsch & Smale "Differential Equations, Dynamical Systems" (2013)
 *
 * Course Mapping:
 *   MIT 6.241 Dynamic Systems – State-space representation
 *   Stanford ENGR 207B – Linear system structures
 *   Berkeley ME 132 – ODE classification
 *   ETH 151-0591 – Vector field definitions
 *   Tsinghua 自动控制原理 – System modeling types
 */

#ifndef ODE_TYPES_H
#define ODE_TYPES_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ──────────────────────────── L1: Core Enumerations ──────────────────── */

/** ODE classification by linearity and structure */
typedef enum {
    ODE_LINEAR_CONSTANT_COEFF,    /* y' = Ay + b(t), A constant */
    ODE_LINEAR_VARIABLE_COEFF,    /* y' = A(t)y + b(t) */
    ODE_NONLINEAR_AUTONOMOUS,     /* y' = f(y) */
    ODE_NONLINEAR_NONAUTONOMOUS   /* y' = f(t, y) */
} ODEClass;

/** ODE order classification */
typedef enum {
    ODE_FIRST_ORDER,              /* dy/dt = f(t,y) */
    ODE_SECOND_ORDER,             /* d²y/dt² = f(t,y,y') */
    ODE_HIGHER_ORDER,             /* n-th order ODE */
    ODE_FIRST_ORDER_SYSTEM        /* vector ODE: y' = f(t,y) */
} ODEOrder;

/** Stability classification for equilibrium points (L1: Definitions, L2: Concepts) */
typedef enum {
    STABILITY_ASYMPTOTICALLY_STABLE,   /* All eigenvalues Re(λ) < 0 */
    STABILITY_MARGINALLY_STABLE,       /* Some Re(λ) = 0, none > 0 (center) */
    STABILITY_UNSTABLE,                /* At least one Re(λ) > 0 */
    STABILITY_SADDLE,                  /* Mixed signs, 2D saddle */
    STABILITY_UNKNOWN                  /* Cannot determine */
} StabilityType;

/** Classification of 2D linear equilibria (L2: Core Concepts) */
typedef enum {
    EQUILIBRIUM_STABLE_NODE,           /* λ₁,λ₂ < 0 real */
    EQUILIBRIUM_UNSTABLE_NODE,         /* λ₁,λ₂ > 0 real */
    EQUILIBRIUM_SADDLE_POINT,          /* λ₁ < 0 < λ₂ real */
    EQUILIBRIUM_STABLE_FOCUS,          /* Re(λ) < 0, complex */
    EQUILIBRIUM_UNSTABLE_FOCUS,        /* Re(λ) > 0, complex */
    EQUILIBRIUM_CENTER,                /* Re(λ) = 0, purely imaginary */
    EQUILIBRIUM_DEGENERATE,            /* Repeated eigenvalues, deficient */
    EQUILIBRIUM_ZERO_EIGENVALUE        /* At least one λ = 0 */
} EquilibriaClassification;

/** Numerical integration method type (L5: Computational Methods) */
typedef enum {
    METHOD_EULER_FORWARD,               /* Explicit Euler: O(h) */
    METHOD_EULER_BACKWARD,              /* Implicit Euler: O(h), A-stable */
    METHOD_HEUN,                        /* Improved Euler: O(h²) */
    METHOD_RK4_CLASSICAL,               /* Classical RK4: O(h⁴) */
    METHOD_RK45_DORMAND_PRINCE,         /* Adaptive RK5(4): O(h⁵) */
    METHOD_ADAMS_BASHFORTH_2,           /* AB-2: O(h²) */
    METHOD_ADAMS_BASHFORTH_3,           /* AB-3: O(h³) */
    METHOD_ADAMS_BASHFORTH_4,           /* AB-4: O(h⁴) */
    METHOD_ADAMS_MOULTON_2,             /* AM-2 (trapezoidal): O(h²) */
    METHOD_ADAMS_MOULTON_3,             /* AM-3: O(h³) */
    METHOD_BDF_2,                       /* BDF-2: O(h²), A-stable */
    METHOD_BDF_3,                       /* BDF-3: O(h³) */
    METHOD_USER_DEFINED                 /* External method */
} ODEMethod;

/** BVP solution method (L5: Computational Methods) */
typedef enum {
    BVP_SHOOTING_LINEAR,                /* Linear shooting method */
    BVP_SHOOTING_NONLINEAR,             /* Nonlinear shooting (Newton) */
    BVP_FINITE_DIFFERENCE_LINEAR,       /* Linear FD with tridiagonal solve */
    BVP_FINITE_DIFFERENCE_NONLINEAR     /* Nonlinear FD with Newton iteration */
} BVPMethod;

/** Bifurcation type (L8: Advanced Topics) */
typedef enum {
    BIFURCATION_NONE,                   /* No bifurcation detected */
    BIFURCATION_SADDLE_NODE,            /* Fold/Limit point */
    BIFURCATION_PITCHFORK_SUPERCRITICAL,
    BIFURCATION_PITCHFORK_SUBCRITICAL,
    BIFURCATION_HOPF_SUPERCRITICAL,     /* Stable limit cycle emerges */
    BIFURCATION_HOPF_SUBCRITICAL,       /* Unstable limit cycle */
    BIFURCATION_TRANSCRITICAL,
    BIFURCATION_HOMOCLINIC              /* Homoclinic orbit bifurcation */
} BifurcationType;

/* ──────────────────────── L1: Core Data Structures ───────────────────── */

/** A real-valued mathematical vector */
typedef struct {
    double *data;               /* Array of components */
    int     dim;                /* Dimension */
} ODEVector;

/** A real matrix (dense, row-major) */
typedef struct {
    double *data;               /* Flat array: data[i*cols + j] */
    int     rows;
    int     cols;
} Matrix;

/**
 * @brief ODE right-hand side function: dy/dt = f(t, y)
 *
 * @param t     Current time
 * @param y     Current state vector (length: dim)
 * @param f_out Output: derivative vector (length: dim)
 * @param dim   Dimension of state
 * @param ctx   User-defined context (parameters, etc.)
 * @return      0 on success, nonzero on error
 */
typedef int (*ODERHSFunc)(double t, const double *y, double *f_out,
                           int dim, void *ctx);

/**
 * @brief Jacobian of the vector field: J[i][j] = ∂f_i/∂y_j
 *
 * Reference: Wiggins "Introduction to Applied Nonlinear Dynamical Systems" (2003)
 *
 * @param t      Current time
 * @param y      Current state vector (length: dim)
 * @param J_out  Output: Jacobian matrix (dim × dim, row-major)
 * @param dim    Dimension
 * @param ctx    User context
 * @return       0 on success
 */
typedef int (*ODEJacobianFunc)(double t, const double *y, double *J_out,
                                int dim, void *ctx);

/**
 * @brief Event/root-finding function g(t,y) for detecting crossings.
 *
 * Used for event detection in integration (Poincaré sections, impacts).
 */
typedef int (*ODEEventFunc)(double t, const double *y, double *g_out,
                             int dim, void *ctx);

/* ──────────────────────── L2-L3: System Definitions ──────────────────── */

/**
 * @brief Complete specification of an initial value problem (IVP).
 *
 * ẏ = f(t, y),  y(t₀) = y₀,  t ∈ [t₀, t_end]
 *
 * Reference: Ascher & Petzold "Computer Methods for ODEs" (1998)
 */
typedef struct {
    ODERHSFunc      f;          /* RHS function */
    ODEJacobianFunc jac;        /* Jacobian (may be NULL for explicit methods) */
    ODEEventFunc    event;      /* Event function (may be NULL) */
    void           *ctx;        /* User context */
    double          t0;         /* Initial time */
    double          t_end;      /* Final time */
    double         *y0;         /* Initial condition (length: dim) */
    int             dim;        /* System dimension */
    ODEClass        ode_class;  /* Classification */
    ODEOrder        ode_order;  /* Order type */
} ODEIVP;

/**
 * @brief Simulation configuration for numerical integration.
 */
typedef struct {
    ODEMethod   method;         /* Integration method */
    double      h_init;         /* Initial step size */
    double      h_min;          /* Minimum allowed step (for adaptive methods) */
    double      h_max;          /* Maximum allowed step */
    double      rtol;           /* Relative tolerance */
    double      atol;           /* Absolute tolerance */
    int         max_steps;      /* Maximum number of integration steps */
    bool        dense_output;   /* Produce dense/interpolated output */
    double      dense_dt;       /* Output spacing for dense output (if enabled) */
} ODESimConfig;

/**
 * @brief Result of a single ODE integration step.
 */
typedef struct {
    double  t;                  /* Time at end of step */
    double *y;                  /* State at end of step (length: dim) */
    double  h_used;             /* Step size actually used */
    double  err_est;            /* Estimated local truncation error */
    int     fevals;             /* Number of RHS evaluations this step */
    int     jevals;             /* Number of Jacobian evaluations this step */
    bool    rejected;           /* Was this step rejected? */
} ODEStep;

/**
 * @brief Complete solution trajectory from numerical integration.
 */
typedef struct {
    double  *t_vals;            /* Time points (length: n_points) */
    double  *y_vals;            /* State values (n_points × dim, row-major) */
    int      n_points;          /* Number of saved points */
    int      dim;               /* Dimension */
    int      total_steps;       /* Total steps attempted */
    int      rejected_steps;    /* Steps rejected */
    int      total_fevals;      /* Total RHS evaluations */
    double   final_err;         /* Estimated error at final time */
    bool     success;           /* Integration completed successfully? */
    char     message[256];      /* Status message */
} ODESolution;

/* ──────────────────────── L2: Linear System ──────────────────────── */

/**
 * @brief Linear ODE system: ẏ = A y + b(t)
 *
 * A is a dim×dim matrix, b(t) is a time-dependent forcing vector.
 *
 * Reference: Chen, "Linear System Theory and Design" (1999)
 *            MIT 6.241 – Linear time-invariant systems
 */
typedef struct {
    double  *A;                 /* System matrix (dim × dim, row-major) */
    double  *b;                 /* Forcing vector (length: dim) or NULL if homogeneous */
    int      dim;               /* Dimension */
    bool     time_varying;      /* Is A(t) time-varying? */
} LinearODESystem;

/* ──────────────────────── L2: Phase Plane (2D) ─────────────────────── */

/**
 * @brief Phase plane point for 2D systems.
 *
 * Reference: Strogatz "Nonlinear Dynamics and Chaos" (2015)
 */
typedef struct {
    double  x;                  /* State variable 1 */
    double  y;                  /* State variable 2 */
} PhasePoint;

/**
 * @brief Phase portrait of a 2D system.
 */
typedef struct {
    PhasePoint *trajectories;   /* Array of trajectory points */
    int        *lengths;        /* Number of points per trajectory */
    int         n_trajectories; /* Number of trajectories */
    PhasePoint *equilibria;     /* Equilibrium points */
    int         n_equilibria;
    PhasePoint *nullcline_x;    /* x-nullcline points (where dx/dt = 0) */
    int         n_nullcline_x;
    PhasePoint *nullcline_y;    /* y-nullcline points (where dy/dt = 0) */
    int         n_nullcline_y;
} PhasePortrait;

/* ──────────────────────── L3: Vector Field ──────────────────────── */

/**
 * @brief Vector field for autonomous system ẏ = f(y).
 *
 * Reference: Arnold "Ordinary Differential Equations" (1992)
 *            MIT 6.243 – Nonlinear dynamics
 */
typedef struct {
    ODERHSFunc  f;              /* Vector field function */
    ODEJacobianFunc jac;        /* Jacobian of vector field */
    void       *ctx;            /* Parameters context */
    int         dim;            /* Dimension of phase space */
} VectorField;

/**
 * @brief Flow of a vector field: φ^t(y₀) is the solution at time t.
 */
typedef struct {
    VectorField *vf;
    double       t;             /* Flow time */
    double      *y0;            /* Starting point (length: dim) */
    double      *y_end;         /* End point after flow (length: dim) */
    int          dim;
} Flow;

/* ──────────────────────── L4: Wronskian ──────────────────────── */

/**
 * @brief Wronskian structure for analyzing linear independence of solutions.
 *
 * For solutions y₁(t), ..., yₙ(t) of an n-th order linear ODE:
 *   W(t) = det[ y₁  y₂  ...  yₙ;  y₁' y₂' ... yₙ';  ...;  y₁^{(n-1)} ... yₙ^{(n-1)} ]
 *
 * Reference: Coddington & Levinson "Theory of Ordinary Differential Equations" (1955)
 */
typedef struct {
    double  *fundamental_matrix; /* Φ(t): dim×dim fundamental matrix */
    double   wronskian;          /* det(Φ(t)) */
    double   t;                  /* Evaluation time */
    int      dim;                /* Dimension */
} WronskianInfo;

/* ──────────────────────── L6: Canonical ODE System Parameters ──────── */

/** Harmonic oscillator parameters: m·ẍ + c·ẋ + k·x = F(t) */
typedef struct {
    double  m;                  /* Mass (kg) */
    double  c;                  /* Damping coefficient (N·s/m) */
    double  k;                  /* Spring constant (N/m) */
    double  F0;                 /* Forcing amplitude */
    double  omega;              /* Forcing frequency (rad/s) */
} OscillatorParams;

/** Pendulum parameters: θ̈ + (g/L)·sin(θ) = 0 or with damping/forcing */
typedef struct {
    double  g;                  /* Gravitational acceleration (m/s²) */
    double  L;                  /* Pendulum length (m) */
    double  b;                  /* Damping coefficient (1/s) */
    double  A;                  /* Driving amplitude */
    double  omega_d;            /* Driving frequency (rad/s) */
} PendulumParams;

/** Van der Pol oscillator: ẍ - μ(1 - x²)ẋ + x = 0 */
typedef struct {
    double  mu;                 /* Nonlinearity parameter */
    double  omega;              /* Forcing frequency (0 if autonomous) */
    double  A;                  /* Forcing amplitude */
} VanDerPolParams;

/** Lotka-Volterra (predator-prey): ẋ = αx - βxy, ẏ = δxy - γy */
typedef struct {
    double  alpha;              /* Prey growth rate */
    double  beta;               /* Predation rate */
    double  gamma;              /* Predator death rate */
    double  delta;              /* Predator growth per prey eaten */
} LotkaVolterraParams;

/** Lorenz system: ẋ = σ(y - x), ẏ = x(ρ - z) - y, ż = xy - βz */
typedef struct {
    double  sigma;              /* Prandtl number */
    double  rho;                /* Rayleigh number */
    double  beta;               /* Geometric factor */
} LorenzParams;

/* ──────────────────────── L7: Application Parameters ────────────────── */

/** DC motor model: J·ω̇ + b·ω = K·i, L·i̇ + R·i = V - K·ω */
typedef struct {
    double  J;                  /* Rotor inertia (kg·m²) */
    double  b;                  /* Viscous friction (N·m·s) */
    double  K;                  /* Torque/EMF constant (N·m/A = V·s/rad) */
    double  L;                  /* Armature inductance (H) */
    double  R;                  /* Armature resistance (Ω) */
    double  V;                  /* Applied voltage (V) */
} DCMotorParams;

/** RLC circuit: L·ï + R·i̇ + (1/C)·i = v̇(t),  or state-space form */
typedef struct {
    double  L;                  /* Inductance (H) */
    double  R;                  /* Resistance (Ω) */
    double  C;                  /* Capacitance (F) */
    double  v0;                 /* Source voltage amplitude */
    double  omega;              /* Source frequency (rad/s) */
} RLCCircuitParams;

/* ──────────────────────── L8: Advanced Structures ───────────────────── */

/**
 * @brief Lyapunov function candidate V(x) for stability analysis.
 *
 * V: Rⁿ → R, must be positive definite with V(0) = 0.
 * Reference: Khalil "Nonlinear Systems" (2002)
 */
typedef struct {
    int     dim;                /* State dimension */
    int     (*V)(const double *x, int dim, void *ctx, double *value);
    int     (*grad_V)(const double *x, int dim, void *ctx, double *grad);
    void   *ctx;
    double  value;              /* Cached V(x) */
    double *gradient;           /* Cached ∇V(x) */
} LyapunovFunction;

/**
 * @brief Bifurcation analysis result for a parameter-dependent ODE.
 *
 * Analysis of ẋ = f(x; μ) as parameter μ varies.
 */
typedef struct {
    BifurcationType type;
    double   mu_crit;            /* Critical parameter value */
    double  *x_crit;             /* State at bifurcation point (length: dim) */
    double  *eigenvalues;        /* Eigenvalues at bifurcation */
    int      dim;
    double   lyapunov_coeff;     /* First Lyapunov coefficient (Hopf) */
} BifurcationResult;

#ifdef __cplusplus
}
#endif

#endif /* ODE_TYPES_H */
