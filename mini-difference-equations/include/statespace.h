/**
 * statespace.h — State-Space Representation of Discrete-Time Systems
 *
 * Domain: State-space models provide an internal description of
 * discrete-time dynamical systems, complementing the input-output
 * transfer function approach.
 *
 * Standard form:
 *   x[k+1] = A x[k] + B u[k]   (state equation)
 *   y[k]   = C x[k] + D u[k]   (output equation)
 *
 * Knowledge Level Coverage:
 *   L1: State-space quadruple (A,B,C,D), state vector, dimensions
 *   L2: Controllability, observability, state transition Φ(k) = A^k
 *   L3: Canonical forms (controllable, observable, diagonal, Jordan)
 *   L4: Cayley-Hamilton theorem (discrete), state transition properties
 *   L5: Ackermann pole placement, deadbeat control, similarity transforms
 *
 * Reference: Ogata, K. "Discrete-Time Control Systems" (1995)
 * Reference: Kailath, T. "Linear Systems" (1980)
 * Reference: Åström & Wittenmark "Computer-Controlled Systems" (1997)
 *
 * Course Mapping:
 *   MIT 6.241 — Dynamic Systems (state-space discrete)
 *   Stanford ENGR207B — Linear Control (controllability/observability)
 *   Berkeley ME 232 — Advanced Control (canonical forms)
 *   ETH 151-0555 — Linear Systems Theory
 */

#ifndef STATESPACE_H
#define STATESPACE_H

#include <stddef.h>
#include <complex.h>
#include "ztransform.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ── L1: State-Space Data Structures ─────────────────────────────────── */

/**
 * Discrete-time state-space model:
 *   x[k+1] = A x[k] + B u[k]
 *   y[k]   = C x[k] + D u[k]
 *
 * Dimensions:
 *   n = number of states
 *   m = number of inputs
 *   p = number of outputs
 */
typedef struct {
    int      n;          /* State dimension */
    int      m;          /* Input dimension */
    int      p;          /* Output dimension */
    double **A;          /* n×n state matrix */
    double **B;          /* n×m input matrix */
    double **C;          /* p×n output matrix */
    double **D;          /* p×m feedthrough matrix */
} StateSpace;

/**
 * State trajectory: sequence of state vectors x[0], x[1], ..., x[N].
 */
typedef struct {
    double **x;          /* x[k][i] = i-th state at time k */
    int      num_states; /* n */
    int      num_steps;  /* N+1 time steps */
} StateTrajectory;

/**
 * Controllability result: (A,B) is controllable if
 * rank(Co) = n where Co = [B, AB, A²B, ..., A^{n-1}B]
 */
typedef struct {
    double **Co;         /* n×(n·m) controllability matrix */
    int      rank;       /* Rank of Co */
    int      is_controllable; /* 1 if fully controllable */
} ControllabilityResult;

/**
 * Observability result: (A,C) is observable if
 * rank(Ob) = n where Ob = [C; CA; CA²; ...; CA^{n-1}]
 */
typedef struct {
    double **Ob;         /* (n·p)×n observability matrix */
    int      rank;       /* Rank of Ob */
    int      is_observable; /* 1 if fully observable */
} ObservabilityResult;

/**
 * Canonical form type.
 */
typedef enum {
    CANON_CONTROLLABLE,  /* Controllable canonical form */
    CANON_OBSERVABLE,    /* Observable canonical form */
    CANON_DIAGONAL,      /* Diagonal (modal) form */
    CANON_JORDAN         /* Jordan canonical form */
} CanonicalForm;

/* ── API: L3/L5 — State Transition ───────────────────────────────────── */

/**
 * Compute state transition matrix A^k for k ≥ 0.
 * Uses binary exponentiation for O(log k · n³) complexity.
 * Complexity: O(log k · n³)
 */
int state_transition(const double *const *A, int n, int k,
                     double **A_to_k);

/**
 * Simulate state-space system for N steps given initial state x0
 * and input sequence u[0..N-1].
 * x[k+1] = A x[k] + B u[k]
 * y[k]   = C x[k] + D u[k]
 * Complexity: O(N · n²)
 */
int ss_simulate(const StateSpace *sys,
                const double *x0,
                const double *const *u, int N,
                StateTrajectory *traj,
                double **y);

/**
 * Compute state vector at time k analytically:
 *   x[k] = A^k x[0] + Σ_{i=0}^{k-1} A^{k-1-i} B u[i]
 * Complexity: O(k · n²)
 */
int ss_state_at_k(const StateSpace *sys,
                  const double *x0,
                  const double *const *u, int k,
                  double *x_k);

/* ── API: L2/L5 — Controllability & Observability ────────────────────── */

/**
 * Build the controllability matrix:
 *   Co = [B, AB, A²B, ..., A^{n-1}B]  (n × n·m)
 * Complexity: O(n² · m · n) = O(n³ · m)
 */
int controllability_matrix(const StateSpace *sys,
                           ControllabilityResult *result);

/**
 * Check if (A,B) is controllable.
 * Uses rank check of controllability matrix via Gaussian elimination.
 * Complexity: O(n³ · m)
 */
int is_controllable(const StateSpace *sys);

/**
 * Build the observability matrix:
 *   Ob = [C; CA; CA²; ...; CA^{n-1}]  (n·p × n)
 * Complexity: O(n² · p · n) = O(n³ · p)
 */
int observability_matrix(const StateSpace *sys,
                         ObservabilityResult *result);

/**
 * Check if (A,C) is observable.
 * Complexity: O(n³ · p)
 */
int is_observable(const StateSpace *sys);

/* ── API: L3 — Canonical Forms ───────────────────────────────────────── */

/**
 * Convert transfer function to controllable canonical form:
 *   H(z) = (b_0 z^{n-1} + ... + b_{n-1}) / (z^n + a_1 z^{n-1} + ... + a_n)
 *
 *   A = [ -a_1  -a_2  ...  -a_n ]
 *       [   1     0   ...    0  ]
 *       [   0     1   ...    0  ]
 *       [  ...   ...  ...  ...  ]
 *       [   0     0   ...    0  ]
 *
 *   B = [1, 0, ..., 0]^T
 *   C = [b_0, b_1, ..., b_{n-1}]
 *   D = [0]
 *
 * Complexity: O(n²)
 */
int tf_to_ss_controllable(const ZTransferFunc *tf, StateSpace *sys);

/**
 * Convert transfer function to observable canonical form.
 * Complexity: O(n²)
 */
int tf_to_ss_observable(const ZTransferFunc *tf, StateSpace *sys);

/**
 * Convert state-space to transfer function:
 *   H(z) = C(zI - A)^{-1}B + D
 * Uses Faddeev-Leverrier algorithm.
 * Complexity: O(n⁴)
 */
int ss_to_tf(const StateSpace *sys, ZTransferFunc *tf);

/**
 * Apply similarity transformation: x̄ = T x
 *   Ā = T A T^{-1}, B̄ = T B, C̄ = C T^{-1}, D̄ = D
 * Complexity: O(n³)
 */
int ss_similarity_transform(const StateSpace *sys,
                             const double *const *T,
                             StateSpace *transformed);

/**
 * Diagonalize system matrix A via eigen-decomposition.
 * Only possible if A has n linearly independent eigenvectors.
 * Complexity: O(n³)
 */
int ss_diagonalize(StateSpace *sys);

/* ── API: L5 — Pole Placement ────────────────────────────────────────── */

/**
 * Ackermann formula for discrete-time pole placement.
 * Given controllable (A,B) and desired poles {p_1,...,p_n},
 * compute state feedback K such that A-BK has eigenvalues {p_i}.
 *
 * K = [0,...,0,1] · Co^{-1} · p(A)
 * where p(A) is the desired characteristic polynomial evaluated at A.
 *
 * Complexity: O(n³)
 */
int ackermann_pole_placement(const StateSpace *sys,
                              const double complex *desired_poles,
                              double **K);

/**
 * Deadbeat controller design: place all closed-loop poles at z = 0.
 * The system reaches equilibrium in at most n steps.
 * Complexity: O(n³)
 */
int deadbeat_controller(const StateSpace *sys, double **K);

/**
 * Deadbeat observer design: place all observer poles at z = 0.
 * Estimation error converges to zero in at most n steps.
 * Complexity: O(n³)
 */
int deadbeat_observer(const StateSpace *sys, double **L);

/* ── API: L3 — Discretization of State-Space ─────────────────────────── */

/**
 * Discretize continuous state-space (A_c, B_c, C_c, D_c) with ZOH:
 *   A = e^{A_c T}
 *   B = (∫_0^T e^{A_c τ} dτ) B_c
 *   C = C_c
 *   D = D_c
 *
 * Uses matrix exponential via Taylor series / scaling-and-squaring.
 * Complexity: O(n³) per matrix exponential step
 */
int ss_continuous_to_discrete(const double *const *Ac, int n,
                               const double *const *Bc, int m,
                               const double *const *Cc, int p,
                               const double *const *Dc,
                               double T,
                               StateSpace *disc_sys);

/* ── Utility ─────────────────────────────────────────────────────────── */

void statespace_free(StateSpace *sys);
StateSpace* statespace_alloc(int n, int m, int p);
void state_trajectory_free(StateTrajectory *traj);
void controllability_result_free(ControllabilityResult *cr);
void observability_result_free(ObservabilityResult *or);

/**
 * Print state-space matrices in readable format.
 */
void statespace_print(const StateSpace *sys);

/**
 * Matrix exponential: e^{At} via truncated Taylor series.
 * Complexity: O(terms · n³)
 */
int matrix_exponential(const double *const *A, int n, double t,
                       int num_terms, double **expAt);

#ifdef __cplusplus
}
#endif

#endif /* STATESPACE_H */
