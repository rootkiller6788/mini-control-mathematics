/**
 * optimization_core.h — Core Definitions for Optimization
 *
 * Defines fundamental data structures and types used throughout
 * the mini-optimization-basics module. Covers L1 (Definitions):
 * convex sets, convex functions, gradient/Hessian structures,
 * optimality conditions structures.
 *
 * Reference: Boyd & Vandenberghe, "Convex Optimization" (2004)
 *            Nocedal & Wright, "Numerical Optimization" (2006)
 *
 * @module mini-optimization-basics
 * @level L1 (Definitions)
 */

#ifndef OPTIMIZATION_CORE_H
#define OPTIMIZATION_CORE_H

#include <stdlib.h>
#include <math.h>
#include <stdint.h>

/* ─── Constants ─────────────────────────────────────────────────────── */

#define OPT_EPSILON      1e-8
#define OPT_MAX_ITERS    10000
#define GRAD_TOL         1e-6
#define FEAS_TOL         1e-7
#define DUAL_GAP_TOL     1e-6

/* ─── Vector Type ─────────────────────────────────────────────────────── */

typedef struct {
    double *data;
    size_t  n;
} vector_t;

typedef struct {
    double *data;
    size_t  rows;
    size_t  cols;
} matrix_t;

/* Matrix element access: row-major (i*cols + j) */
#define MAT_ELT(M, i, j) ((M)->data[(i) * (M)->cols + (j)])

/* ─── Convex Set Definition (L1) ─────────────────────────────────────── */

typedef enum {
    CONVEX_SET_HALFSPACE,
    CONVEX_SET_POLYHEDRON,
    CONVEX_SET_L2_BALL,
    CONVEX_SET_L1_BALL,
    CONVEX_SET_LINF_BALL,
    CONVEX_SET_SIMPLEX,
    CONVEX_SET_PSD_CONE,
    CONVEX_SET_SOC,
    CONVEX_SET_BOX,
    CONVEX_SET_RAY,
    CONVEX_SET_USER_DEFINED
} convex_set_type_t;

typedef struct convex_set convex_set_t;

struct convex_set {
    convex_set_type_t type;
    double *param;
    size_t  param_len;
    int   (*membership)(const convex_set_t *C, const vector_t *x);
    void  (*project)(const convex_set_t *C, const vector_t *x, vector_t *proj);
    double (*support)(const convex_set_t *C, const vector_t *y);
};

/* ─── Convex Function Definition (L1) ─────────────────────────────────── */

typedef enum {
    CONVEX_FN_LINEAR,
    CONVEX_FN_QUADRATIC,
    CONVEX_FN_LOGISTIC,
    CONVEX_FN_EXPONENTIAL,
    CONVEX_FN_LOG_SUM_EXP,
    CONVEX_FN_NEG_ENTROPY,
    CONVEX_FN_NORM,
    CONVEX_FN_HUBER,
    CONVEX_FN_LOG_BARRIER,
    CONVEX_FN_INDICATOR,
    CONVEX_FN_USER
} convex_fn_type_t;

typedef struct convex_function convex_function_t;

struct convex_function {
    convex_fn_type_t type;
    size_t  n;
    double *param;
    size_t  param_len;
    double (*eval)(const convex_function_t *f, const vector_t *x);
    void   (*grad)(const convex_function_t *f, const vector_t *x, vector_t *g);
    void   (*hess)(const convex_function_t *f, const vector_t *x, matrix_t *H);
    double (*conjugate)(const convex_function_t *f, const vector_t *y);
    void   (*prox)(const convex_function_t *f, double lambda,
                   const vector_t *v, vector_t *prox_out);
};

/* ─── Gradient and Hessian (L3) ──────────────────────────────────────── */

typedef struct {
    vector_t grad;
    double   norm_l2;
    double   norm_linf;
} gradient_t;

typedef struct {
    matrix_t H;
    double   cond_number;
    double   min_eigval;
    double   max_eigval;
    int      is_psd;
} hessian_t;

/* ─── Optimality Conditions (L4) ──────────────────────────────────────── */

typedef struct {
    vector_t x_opt;
    vector_t grad_at_x;
    double   grad_norm;
    int      first_order_ok;
    double   min_eig_hess;
    int      second_order_ok;
} optimality_cert_t;

/* ─── Optimization Problem Types (L1/L6) ─────────────────────────────── */

typedef struct {
    size_t    n;
    convex_function_t f;
    vector_t  x0;
    double    opt_value;
    int       has_analytic;
} unconstr_problem_t;

typedef struct {
    size_t    n;
    convex_function_t f;
    matrix_t  A;
    vector_t  b;
    size_t    m;
    int       feasible;
} eq_constrained_problem_t;

typedef struct {
    size_t    n;
    convex_function_t f;
    convex_function_t *g;
    size_t    p;
    double    *lambdas;
} ineq_constrained_problem_t;

typedef struct {
    size_t    n;
    convex_function_t f;
    convex_function_t *g;
    size_t    p;
    double    *lambdas;
    double    *nus;
    matrix_t  A_eq;
    vector_t  b_eq;
    size_t    q;
    int       slater_holds;
} general_nlp_t;

/* ─── Linear Programming (L6) ────────────────────────────────────────── */

typedef struct {
    size_t    n;
    size_t    m;
    vector_t  c;
    matrix_t  A;
    vector_t  b;
    vector_t  x_opt;
    double    opt_value;
    int       bounded;
    int       feasible;
    vector_t  y_opt;
} linear_program_t;

/* ─── Quadratic Programming (L6) ──────────────────────────────────────── */

typedef struct {
    size_t    n;
    size_t    m;
    matrix_t  P;
    vector_t  q;
    double    r;
    matrix_t  A;
    vector_t  b;
    vector_t  x_opt;
    double    opt_value;
    int       psd;
} quadratic_program_t;

/* ─── Optimization Result ─────────────────────────────────────────────── */

typedef enum {
    OPT_SUCCESS = 0,
    OPT_MAXITER_REACHED,
    OPT_UNBOUNDED,
    OPT_INFEASIBLE,
    OPT_NUMERICAL_ERROR,
    OPT_NOT_CONVERGED,
    OPT_LINE_SEARCH_FAILED
} opt_status_t;

typedef struct {
    opt_status_t status;
    vector_t     x;
    double       f_value;
    double       grad_norm;
    size_t       iterations;
    double       primal_obj;
    double       dual_obj;
    double       dual_gap;
    int          kkt_satisfied;
    double       time_elapsed;
} opt_result_t;

/* ─── Line Search Structures ────────────────────────────────────────── */

typedef enum {
    LS_EXACT,
    LS_BACKTRACKING,
    LS_WOLFE,
    LS_STRONG_WOLFE,
    LS_GOLDSTEIN,
    LS_FIXED,
    LS_DIMINISHING
} linesearch_type_t;

typedef struct {
    double alpha_init;
    double rho;
    double c1;
    double c2;
    size_t max_iter;
} backtracking_params_t;

/* ─── Core API: Vector/Matrix Operations ────────────────────────────── */

vector_t *vector_alloc(size_t n);
vector_t *vector_from_array(const double *data, size_t n);
void vector_free(vector_t *v);
void vector_copy(const vector_t *src, vector_t *dst);
void vector_scale(vector_t *v, double alpha);
void vector_add(const vector_t *a, const vector_t *b, vector_t *out);
void vector_sub(const vector_t *a, const vector_t *b, vector_t *out);
double vector_dot(const vector_t *a, const vector_t *b);
double vector_norm_l2(const vector_t *v);
double vector_norm_linf(const vector_t *v);
double vector_norm_l1(const vector_t *v);
void vector_saxpy(double alpha, const vector_t *x, vector_t *y);

matrix_t *matrix_alloc(size_t rows, size_t cols);
void matrix_free(matrix_t *M);
void matrix_vector_mul(const matrix_t *A, const vector_t *x, vector_t *y);
matrix_t *matrix_transpose(const matrix_t *A);
double quadratic_form(const matrix_t *A, const vector_t *x);
int matrix_is_psd(const matrix_t *A);
int cholesky_decomp(matrix_t *A);
void cholesky_solve(const matrix_t *A, const vector_t *b, vector_t *x);
int symmetric_eig_range(const matrix_t *A, double *lambda_min, double *lambda_max);

#endif /* OPTIMIZATION_CORE_H */
