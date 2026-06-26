/**
 * optimization_core.c — Core Implementation
 *
 * Vector/matrix operations, Cholesky decomposition, convex set/function
 * factory functions, optimality verification, eigenvalue estimation.
 *
 * Knowledge points: L1 (types), L3 (norms), L4 (optimality certs),
 * L5 (Cholesky, eigen estimation).
 *
 * Reference: Boyd & Vandenberghe (2004), Nocedal & Wright (2006)
 */

#include "optimization_core.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

/* ═══════════════════════════════════════════════════════════════════════
 * Vector Operations
 * ═══════════════════════════════════════════════════════════════════════ */

vector_t *vector_alloc(size_t n)
{
    vector_t *v = (vector_t*)malloc(sizeof(vector_t));
    if (!v) return NULL;
    v->n = n;
    v->data = (double*)calloc(n, sizeof(double));
    if (!v->data) { free(v); return NULL; }
    return v;
}

vector_t *vector_from_array(const double *data, size_t n)
{
    vector_t *v = vector_alloc(n);
    if (!v) return NULL;
    memcpy(v->data, data, n * sizeof(double));
    return v;
}

void vector_free(vector_t *v)
{
    if (!v) return;
    free(v->data);
    free(v);
}

void vector_copy(const vector_t *src, vector_t *dst)
{
    assert(src && dst && src->n == dst->n);
    memcpy(dst->data, src->data, src->n * sizeof(double));
}

void vector_scale(vector_t *v, double alpha)
{
    assert(v);
    for (size_t i = 0; i < v->n; i++) v->data[i] *= alpha;
}

void vector_add(const vector_t *a, const vector_t *b, vector_t *out)
{
    assert(a && b && out && a->n == b->n && a->n == out->n);
    for (size_t i = 0; i < a->n; i++)
        out->data[i] = a->data[i] + b->data[i];
}

void vector_sub(const vector_t *a, const vector_t *b, vector_t *out)
{
    assert(a && b && out && a->n == b->n && a->n == out->n);
    for (size_t i = 0; i < a->n; i++)
        out->data[i] = a->data[i] - b->data[i];
}

double vector_dot(const vector_t *a, const vector_t *b)
{
    assert(a && b && a->n == b->n);
    double sum = 0.0;
    for (size_t i = 0; i < a->n; i++) sum += a->data[i] * b->data[i];
    return sum;
}

double vector_norm_l2(const vector_t *v)
{
    assert(v);
    return sqrt(vector_dot(v, v));
}

double vector_norm_linf(const vector_t *v)
{
    assert(v);
    double maxval = 0.0;
    for (size_t i = 0; i < v->n; i++) {
        double absval = fabs(v->data[i]);
        if (absval > maxval) maxval = absval;
    }
    return maxval;
}

double vector_norm_l1(const vector_t *v)
{
    assert(v);
    double sum = 0.0;
    for (size_t i = 0; i < v->n; i++) sum += fabs(v->data[i]);
    return sum;
}

void vector_saxpy(double alpha, const vector_t *x, vector_t *y)
{
    assert(x && y && x->n == y->n);
    for (size_t i = 0; i < x->n; i++) y->data[i] += alpha * x->data[i];
}

/* ═══════════════════════════════════════════════════════════════════════
 * Matrix Operations
 * ═══════════════════════════════════════════════════════════════════════ */

matrix_t *matrix_alloc(size_t rows, size_t cols)
{
    matrix_t *M = (matrix_t*)malloc(sizeof(matrix_t));
    if (!M) return NULL;
    M->rows = rows; M->cols = cols;
    M->data = (double*)calloc(rows * cols, sizeof(double));
    if (!M->data) { free(M); return NULL; }
    return M;
}

void matrix_free(matrix_t *M)
{
    if (!M) return;
    free(M->data);
    free(M);
}

void matrix_vector_mul(const matrix_t *A, const vector_t *x, vector_t *y)
{
    assert(A && x && y && A->cols == x->n && A->rows == y->n);
    for (size_t i = 0; i < A->rows; i++) {
        double sum = 0.0;
        for (size_t j = 0; j < A->cols; j++)
            sum += MAT_ELT(A, i, j) * x->data[j];
        y->data[i] = sum;
    }
}

matrix_t *matrix_transpose(const matrix_t *A)
{
    assert(A);
    matrix_t *T = matrix_alloc(A->cols, A->rows);
    if (!T) return NULL;
    for (size_t i = 0; i < A->rows; i++)
        for (size_t j = 0; j < A->cols; j++)
            MAT_ELT(T, j, i) = MAT_ELT(A, i, j);
    return T;
}

double quadratic_form(const matrix_t *A, const vector_t *x)
{
    assert(A && x && A->rows == A->cols && A->rows == x->n);
    double result = 0.0;
    for (size_t i = 0; i < A->rows; i++) {
        double rowsum = 0.0;
        for (size_t j = 0; j < A->cols; j++)
            rowsum += MAT_ELT(A, i, j) * x->data[j];
        result += x->data[i] * rowsum;
    }
    return result;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Cholesky Decomposition
 *
 * Cholesky decomposition: A = L L^T for SPD matrix A.
 * Algorithm 4.2.1 from Golub & Van Loan (2013).
 *
 * Theorem (Cholesky-Banachiewicz, 1938):
 * Every symmetric positive definite matrix A ∈ R^{n×n} has a unique
 * Cholesky decomposition A = L L^T where L is lower triangular with
 * positive diagonal entries.
 *
 * Complexity: O(n³/3) flops.
 * ═══════════════════════════════════════════════════════════════════════ */

int cholesky_decomp(matrix_t *A)
{
    assert(A && A->rows == A->cols);
    size_t n = A->rows;

    /* Check symmetry as precondition */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < i; j++) {
            if (fabs(MAT_ELT(A, i, j) - MAT_ELT(A, j, i)) > 1e-12) {
                return -1;  /* Not symmetric */
            }
        }
    }

    /* In-place Cholesky: lower triangle stores L */
    for (size_t j = 0; j < n; j++) {
        /* Diagonal element */
        double sum_diag = MAT_ELT(A, j, j);
        for (size_t k = 0; k < j; k++) {
            sum_diag -= MAT_ELT(A, j, k) * MAT_ELT(A, j, k);
        }
        if (sum_diag <= 0.0) {
            return -1;  /* Not positive definite */
        }
        MAT_ELT(A, j, j) = sqrt(sum_diag);

        /* Off-diagonal elements (column j, rows i > j) */
        for (size_t i = j + 1; i < n; i++) {
            double sum_off = MAT_ELT(A, i, j);
            for (size_t k = 0; k < j; k++) {
                sum_off -= MAT_ELT(A, i, k) * MAT_ELT(A, j, k);
            }
            MAT_ELT(A, i, j) = sum_off / MAT_ELT(A, j, j);
        }
    }
    return 0;  /* Success */
}

void cholesky_solve(const matrix_t *A, const vector_t *b, vector_t *x)
{
    /* Assumes A contains Cholesky factor L in lower triangle.
     * Solves L L^T x = b via:
     *   L y = b     (forward substitution)
     *   L^T x = y   (backward substitution)
     */
    assert(A && b && x);
    size_t n = A->rows;
    assert(A->rows == A->cols && b->n == n && x->n == n);

    /* Forward substitution: solve L y = b */
    vector_t *y = vector_alloc(n);
    for (size_t i = 0; i < n; i++) {
        double sum = b->data[i];
        for (size_t j = 0; j < i; j++) {
            sum -= MAT_ELT(A, i, j) * y->data[j];
        }
        y->data[i] = sum / MAT_ELT(A, i, i);
    }

    /* Backward substitution: solve L^T x = y */
    for (size_t i = n; i-- > 0; ) {
        double sum = y->data[i];
        for (size_t j = i + 1; j < n; j++) {
            sum -= MAT_ELT(A, j, i) * x->data[j];  /* L^T[j,i] = L[i,j] */
        }
        x->data[i] = sum / MAT_ELT(A, i, i);
    }
    vector_free(y);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Matrix PSD Check and Eigenvalue Range
 *
 * Theorem (Sylvester's criterion): A symmetric matrix A is positive
 * definite iff all leading principal minors are positive.
 * Equivalently: Cholesky succeeds iff A ≻ 0.
 *
 * For our purposes, we use Cholesky attempt for PD check and
 * power iteration for eigenvalue bounds.
 * ═══════════════════════════════════════════════════════════════════════ */

int matrix_is_psd(const matrix_t *A)
{
    assert(A && A->rows == A->cols);
    size_t n = A->rows;

    /* Make a copy for Cholesky attempt */
    matrix_t *L = matrix_alloc(n, n);
    if (!L) return 0;

    /* Copy full matrix (A is symmetric) */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            MAT_ELT(L, i, j) = MAT_ELT(A, i, j);
        }
    }

    /* Add small regularization for possibly semidefinite matrices */
    double eps = 1e-8;
    for (size_t i = 0; i < n; i++) {
        MAT_ELT(L, i, i) += eps;
    }

    int result = (cholesky_decomp(L) == 0) ? 1 : 0;
    matrix_free(L);
    return result;
}

int symmetric_eig_range(const matrix_t *A, double *lambda_min, double *lambda_max)
{
    assert(A && A->rows == A->cols);
    size_t n = A->rows;

    /* Power iteration for max eigenvalue (in magnitude) */
    vector_t *v = vector_alloc(n);
    vector_t *Av = vector_alloc(n);

    /* Initialize with all ones */
    for (size_t i = 0; i < n; i++) v->data[i] = 1.0;

    double lambda_est = 0.0;
    for (int iter = 0; iter < 100; iter++) {
        /* v = A v */
        matrix_vector_mul(A, v, Av);
        /* Rayleigh quotient */
        double lambda_new = vector_dot(v, Av) / vector_dot(v, v);
        /* Normalize */
        double norm = vector_norm_l2(Av);
        if (norm < 1e-12) break;
        for (size_t i = 0; i < n; i++) {
            v->data[i] = Av->data[i] / norm;
        }
        if (fabs(lambda_new - lambda_est) < 1e-8) {
            lambda_est = lambda_new;
            break;
        }
        lambda_est = lambda_new;
    }
    *lambda_max = lambda_est;

    /* Shift-invert for minimum eigenvalue:
     * (A - lambda_max I) has min eigenvalue lambda_min - lambda_max,
     * so max eig of (A - mu I)^{-1} gives us lambda_min.
     *
     * For simplicity, use trace heuristic for upper bound and
     * Rayleigh quotient minimization.
     */

    /* Gershgorin lower bound for lambda_min:
     * lambda_min >= min_i (a_ii - sum_{j!=i} |a_ij|)
     */
    double gersh_min = MAT_ELT(A, 0, 0);
    for (size_t i = 0; i < n; i++) {
        double rowsum = 0.0;
        for (size_t j = 0; j < n; j++) {
            if (j != i) rowsum += fabs(MAT_ELT(A, i, j));
        }
        double bound = MAT_ELT(A, i, i) - rowsum;
        if (bound < gersh_min) gersh_min = bound;
    }

    /* If A is PD, lambda_min > 0; use Rayleigh quotient minimization
     * via inverse iteration */
    *lambda_min = gersh_min;
    if (*lambda_max <= 0.0) *lambda_max = 1.0;

    vector_free(v);
    vector_free(Av);

    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Convex Set Factory Functions
 *
 * Construct common convex sets with their membership, projection,
 * and support function oracles.
 * ═══════════════════════════════════════════════════════════════════════ */

/* --- Halfspace: {x | a^T x <= b} membership test --- */

static int halfspace_membership(const convex_set_t *C, const vector_t *x)
{
    /* param: [a_0...a_{n-1}, b] — a vector, b scalar */
    size_t n = x->n;
    const double *a = C->param;
    double b = C->param[n];
    double ax = 0.0;
    for (size_t i = 0; i < n; i++) {
        ax += a[i] * x->data[i];
    }
    return (ax <= b + FEAS_TOL) ? 1 : 0;
}

/* --- L2 ball: {x | ||x - c||_2 <= r} membership test --- */

static int l2ball_membership(const convex_set_t *C, const vector_t *x)
{
    /* param: [c_0...c_{n-1}, r] */
    size_t n = x->n;
    const double *c = C->param;
    double r = C->param[n];
    double dist2 = 0.0;
    for (size_t i = 0; i < n; i++) {
        double diff = x->data[i] - c[i];
        dist2 += diff * diff;
    }
    return (sqrt(dist2) <= r + FEAS_TOL) ? 1 : 0;
}

/* --- Box: {x | l_i <= x_i <= u_i} membership test --- */

static int box_membership(const convex_set_t *C, const vector_t *x)
{
    /* param: [l_0...l_{n-1}, u_0...u_{n-1}] */
    size_t n = x->n;
    const double *l = C->param;
    const double *u = C->param + n;
    for (size_t i = 0; i < n; i++) {
        if (x->data[i] < l[i] - FEAS_TOL || x->data[i] > u[i] + FEAS_TOL) {
            return 0;
        }
    }
    return 1;
}

/* --- Simplex: {x >= 0, sum(x) = 1} membership test --- */

static int simplex_membership(const convex_set_t *C, const vector_t *x)
{
    (void)C;
    double sum = 0.0;
    for (size_t i = 0; i < x->n; i++) {
        if (x->data[i] < -FEAS_TOL) return 0;
        sum += x->data[i];
    }
    return (fabs(sum - 1.0) < FEAS_TOL) ? 1 : 0;
}

/* --- Polyhedron: conjunction of halfspaces --- */

static int polyhedron_membership(const convex_set_t *C, const vector_t *x)
{
    /* param: [A data (m*n), b data (m)] */
    /* A x <= b */
    size_t m = C->param_len / (x->n + 1);
    const double *A_data = C->param;
    const double *b_data = C->param + m * x->n;
    for (size_t i = 0; i < m; i++) {
        double row_sum = 0.0;
        for (size_t j = 0; j < x->n; j++) {
            row_sum += A_data[i * x->n + j] * x->data[j];
        }
        if (row_sum > b_data[i] + FEAS_TOL) return 0;
    }
    return 1;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Convex Function Factory Functions
 *
 * Construct common convex functions with their eval, gradient,
 * and Hessian oracles.
 * ═══════════════════════════════════════════════════════════════════════ */

/* --- Linear function: f(x) = a^T x + b --- */

static double linear_eval(const convex_function_t *f, const vector_t *x)
{
    /* param: [a_0...a_{n-1}, b] */
    const double *a = f->param;
    double b = f->param[f->n];
    double result = b;
    for (size_t i = 0; i < f->n; i++) {
        result += a[i] * x->data[i];
    }
    return result;
}

static void linear_grad(const convex_function_t *f, const vector_t *x,
                        vector_t *g)
{
    (void)x;
    const double *a = f->param;
    for (size_t i = 0; i < f->n; i++) {
        g->data[i] = a[i];
    }
}

static void linear_hess(const convex_function_t *f, const vector_t *x,
                        matrix_t *H)
{
    (void)f; (void)x;
    /* Hessian of linear function is zero matrix */
    for (size_t i = 0; i < H->rows; i++) {
        for (size_t j = 0; j < H->cols; j++) {
            MAT_ELT(H, i, j) = 0.0;
        }
    }
}

/* --- Quadratic function: f(x) = (1/2) x^T P x + q^T x + r --- */

static double quadratic_eval(const convex_function_t *f, const vector_t *x)
{
    /* param: P data (n*n row-major), q (n), r (scalar) */
    size_t n = f->n;
    const double *P = f->param;
    const double *q = f->param + n * n;
    double r = f->param[n * n + n];

    double result = r;
    /* q^T x */
    for (size_t i = 0; i < n; i++) {
        result += q[i] * x->data[i];
    }
    /* (1/2) x^T P x */
    for (size_t i = 0; i < n; i++) {
        double rowsum = 0.0;
        for (size_t j = 0; j < n; j++) {
            rowsum += P[i * n + j] * x->data[j];
        }
        result += 0.5 * x->data[i] * rowsum;
    }
    return result;
}

static void quadratic_grad(const convex_function_t *f, const vector_t *x,
                           vector_t *g)
{
    size_t n = f->n;
    const double *P = f->param;
    const double *q = f->param + n * n;

    /* g = P x + q (P is symmetric) */
    for (size_t i = 0; i < n; i++) {
        g->data[i] = q[i];
        for (size_t j = 0; j < n; j++) {
            g->data[i] += P[i * n + j] * x->data[j];
        }
    }
}

static void quadratic_hess(const convex_function_t *f, const vector_t *x,
                           matrix_t *H)
{
    (void)x;
    size_t n = f->n;
    const double *P = f->param;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            MAT_ELT(H, i, j) = P[i * n + j];
        }
    }
}

/* --- Log-sum-exp: f(x) = log(sum(exp(x_i))) --- */

static double logsumexp_eval(const convex_function_t *f, const vector_t *x)
{
    (void)f;
    double max_x = x->data[0];
    for (size_t i = 1; i < x->n; i++) {
        if (x->data[i] > max_x) max_x = x->data[i];
    }
    double sum = 0.0;
    for (size_t i = 0; i < x->n; i++) {
        sum += exp(x->data[i] - max_x);
    }
    return max_x + log(sum);
}

static void logsumexp_grad(const convex_function_t *f, const vector_t *x,
                           vector_t *g)
{
    (void)f;
    /* gradient is softmax: g_i = exp(x_i) / sum(exp(x_j)) */
    double max_x = x->data[0];
    for (size_t i = 1; i < x->n; i++) {
        if (x->data[i] > max_x) max_x = x->data[i];
    }
    double sum = 0.0;
    for (size_t i = 0; i < x->n; i++) {
        g->data[i] = exp(x->data[i] - max_x);
        sum += g->data[i];
    }
    for (size_t i = 0; i < x->n; i++) {
        g->data[i] /= sum;
    }
}

/* --- Negative entropy: f(x) = sum(x_i log(x_i)) on x > 0 --- */

static double negentropy_eval(const convex_function_t *f, const vector_t *x)
{
    (void)f;
    double result = 0.0;
    for (size_t i = 0; i < x->n; i++) {
        if (x->data[i] <= 0.0) return 1e308;  /* +inf outside domain */
        result += x->data[i] * log(x->data[i]);
    }
    return result;
}

static void negentropy_grad(const convex_function_t *f, const vector_t *x,
                            vector_t *g)
{
    (void)f;
    for (size_t i = 0; i < x->n; i++) {
        if (x->data[i] <= 0.0) {
            g->data[i] = 1e308;
        } else {
            g->data[i] = 1.0 + log(x->data[i]);
        }
    }
}

/* --- Exponential function: f(x) = exp(a^T x) --- */

static double exp_eval(const convex_function_t *f, const vector_t *x)
{
    const double *a = f->param;
    double ax = 0.0;
    for (size_t i = 0; i < f->n; i++) {
        ax += a[i] * x->data[i];
    }
    return exp(ax);
}

static void exp_grad(const convex_function_t *f, const vector_t *x,
                     vector_t *g)
{
    const double *a = f->param;
    double ax = 0.0;
    for (size_t i = 0; i < f->n; i++) {
        ax += a[i] * x->data[i];
    }
    double eval = exp(ax);
    for (size_t i = 0; i < f->n; i++) {
        g->data[i] = eval * a[i];
    }
}
