/**
 * @file matrix.c
 * @brief Core matrix data structure implementation and fundamental operations
 * @module mini-linear-algebra-control
 *
 * Knowledge coverage:
 *   L1: Matrix create/free, element access, trace, determinant, rank
 *   L2: Matrix addition, multiplication, transpose, inverse, pseudoinverse
 *   L4: Cayley-Hamilton verification (matrix polynomial)
 *
 * Reference: Golub & Van Loan, Matrix Computations (4th Ed, 2013)
 */

#include "matrix.h"
#include "vector.h"
#include "decompositions.h"
#include <stdio.h>
#include <assert.h>

/* =========================================================================
 * L1: Construction and Destruction
 * ========================================================================= */

Matrix* matrix_create(int rows, int cols) {
    if (rows <= 0 || cols <= 0) return NULL;
    Matrix *m = (Matrix*)malloc(sizeof(Matrix));
    if (!m) return NULL;
    m->data = (double*)calloc((size_t)rows * cols, sizeof(double));
    if (!m->data) { free(m); return NULL; }
    m->rows = rows;
    m->cols = cols;
    m->owner = 1;
    return m;
}

Matrix* matrix_create_zero(int rows, int cols) {
    return matrix_create(rows, cols);
}

Matrix* matrix_create_identity(int n) {
    Matrix *m = matrix_create(n, n);
    if (!m) return NULL;
    for (int i = 0; i < n; i++) {
        matrix_set(m, i, i, 1.0);
    }
    return m;
}

Matrix* matrix_create_ones(int rows, int cols) {
    Matrix *m = matrix_create(rows, cols);
    if (!m) return NULL;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix_set(m, i, j, 1.0);
        }
    }
    return m;
}

Matrix* matrix_create_diag(const double *diag_vals, int n) {
    if (!diag_vals || n <= 0) return NULL;
    Matrix *m = matrix_create_zero(n, n);
    if (!m) return NULL;
    for (int i = 0; i < n; i++) {
        matrix_set(m, i, i, diag_vals[i]);
    }
    return m;
}

Matrix* matrix_create_copy(const Matrix *src) {
    if (!src) return NULL;
    Matrix *m = matrix_create(src->rows, src->cols);
    if (!m) return NULL;
    memcpy(m->data, src->data, (size_t)src->rows * src->cols * sizeof(double));
    return m;
}

Matrix* matrix_create_from_array(const double *vals, int rows, int cols) {
    if (!vals || rows <= 0 || cols <= 0) return NULL;
    Matrix *m = matrix_create(rows, cols);
    if (!m) return NULL;
    memcpy(m->data, vals, (size_t)rows * cols * sizeof(double));
    return m;
}

void matrix_free(Matrix *m) {
    if (!m) return;
    if (m->owner && m->data) {
        free(m->data);
    }
    free(m);
}

/* =========================================================================
 * L1: Element Access
 * ========================================================================= */

double matrix_get(const Matrix *m, int i, int j) {
    assert(m && m->data);
    assert(i >= 0 && i < m->rows);
    assert(j >= 0 && j < m->cols);
    return m->data[(size_t)i * m->cols + j];
}

void matrix_set(Matrix *m, int i, int j, double val) {
    assert(m && m->data);
    assert(i >= 0 && i < m->rows);
    assert(j >= 0 && j < m->cols);
    m->data[(size_t)i * m->cols + j] = val;
}

double* matrix_row_ptr(Matrix *m, int i) {
    assert(m && m->data);
    assert(i >= 0 && i < m->rows);
    return &m->data[(size_t)i * m->cols];
}

/* =========================================================================
 * L1: Dimension Queries and Matrix Properties
 * ========================================================================= */

int matrix_is_square(const Matrix *m) {
    if (!m) return 0;
    return m->rows == m->cols;
}

int matrix_is_symmetric(const Matrix *m, double tol) {
    if (!m || !matrix_is_square(m)) return 0;
    for (int i = 0; i < m->rows; i++) {
        for (int j = i + 1; j < m->cols; j++) {
            if (fabs(matrix_get(m, i, j) - matrix_get(m, j, i)) > tol) {
                return 0;
            }
        }
    }
    return 1;
}

int matrix_is_diagonal(const Matrix *m, double tol) {
    if (!m) return 0;
    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < m->cols; j++) {
            if (i != j && fabs(matrix_get(m, i, j)) > tol) return 0;
        }
    }
    return 1;
}

int matrix_is_upper_triangular(const Matrix *m, double tol) {
    if (!m) return 0;
    for (int i = 1; i < m->rows; i++) {
        for (int j = 0; j < i && j < m->cols; j++) {
            if (fabs(matrix_get(m, i, j)) > tol) return 0;
        }
    }
    return 1;
}

int matrix_is_lower_triangular(const Matrix *m, double tol) {
    if (!m) return 0;
    for (int i = 0; i < m->rows; i++) {
        for (int j = i + 1; j < m->cols; j++) {
            if (fabs(matrix_get(m, i, j)) > tol) return 0;
        }
    }
    return 1;
}

int matrix_is_orthogonal(const Matrix *m, double tol) {
    if (!m || !matrix_is_square(m)) return 0;
    int n = m->rows;
    Matrix *mt = matrix_create(n, n);
    Matrix *prod = matrix_create(n, n);
    if (!mt || !prod) { matrix_free(mt); matrix_free(prod); return 0; }
    matrix_transpose(m, mt);
    matrix_multiply(mt, m, prod);
    int is_orth = 1;
    for (int i = 0; i < n && is_orth; i++) {
        for (int j = 0; j < n && is_orth; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            if (fabs(matrix_get(prod, i, j) - expected) > tol) {
                is_orth = 0;
            }
        }
    }
    matrix_free(mt);
    matrix_free(prod);
    return is_orth;
}

/* =========================================================================
 * L2: Matrix Arithmetic
 * ========================================================================= */

void matrix_scale(Matrix *m, double alpha) {
    if (!m || !m->data) return;
    int n = m->rows * m->cols;
    for (int k = 0; k < n; k++) {
        m->data[k] *= alpha;
    }
}

void matrix_add_scalar(Matrix *m, double alpha) {
    if (!m || !m->data) return;
    int n = m->rows * m->cols;
    for (int k = 0; k < n; k++) {
        m->data[k] += alpha;
    }
}

void matrix_negate(Matrix *m) {
    matrix_scale(m, -1.0);
}

int matrix_add(const Matrix *A, const Matrix *B, Matrix *result) {
    if (!A || !B || !result) return -1;
    if (A->rows != B->rows || A->cols != B->cols) return -1;
    if (result->rows != A->rows || result->cols != A->cols) return -1;
    int n = A->rows * A->cols;
    for (int k = 0; k < n; k++) {
        result->data[k] = A->data[k] + B->data[k];
    }
    return 0;
}

int matrix_subtract(const Matrix *A, const Matrix *B, Matrix *result) {
    if (!A || !B || !result) return -1;
    if (A->rows != B->rows || A->cols != B->cols) return -1;
    if (result->rows != A->rows || result->cols != A->cols) return -1;
    int n = A->rows * A->cols;
    for (int k = 0; k < n; k++) {
        result->data[k] = A->data[k] - B->data[k];
    }
    return 0;
}

int matrix_multiply(const Matrix *A, const Matrix *B, Matrix *result) {
    if (!A || !B || !result) return -1;
    if (A->cols != B->rows) return -1;
    if (result->rows != A->rows || result->cols != B->cols) return -1;

    int m = A->rows, n = A->cols, p = B->cols;
    memset(result->data, 0, (size_t)m * p * sizeof(double));

    /* i-k-j loop order for cache efficiency */
    for (int i = 0; i < m; i++) {
        for (int k = 0; k < n; k++) {
            double aik = matrix_get(A, i, k);
            if (fabs(aik) < 1e-300) continue;
            double *res_row = &result->data[(size_t)i * p];
            double *b_row = &B->data[(size_t)k * p];
            for (int j = 0; j < p; j++) {
                res_row[j] += aik * b_row[j];
            }
        }
    }
    return 0;
}

int matrix_hadamard(const Matrix *A, const Matrix *B, Matrix *result) {
    if (!A || !B || !result) return -1;
    if (A->rows != B->rows || A->cols != B->cols) return -1;
    if (result->rows != A->rows || result->cols != A->cols) return -1;
    int n = A->rows * A->cols;
    for (int k = 0; k < n; k++) {
        result->data[k] = A->data[k] * B->data[k];
    }
    return 0;
}

/* =========================================================================
 * L2: Transpose
 * ========================================================================= */

int matrix_transpose(const Matrix *A, Matrix *result) {
    if (!A || !result) return -1;
    if (result->rows != A->cols || result->cols != A->rows) return -1;
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            matrix_set(result, j, i, matrix_get(A, i, j));
        }
    }
    return 0;
}

void matrix_transpose_inplace_square(Matrix *A) {
    if (!A || !matrix_is_square(A)) return;
    for (int i = 0; i < A->rows; i++) {
        for (int j = i + 1; j < A->cols; j++) {
            double tmp = matrix_get(A, i, j);
            matrix_set(A, i, j, matrix_get(A, j, i));
            matrix_set(A, j, i, tmp);
        }
    }
}

/* =========================================================================
 * L1: Trace - sum of diagonal entries
 * Linearity: tr(A+B) = tr(A) + tr(B); Cyclic: tr(AB) = tr(BA)
 * ========================================================================= */

double matrix_trace(const Matrix *A) {
    if (!A) return 0.0;
    int n = (A->rows < A->cols) ? A->rows : A->cols;
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += matrix_get(A, i, i);
    }
    return sum;
}

/* =========================================================================
 * L1: Determinant - recursive cofactor expansion for small n
 * For n <= 8: direct Laplace expansion (exact).
 * For n > 8: delegates to LU decomposition.
 * ========================================================================= */

static double det_recursive(const double *data, int n, int stride) {
    if (n == 1) return data[0];
    if (n == 2) return data[0] * data[stride + 1] - data[1] * data[stride];

    double det = 0.0;
    double *sub = (double*)malloc((size_t)(n - 1) * (n - 1) * sizeof(double));
    if (!sub) return 0.0;

    for (int j = 0; j < n; j++) {
        int idx = 0;
        for (int r = 1; r < n; r++) {
            for (int c = 0; c < n; c++) {
                if (c == j) continue;
                sub[idx++] = data[(size_t)r * stride + c];
            }
        }
        double cofactor = det_recursive(sub, n - 1, n - 1);
        det += (j % 2 == 0) ? data[j] * cofactor : -data[j] * cofactor;
    }
    free(sub);
    return det;
}

double matrix_determinant(const Matrix *A) {
    if (!A || !matrix_is_square(A)) return 0.0;
    int n = A->rows;

    if (n <= 8) {
        return det_recursive(A->data, n, n);
    }

    /* For larger matrices, use LU decomposition */
    Matrix *LU = matrix_create_copy(A);
    int *piv = (int*)malloc((size_t)n * sizeof(int));
    if (!LU || !piv) {
        matrix_free(LU);
        free(piv);
        return 0.0;
    }
    if (lu_decompose(LU, piv) == 0) {
        double det_sign = 1.0;
        for (int i = 0; i < n; i++) {
            if (piv[i] != i) det_sign = -det_sign;
        }
        double det = det_sign;
        for (int i = 0; i < n; i++) {
            det *= matrix_get(LU, i, i);
        }
        matrix_free(LU);
        free(piv);
        return det;
    }

    matrix_free(LU);
    free(piv);
    return 0.0;
}

/* =========================================================================
 * L1: Rank - via SVD (number of singular values above tolerance)
 * ========================================================================= */

int matrix_rank(const Matrix *A, double tol) {
    if (!A) return -1;
    if (tol <= 0) tol = 1e-10;

    int min_dim = (A->rows < A->cols) ? A->rows : A->cols;
    Matrix *Acopy = matrix_create_copy(A);
    Matrix *U = matrix_create(A->rows, A->rows);
    Vector *S = vector_create(min_dim);
    Matrix *Vt = matrix_create(A->cols, A->cols);

    if (!Acopy || !U || !S || !Vt) {
        matrix_free(Acopy); matrix_free(U);
        vector_free(S); matrix_free(Vt);
        return -1;
    }

    int result_rank = 0;
    if (svd_jacobi(Acopy, U, S, Vt, 1e-10, 100) >= 0) {
        result_rank = svd_effective_rank(S, tol);
    } else {
        /* Fallback: count non-negligible rows in R from QR-reduced form */
        for (int i = 0; i < min_dim; i++) {
            double norm = 0.0;
            for (int j = 0; j < Acopy->cols; j++) {
                double val = matrix_get(Acopy, i, j);
                norm += val * val;
            }
            if (sqrt(norm) > tol) result_rank++;
        }
    }

    matrix_free(Acopy); matrix_free(U);
    vector_free(S); matrix_free(Vt);
    return result_rank;
}

/* =========================================================================
 * L1: Matrix Power - A^k via exponentiation by squaring
 * ========================================================================= */

int matrix_power(const Matrix *A, int k, Matrix *result) {
    if (!A || !result || k < 0) return -1;
    if (!matrix_is_square(A)) return -1;
    int n = A->rows;

    Matrix *base = matrix_create_copy(A);
    Matrix *temp = matrix_create(n, n);
    if (!base || !temp) { matrix_free(base); matrix_free(temp); return -1; }

    /* Initialize result to identity */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix_set(result, i, j, (i == j) ? 1.0 : 0.0);
        }
    }

    int exp = k;
    while (exp > 0) {
        if (exp & 1) {
            matrix_multiply(result, base, temp);
            memcpy(result->data, temp->data, (size_t)n * n * sizeof(double));
        }
        matrix_multiply(base, base, temp);
        memcpy(base->data, temp->data, (size_t)n * n * sizeof(double));
        exp >>= 1;
    }

    matrix_free(base);
    matrix_free(temp);
    return 0;
}

/* =========================================================================
 * L4: Matrix Polynomial - p(A) = c0*I + c1*A + ... + cd*A^d
 * Cayley-Hamilton verification: for characteristic polynomial p,
 * p(A) = 0 (zero matrix). Every square matrix satisfies its own
 * characteristic equation.
 * ========================================================================= */

int matrix_polynomial(const Matrix *A, const double *coeffs,
                       int degree, Matrix *result) {
    if (!A || !coeffs || !result || degree < 0) return -1;
    if (!matrix_is_square(A)) return -1;
    int n = A->rows;
    if (result->rows != n || result->cols != n) return -1;

    /* result = c0 * I */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix_set(result, i, j, (i == j) ? coeffs[0] : 0.0);
        }
    }

    Matrix *Ak = matrix_create_identity(n);
    Matrix *temp = matrix_create(n, n);
    if (!Ak || !temp) { matrix_free(Ak); matrix_free(temp); return -1; }

    for (int k = 1; k <= degree; k++) {
        matrix_multiply(Ak, A, temp);
        memcpy(Ak->data, temp->data, (size_t)n * n * sizeof(double));

        if (fabs(coeffs[k]) > 1e-300) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double val = matrix_get(result, i, j)
                               + coeffs[k] * matrix_get(Ak, i, j);
                    matrix_set(result, i, j, val);
                }
            }
        }
    }

    matrix_free(Ak);
    matrix_free(temp);
    return 0;
}

/* =========================================================================
 * L2: Matrix Inverse via LU decomposition
 * ========================================================================= */

int matrix_inverse(const Matrix *A, Matrix *result) {
    if (!A || !result) return -1;
    if (!matrix_is_square(A)) return -1;
    int n = A->rows;
    if (result->rows != n || result->cols != n) return -1;
    return lu_inverse(A, result);
}

/* =========================================================================
 * L2: Moore-Penrose Pseudoinverse via SVD
 * Four Penrose conditions: AA^+A = A, A^+AA^+ = A^+, (AA^+)^T = AA^+, (A^+A)^T = A^+A
 * ========================================================================= */

int matrix_pseudoinverse(const Matrix *A, Matrix *result, double tol) {
    if (!A || !result) return -1;
    if (result->rows != A->cols || result->cols != A->rows) return -1;
    if (tol <= 0) tol = 1e-10;

    Matrix *Acopy = matrix_create_copy(A);
    Matrix *U = matrix_create(A->rows, A->rows);
    int min_dim = (A->rows < A->cols) ? A->rows : A->cols;
    Vector *S = vector_create(min_dim);
    Matrix *Vt = matrix_create(A->cols, A->cols);

    if (!Acopy || !U || !S || !Vt) {
        matrix_free(Acopy); matrix_free(U);
        vector_free(S); matrix_free(Vt);
        return -1;
    }

    if (svd_jacobi(Acopy, U, S, Vt, 1e-10, 100) < 0) {
        matrix_free(Acopy); matrix_free(U);
        vector_free(S); matrix_free(Vt);
        return -1;
    }

    /* Build V * S^+ * U^T: result = V * (S^+ * U^T) */
    Matrix *SuUt = matrix_create_zero(min_dim, A->rows);
    Matrix *V = matrix_create(A->cols, A->cols);
    if (!SuUt || !V) {
        matrix_free(Acopy); matrix_free(U); vector_free(S); matrix_free(Vt);
        matrix_free(SuUt); matrix_free(V);
        return -1;
    }

    for (int i = 0; i < min_dim; i++) {
        double sigma = vector_get(S, i);
        if (fabs(sigma) > tol) {
            double inv_sigma = 1.0 / sigma;
            for (int j = 0; j < A->rows; j++) {
                matrix_set(SuUt, i, j, inv_sigma * matrix_get(U, j, i));
            }
        }
    }

    matrix_transpose(Vt, V);
    matrix_multiply(V, SuUt, result);

    matrix_free(V); matrix_free(SuUt);
    matrix_free(Acopy); matrix_free(U);
    vector_free(S); matrix_free(Vt);
    return 0;
}

/* =========================================================================
 * L2: Kronecker Product A (x) B
 * ========================================================================= */

int matrix_kronecker(const Matrix *A, const Matrix *B, Matrix *result) {
    if (!A || !B || !result) return -1;
    int m = A->rows, n = A->cols, p = B->rows, q = B->cols;
    if (result->rows != m * p || result->cols != n * q) return -1;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double aij = matrix_get(A, i, j);
            for (int r = 0; r < p; r++) {
                for (int c = 0; c < q; c++) {
                    matrix_set(result, i * p + r, j * q + c,
                              aij * matrix_get(B, r, c));
                }
            }
        }
    }
    return 0;
}

/* =========================================================================
 * L2: Vectorization - vec(A) stacks columns into column vector
 * ========================================================================= */

int matrix_vec(const Matrix *A, Matrix *result) {
    if (!A || !result) return -1;
    if (result->rows != A->rows * A->cols || result->cols != 1) return -1;
    int idx = 0;
    for (int j = 0; j < A->cols; j++) {
        for (int i = 0; i < A->rows; i++) {
            matrix_set(result, idx++, 0, matrix_get(A, i, j));
        }
    }
    return 0;
}

/* =========================================================================
 * L2: Submatrix and Block Operations
 * ========================================================================= */

Matrix* matrix_submatrix(const Matrix *A, int row_start, int col_start,
                          int num_rows, int num_cols) {
    if (!A || row_start < 0 || col_start < 0 ||
        row_start + num_rows > A->rows || col_start + num_cols > A->cols) {
        return NULL;
    }
    Matrix *sub = matrix_create(num_rows, num_cols);
    if (!sub) return NULL;
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            matrix_set(sub, i, j, matrix_get(A, row_start + i, col_start + j));
        }
    }
    return sub;
}

int matrix_block_multiply(const Matrix *A11, const Matrix *A12,
                          const Matrix *A21, const Matrix *A22,
                          const Matrix *B11, const Matrix *B12,
                          const Matrix *B21, const Matrix *B22,
                          Matrix *C11, Matrix *C12,
                          Matrix *C21, Matrix *C22) {
    if (!A11 || !A12 || !A21 || !A22) return -1;
    if (!B11 || !B12 || !B21 || !B22) return -1;
    if (!C11 || !C12 || !C21 || !C22) return -1;

    Matrix *T1 = matrix_create(C11->rows, C11->cols);
    Matrix *T2 = matrix_create(C11->rows, C11->cols);
    if (!T1 || !T2) { matrix_free(T1); matrix_free(T2); return -1; }

    /* C11 = A11*B11 + A12*B21 */
    matrix_multiply(A11, B11, T1);
    matrix_multiply(A12, B21, T2);
    matrix_add(T1, T2, C11);

    /* C12 = A11*B12 + A12*B22 */
    matrix_multiply(A11, B12, T1);
    matrix_multiply(A12, B22, T2);
    matrix_add(T1, T2, C12);

    /* C21 = A21*B11 + A22*B21 */
    matrix_multiply(A21, B11, T1);
    matrix_multiply(A22, B21, T2);
    matrix_add(T1, T2, C21);

    /* C22 = A21*B12 + A22*B22 */
    matrix_multiply(A21, B12, T1);
    matrix_multiply(A22, B22, T2);
    matrix_add(T1, T2, C22);

    matrix_free(T1); matrix_free(T2);
    return 0;
}

/* =========================================================================
 * L2: Matrix Concatenation
 * ========================================================================= */

int matrix_horzcat(const Matrix *A, const Matrix *B, Matrix *result) {
    if (!A || !B || !result) return -1;
    if (A->rows != B->rows) return -1;
    if (result->rows != A->rows || result->cols != A->cols + B->cols) return -1;

    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            matrix_set(result, i, j, matrix_get(A, i, j));
        }
        for (int j = 0; j < B->cols; j++) {
            matrix_set(result, i, A->cols + j, matrix_get(B, i, j));
        }
    }
    return 0;
}

int matrix_vertcat(const Matrix *A, const Matrix *B, Matrix *result) {
    if (!A || !B || !result) return -1;
    if (A->cols != B->cols) return -1;
    if (result->rows != A->rows + B->rows || result->cols != A->cols) return -1;

    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            matrix_set(result, i, j, matrix_get(A, i, j));
        }
    }
    for (int i = 0; i < B->rows; i++) {
        for (int j = 0; j < B->cols; j++) {
            matrix_set(result, A->rows + i, j, matrix_get(B, i, j));
        }
    }
    return 0;
}

/* =========================================================================
 * L2: Matrix Equality and Printing
 * ========================================================================= */

int matrix_equal(const Matrix *A, const Matrix *B, double tol) {
    if (!A || !B) return 0;
    if (A->rows != B->rows || A->cols != B->cols) return 0;
    int n = A->rows * A->cols;
    for (int k = 0; k < n; k++) {
        if (fabs(A->data[k] - B->data[k]) > tol) return 0;
    }
    return 1;
}

void matrix_print(const Matrix *A, const char *name) {
    if (!A) { printf("%s = NULL\n", name ? name : "matrix"); return; }
    printf("%s (%dx%d):\n", name ? name : "matrix", A->rows, A->cols);
    for (int i = 0; i < A->rows; i++) {
        printf("  ");
        for (int j = 0; j < A->cols; j++) {
            printf("%12.6f ", matrix_get(A, i, j));
        }
        printf("\n");
    }
}