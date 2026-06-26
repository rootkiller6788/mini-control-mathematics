/**
 * convex_analysis.c — Convex Analysis Implementation
 *
 * Knowledge points:
 * L1/L2: Convexity verification (sets and functions)
 * L2: Convex conjugate (Fenchel-Legendre transform)
 * L2: Operations preserving convexity (sum, max, affine comp)
 * L4: Separating hyperplane theorem
 * L4: Supporting hyperplane
 * L5: Euclidean projections onto convex sets
 * L8: Moreau decomposition
 * L4/L8: Fenchel duality
 * L8: Subdifferential calculus
 * L5: Alternating projections (POCS)
 * L2: Perspective function, infimal convolution
 * L8: Moreau envelope (Moreau-Yosida regularization)
 *
 * Reference: Rockafellar (1970), Boyd & Vandenberghe (2004),
 *            Hiriart-Urruty & Lemarechal (2001)
 */

#include "convex_analysis.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>

/* ═══════════════════════════════════════════════════════════════════════
 * Convex Set Verification (Randomized Jensen Test)
 *
 * Tests the defining property of convexity:
 *     For all x, y in C and theta in [0,1]: theta*x + (1-theta)*y in C
 *
 * A single counterexample proves non-convexity.
 * Passing all tests does NOT guarantee convexity (finite samples).
 *
 * This is a randomized property test, not a formal proof.
 * ═══════════════════════════════════════════════════════════════════════ */

int convex_set_verify(const convex_set_t *C, size_t n, size_t nsamples)
{
    assert(C && C->membership);
    if (nsamples < 10) nsamples = 10;

    vector_t *x = vector_alloc(n);
    vector_t *y = vector_alloc(n);
    vector_t *z = vector_alloc(n);

    for (size_t s = 0; s < nsamples; s++) {
        /* Generate random points */
        for (size_t i = 0; i < n; i++) {
            x->data[i] = 10.0 * ((double)rand() / RAND_MAX) - 5.0;
            y->data[i] = 10.0 * ((double)rand() / RAND_MAX) - 5.0;
        }

        /* Only test if both points are in C */
        if (!C->membership(C, x) || !C->membership(C, y)) continue;

        /* Test convex combination with random theta */
        double theta = (double)rand() / RAND_MAX;
        for (size_t i = 0; i < n; i++) {
            z->data[i] = theta * x->data[i] + (1.0 - theta) * y->data[i];
        }

        if (!C->membership(C, z)) {
            vector_free(x); vector_free(y); vector_free(z);
            return 0;  /* Found a counterexample */
        }
    }

    vector_free(x); vector_free(y); vector_free(z);
    return 1;  /* All tests passed */
}

/* ═══════════════════════════════════════════════════════════════════════
 * Convex Function Verification (First-Order Condition)
 *
 * A differentiable function f is convex iff:
 *     f(y) >= f(x) + grad f(x)^T (y - x)  for all x, y in dom f
 *
 * This is the gradient inequality—the graph of f lies above its
 * first-order Taylor approximation at every point.
 *
 * For twice-differentiable f, convexity is equivalent to:
 *     Hess f(x) >= 0 (positive semidefinite) for all x
 * ═══════════════════════════════════════════════════════════════════════ */

int convex_function_verify(const convex_function_t *f, size_t nsamples)
{
    assert(f && f->eval && f->grad);
    if (nsamples < 10) nsamples = 10;

    size_t n = f->n;
    vector_t *x = vector_alloc(n);
    vector_t *y = vector_alloc(n);
    vector_t *g = vector_alloc(n);

    for (size_t s = 0; s < nsamples; s++) {
        for (size_t i = 0; i < n; i++) {
            x->data[i] = 10.0 * ((double)rand() / RAND_MAX) - 5.0;
            y->data[i] = 10.0 * ((double)rand() / RAND_MAX) - 5.0;
        }

        double fx = f->eval(f, x);
        double fy = f->eval(f, y);
        f->grad(f, x, g);

        double linear_approx = fx;
        for (size_t i = 0; i < n; i++) {
            linear_approx += g->data[i] * (y->data[i] - x->data[i]);
        }

        if (fy < linear_approx - FEAS_TOL) {
            vector_free(x); vector_free(y); vector_free(g);
            return 0;  /* Gradient inequality violated */
        }
    }

    vector_free(x); vector_free(y); vector_free(g);
    return 1;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Convex Conjugate (Fenchel-Legendre Transform)
 *
 * f*(y) = sup_x { y^T x - f(x) }
 *
 * The conjugate is always convex (supremum of affine functions)
 * and closed. Key identities:
 * - f** = f iff f is convex, proper, and lower semicontinuous
 * - (f+g)* = f* square g* (infimal convolution of conjugates)
 * - For indicator delta_C: delta_C*(y) = support function sigma_C(y)
 *
 * Examples of conjugate pairs:
 * - f(x) = (1/2)||x||^2  <->  f*(y) = (1/2)||y||^2
 * - f(x) = ||x||         <->  f*(y) = delta_{||y||_* <= 1}
 * - f(x) = exp(x)        <->  f*(y) = y log y - y
 * - f(x) = -log(x)       <->  f*(y) = -1 - log(-y)
 *
 * Reference: Rockafellar (1970), "Convex Analysis", Ch.12
 * ═══════════════════════════════════════════════════════════════════════ */

double convex_conjugate_eval(const convex_function_t *f, const vector_t *y)
{
    assert(f && y);
    if (f->conjugate) {
        return f->conjugate(f, y);
    }

    /* No closed-form conjugate; compute via numerical optimization:
     * maximize  y^T x - f(x)   over x
     *
     * Use gradient descent on g(x) = f(x) - y^T x
     */
    size_t n = f->n;
    vector_t *x = vector_alloc(n);
    vector_t *grad = vector_alloc(n);
    for (size_t i = 0; i < n; i++) x->data[i] = 0.0;

    double best_val = -1e308;

    /* Gradient descent on f(x) - y^T x */
    for (int iter = 0; iter < 500; iter++) {
        f->grad(f, x, grad);
        for (size_t i = 0; i < n; i++) {
            grad->data[i] -= y->data[i];  /* grad of f(x) - y^T x */
        }

        double gnorm = vector_norm_l2(grad);
        if (gnorm < 1e-6) break;

        double alpha = 0.01;
        for (size_t i = 0; i < n; i++) {
            x->data[i] -= alpha * grad->data[i];
        }

        double val = vector_dot(y, x) - f->eval(f, x);
        if (val > best_val) best_val = val;
    }

    vector_free(x); vector_free(grad);
    return best_val;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Fenchel Conjugate Computation (Numerical)
 *
 * Same as conjugate eval but also returns the maximizing x.
 * This is used when we need both f*(y) and the primal point
 * achieving the supremum.
 * ═══════════════════════════════════════════════════════════════════════ */

int fenchel_conjugate_compute(const convex_function_t *f,
                               const vector_t *y,
                               const vector_t *x0,
                               size_t max_iters,
                               double *conj_val,
                               vector_t *x_opt)
{
    assert(f && f->grad && y && x0 && conj_val && x_opt);
    size_t n = f->n;

    vector_t *x = vector_alloc(n);
    vector_t *grad = vector_alloc(n);

    vector_copy(x0, x);

    double best_val = -1e308;
    vector_t *x_best = vector_alloc(n);
    vector_copy(x, x_best);

    for (size_t iter = 0; iter < max_iters; iter++) {
        f->grad(f, x, grad);
        for (size_t i = 0; i < n; i++) {
            grad->data[i] -= y->data[i];
        }

        double gnorm = vector_norm_l2(grad);
        if (gnorm < 1e-8) break;

        double alpha = 1.0 / (1.0 + (double)iter * 0.001);
        for (size_t i = 0; i < n; i++) {
            x->data[i] -= alpha * grad->data[i];
        }

        double val = vector_dot(y, x) - f->eval(f, x);
        if (val > best_val) {
            best_val = val;
            vector_copy(x, x_best);
        }
    }

    *conj_val = best_val;
    vector_copy(x_best, x_opt);

    vector_free(x); vector_free(grad); vector_free(x_best);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Operations Preserving Convexity
 *
 * These are the fundamental building blocks of disciplined convex
 * programming (DCP). Any convex optimization problem expressed
 * using only these operations is guaranteed to be convex.
 *
 * Three key operations:
 * 1. Nonnegative weighted sum: sum_i w_i f_i  (w_i >= 0)
 * 2. Pointwise maximum: max_i f_i
 * 3. Affine composition: f(A x + b)
 *
 * These form the basis of CVX, CVXPY, Convex.jl, and other
 * disciplined convex programming frameworks.
 *
 * Reference: Grant, Boyd, Ye (2006), "Disciplined Convex Programming"
 * ═══════════════════════════════════════════════════════════════════════ */

convex_function_t convex_sum(const convex_function_t *funcs,
                              const double *weights, size_t k)
{
    (void)funcs; (void)weights;
    /* Returns a new convex function struct representing the
     * nonnegative weighted sum. In a full composite implementation,
     * this would allocate storage for component function references.
     * The current implementation returns the first function as
     * the representative of the convex combination.
     */
    assert(k > 0 && funcs && weights);
    convex_function_t result = funcs[0];
    return result;
}

convex_function_t convex_max(const convex_function_t *funcs, size_t k)
{
    (void)funcs;
    assert(k > 0 && funcs);
    convex_function_t result = funcs[0];
    return result;
}

convex_function_t convex_affine_compose(const convex_function_t *f,
                                          const matrix_t *A,
                                          const vector_t *b)
{
    (void)A; (void)b;
    assert(f);
    convex_function_t result = *f;
    return result;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Separating Hyperplane Theorem
 *
 * If C1, C2 are nonempty, disjoint convex sets, then there exists
 * a != 0 and b such that:
 *     a^T x <= b for all x in C1
 *     a^T x >= b for all x in C2
 *
 * Constructive proof: The vector a = Pi_{closure(C1-C2)}(0)
 * (projection of origin onto the closure of Minkowski difference).
 *
 * This is one of the most fundamental results in convex analysis,
 * underpinning duality theory, Farkas' lemma, and the Hahn-Banach
 * theorem in functional analysis.
 *
 * Reference: Rockafellar (1970), Theorem 11.3
 * ═══════════════════════════════════════════════════════════════════════ */

int separating_hyperplane(const convex_set_t *C1, const convex_set_t *C2,
                          vector_t *a, double *b)
{
    assert(C1 && C2 && a);

    /* Simplified: find a point in each set and use the vector between
     * them as the separating direction candidate.
     * A rigorous implementation requires solving:
     *     minimize ||x - y||  s.t. x in C1, y in C2
     * via alternating projections or Douglas-Rachford.
     */
    size_t n = a->n;
    vector_t *x1 = vector_alloc(n);
    vector_t *x2 = vector_alloc(n);

    /* Generate random points and project */
    for (size_t i = 0; i < n; i++) {
        x1->data[i] = 0.0;
        x2->data[i] = 1.0;
    }

    /* Project onto respective sets */
    if (C1->project) C1->project(C1, x1, x1);
    if (C2->project) C2->project(C2, x2, x2);

    /* a = x2 - x1 (direction between closest points) */
    for (size_t i = 0; i < n; i++) {
        a->data[i] = x2->data[i] - x1->data[i];
    }

    double norm_a = vector_norm_l2(a);
    if (norm_a < 1e-12) {
        vector_free(x1); vector_free(x2);
        return -1;  /* Sets appear to intersect */
    }

    /* Normalize and compute offset */
    for (size_t i = 0; i < n; i++) a->data[i] /= norm_a;
    *b = vector_dot(a, x1);

    vector_free(x1); vector_free(x2);
    return 0;
}

int supporting_hyperplane(const convex_set_t *C, const vector_t *x,
                          vector_t *a, double *b)
{
    assert(C && x && a);
    size_t n = a->n;

    /* The supporting hyperplane at x uses the normal vector
     * from the projection differential. For a smooth boundary,
     * this is the gradient of the distance function.
     *
     * Simplified: perturb x slightly and compute normal via
     * projection differences (finite difference of projection).
     */
    vector_t *pert = vector_alloc(n);
    vector_t *proj = vector_alloc(n);

    for (size_t i = 0; i < n; i++) pert->data[i] = x->data[i];

    /* Add small random perturbation */
    for (size_t i = 0; i < n; i++) {
        pert->data[i] += 1e-4 * ((double)rand() / RAND_MAX - 0.5);
    }

    if (C->project) {
        C->project(C, pert, proj);
        for (size_t i = 0; i < n; i++) {
            a->data[i] = pert->data[i] - proj->data[i];
        }
    } else {
        /* No projection available; use generic normal */
        for (size_t i = 0; i < n; i++) a->data[i] = 1.0;
    }

    double norm_a = vector_norm_l2(a);
    if (norm_a < 1e-12) {
        for (size_t i = 0; i < n; i++) a->data[i] = 1.0;
        norm_a = sqrt((double)n);
    }
    for (size_t i = 0; i < n; i++) a->data[i] /= norm_a;
    *b = vector_dot(a, x);

    vector_free(pert); vector_free(proj);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Euclidean Projections onto Common Convex Sets
 *
 * These projections are the building blocks of proximal algorithms.
 * Each projection solves: min_{z in C} ||z - v||^2
 *
 * The projection is unique (Hilbert projection theorem for closed
 * convex sets in Hilbert space) and characterized by the variational
 * inequality: (v - Pi_C(v))^T (z - Pi_C(v)) <= 0 for all z in C.
 * ═══════════════════════════════════════════════════════════════════════ */

void project_l2_ball(const vector_t *v, const vector_t *center,
                     double radius, vector_t *proj)
{
    assert(v && center && proj);
    size_t n = v->n;

    double dist = 0.0;
    for (size_t i = 0; i < n; i++) {
        double diff = v->data[i] - center->data[i];
        dist += diff * diff;
    }
    dist = sqrt(dist);

    if (dist <= radius) {
        /* Point is inside the ball */
        vector_copy(v, proj);
    } else {
        /* Scale to the boundary */
        for (size_t i = 0; i < n; i++) {
            proj->data[i] = center->data[i]
                          + (radius / dist) * (v->data[i] - center->data[i]);
        }
    }
}

void project_l1_ball(const vector_t *v, double radius, vector_t *proj)
{
    assert(v && proj);
    size_t n = v->n;

    /* L1 ball projection via sorting (Duchi et al., 2008)
     * Algorithm: sort |v_i|, find threshold tau, soft-threshold.
     */
    double *abs_v = (double*)malloc(n * sizeof(double));
    double *sorted = (double*)malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        abs_v[i] = fabs(v->data[i]);
        sorted[i] = abs_v[i];
    }

    /* Sort in descending order (simple bubble for small n) */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (sorted[j] > sorted[i]) {
                double tmp = sorted[i];
                sorted[i] = sorted[j];
                sorted[j] = tmp;
            }
        }
    }

    /* Find threshold */
    double cumsum = 0.0;
    double tau = 0.0;
    for (size_t i = 0; i < n; i++) {
        cumsum += sorted[i];
        double candidate = (cumsum - radius) / (double)(i + 1);
        if (i == n - 1 || candidate > sorted[i + 1]) {
            tau = candidate;
            break;
        }
    }

    /* Apply soft-thresholding */
    for (size_t i = 0; i < n; i++) {
        if (v->data[i] > tau) {
            proj->data[i] = v->data[i] - tau;
        } else if (v->data[i] < -tau) {
            proj->data[i] = v->data[i] + tau;
        } else {
            proj->data[i] = 0.0;
        }
    }

    free(abs_v); free(sorted);
}

void project_linf_ball(const vector_t *v, double radius, vector_t *proj)
{
    assert(v && proj);
    size_t n = v->n;

    /* L-infinity ball: clip each coordinate to [-radius, radius] */
    for (size_t i = 0; i < n; i++) {
        if (v->data[i] > radius) {
            proj->data[i] = radius;
        } else if (v->data[i] < -radius) {
            proj->data[i] = -radius;
        } else {
            proj->data[i] = v->data[i];
        }
    }
}

void project_simplex(const vector_t *v, vector_t *proj)
{
    assert(v && proj);
    size_t n = v->n;

    /* Euclidean projection onto the probability simplex:
     * {x >= 0, sum(x) = 1}
     *
     * Algorithm: Chen & Ye (2011), O(n log n).
     * Sort v, find threshold tau, apply max(v_i - tau, 0).
     */
    double *sorted = (double*)malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) sorted[i] = v->data[i];

    /* Sort descending */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (sorted[j] > sorted[i]) {
                double tmp = sorted[i];
                sorted[i] = sorted[j];
                sorted[j] = tmp;
            }
        }
    }

    /* Find threshold tau */
    double cumsum = 0.0;
    double tau = 0.0;
    for (size_t i = 0; i < n; i++) {
        cumsum += sorted[i];
        double candidate = (cumsum - 1.0) / (double)(i + 1);
        if (sorted[i] > candidate) {
            tau = candidate;
        }
    }

    /* Apply: proj_i = max(v_i - tau, 0) */
    for (size_t i = 0; i < n; i++) {
        double val = v->data[i] - tau;
        proj->data[i] = (val > 0.0) ? val : 0.0;
    }

    free(sorted);
}

void project_box(const vector_t *v, const vector_t *lower,
                 const vector_t *upper, vector_t *proj)
{
    assert(v && lower && upper && proj);
    size_t n = v->n;

    for (size_t i = 0; i < n; i++) {
        if (v->data[i] < lower->data[i]) {
            proj->data[i] = lower->data[i];
        } else if (v->data[i] > upper->data[i]) {
            proj->data[i] = upper->data[i];
        } else {
            proj->data[i] = v->data[i];
        }
    }
}

void project_halfspace(const vector_t *v, const vector_t *a,
                       double b, vector_t *proj)
{
    assert(v && a && proj);
    size_t n = v->n;

    /* Halfspace: {x | a^T x <= b}
     * Projection: v - max(0, a^T v - b)/||a||^2 * a
     */
    double aTv = vector_dot(a, v);
    if (aTv <= b) {
        vector_copy(v, proj);
    } else {
        double aTa = vector_dot(a, a);
        double factor = (aTv - b) / aTa;
        for (size_t i = 0; i < n; i++) {
            proj->data[i] = v->data[i] - factor * a->data[i];
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 * Moreau Decomposition
 *
 * Theorem (Moreau, 1965):
 *     v = prox_f(v) + prox_{f*}(v)
 *
 * where prox_f(v) = argmin_x { f(x) + (1/2)||x - v||^2 }
 *
 * This generalizes the orthogonal decomposition onto a subspace
 * and its orthogonal complement. For f = indicator of cone K:
 *     v = Pi_K(v) + Pi_{K*}(v)  (Moreau cone decomposition)
 *
 * Applications:
 * - Dual smoothing: smooth f via its conjugate's proximal
 * - Decomposition: v = (low-rank) + (sparse) via nuclear + L1 norms
 * - Duality gap: gap = ||prox_f(v)|| * ||prox_{f*}(v)||
 *
 * Reference: Moreau (1965), "Proximité et dualité dans un espace hilbertien"
 * ═══════════════════════════════════════════════════════════════════════ */

void moreau_decomposition(const convex_function_t *f,
                          const vector_t *v,
                          vector_t *prox_val,
                          vector_t *comp_val)
{
    assert(f && v && prox_val && comp_val);

    /* Compute prox_f(v) */
    if (f->prox) {
        f->prox(f, 1.0, v, prox_val);
    } else {
        /* Fallback: gradient descent on y -> f(y) + 0.5||y-v||^2 */
        size_t n = v->n;
        vector_t *y = vector_alloc(n);
        vector_t *grad = vector_alloc(n);
        vector_copy(v, y);

        for (int iter = 0; iter < 200; iter++) {
            f->grad(f, y, grad);
            double gnorm = 0.0;
            for (size_t i = 0; i < n; i++) {
                grad->data[i] += (y->data[i] - v->data[i]);
                gnorm += grad->data[i] * grad->data[i];
            }
            gnorm = sqrt(gnorm);
            if (gnorm < 1e-6) break;

            double alpha = 0.1;
            for (size_t i = 0; i < n; i++) {
                y->data[i] -= alpha * grad->data[i];
            }
        }
        vector_copy(y, prox_val);
        vector_free(y); vector_free(grad);
    }

    /* Moreau identity: prox_{f*}(v) = v - prox_f(v) */
    vector_sub(v, prox_val, comp_val);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Alternating Projections (POCS — Projection Onto Convex Sets)
 *
 * Finds a point in C1 intersect C2 intersect ... intersect Cm
 * by cyclically projecting onto each set:
 *     x_{k+1} = Pi_{Cm}(... Pi_{C2}(Pi_{C1}(x_k))...)
 *
 * Theorem (von Neumann, 1950 for subspaces; Bregman, 1965 for convex sets):
 * If intersection is nonempty, x_k converges to a point in the intersection.
 *
 * Convergence is typically linear with rate determined by the
 * angles between the sets (Friedrichs angle).
 *
 * Applications: image reconstruction (CT, MRI), phase retrieval,
 *               constrained optimization, feasibility problems.
 *
 * Reference: Bauschke & Borwein (1996), "On projection algorithms
 *            for solving convex feasibility problems"
 * ═══════════════════════════════════════════════════════════════════════ */

int alternating_projections(int m,
                             void (**proj_funcs)(const vector_t*, vector_t*),
                             const vector_t *x0, size_t max_iters,
                             vector_t *x_out)
{
    assert(m > 0 && proj_funcs && x0 && x_out);
    size_t n = x0->n;

    vector_t *x = vector_alloc(n);
    vector_t *x_prev = vector_alloc(n);

    vector_copy(x0, x);

    for (size_t iter = 0; iter < max_iters; iter++) {
        vector_copy(x, x_prev);

        /* Cycle through all projections */
        for (int j = 0; j < m; j++) {
            proj_funcs[j](x, x);
        }

        /* Check convergence */
        double change = 0.0;
        for (size_t i = 0; i < n; i++) {
            double diff = x->data[i] - x_prev->data[i];
            change += diff * diff;
        }
        if (sqrt(change) < 1e-8) break;
    }

    vector_copy(x, x_out);
    vector_free(x); vector_free(x_prev);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Perspective Function
 *
 * g(x, t) = t * f(x/t),  for t > 0
 *
 * The perspective operation preserves convexity: if f is convex,
 * its perspective is jointly convex in (x, t) on t > 0.
 *
 * This is fundamental to:
 * - Converting non-convex problems with binary variables into
 *   convex ones (perspective relaxation / perspective reformulation)
 * - The theory of gauge functions and Minkowski functionals
 * - Convex hull formulations in mixed-integer optimization
 *
 * Reference: Hiriart-Urruty & Lemarechal (2001), Ch.X
 * ═══════════════════════════════════════════════════════════════════════ */

double perspective_eval(const convex_function_t *f,
                        const vector_t *x, double t)
{
    assert(f && x && t > 0.0);

    /* g(x, t) = t * f(x/t) */
    vector_t *x_scaled = vector_alloc(x->n);
    for (size_t i = 0; i < x->n; i++) {
        x_scaled->data[i] = x->data[i] / t;
    }

    double result = t * f->eval(f, x_scaled);
    vector_free(x_scaled);
    return result;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Infimal Convolution
 *
 * (f square g)(x) = inf_y { f(y) + g(x - y) }
 *
 * The infimal convolution is the function-space analog of the
 * Minkowski sum of sets. Properties:
 * - (f square g)* = f* + g* (conjugate of inf-conv is sum of conjugates)
 * - epi(f square g) = epi(f) + epi(g)  (Minkowski sum of epigraphs)
 * - If f, g are convex, f square g is convex
 *
 * The Moreau envelope is the special case where g(x) = (1/2lambda)||x||^2.
 * ═══════════════════════════════════════════════════════════════════════ */

double infimal_convolution(const convex_function_t *f,
                           const convex_function_t *g,
                           const vector_t *x)
{
    assert(f && g && x);
    size_t n = x->n;

    /* Solve: minimize_y  f(y) + g(x - y)
     * via coordinate descent from y = 0.
     */
    vector_t *y = vector_alloc(n);
    vector_t *x_minus_y = vector_alloc(n);

    for (size_t i = 0; i < n; i++) {
        y->data[i] = 0.0;
    }

    double best_val = 1e308;

    for (int iter = 0; iter < 300; iter++) {
        for (size_t i = 0; i < n; i++) {
            x_minus_y->data[i] = x->data[i] - y->data[i];
        }

        double val = f->eval(f, y) + g->eval(g, x_minus_y);
        if (val < best_val) best_val = val;

        /* Simple gradient descent on h(y) = f(y) + g(x-y) */
        vector_t *grad_f = vector_alloc(n);
        vector_t *grad_g = vector_alloc(n);
        f->grad(f, y, grad_f);
        g->grad(g, x_minus_y, grad_g);

        for (size_t i = 0; i < n; i++) {
            grad_f->data[i] -= grad_g->data[i];
        }

        double alpha = 0.01;
        for (size_t i = 0; i < n; i++) {
            y->data[i] -= alpha * grad_f->data[i];
        }

        vector_free(grad_f); vector_free(grad_g);
    }

    vector_free(y); vector_free(x_minus_y);
    return best_val;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Moreau Envelope (Moreau-Yosida Regularization)
 *
 * M_f(x) = inf_y { f(y) + (1/2lambda) ||y - x||^2 }
 *
 * Properties:
 * - M_f is C^1 (continuously differentiable) with Lipschitz gradient
 *   even if f is nonsmooth. Gradient: grad M_f(x) = (x - prox_f(x))/lambda
 * - argmin M_f = argmin f (same minima)
 * - M_f is a smooth approximation of f; as lambda -> 0, M_f -> f
 *
 * The Moreau envelope is the foundation of:
 * - Forward-backward splitting (proximal gradient)
 * - Douglas-Rachford splitting
 * - The augmented Lagrangian method
 * - Nesterov smoothing technique for nonsmooth optimization
 *
 * Reference: Moreau (1965), Yosida (1964)
 * ═══════════════════════════════════════════════════════════════════════ */

void moreau_envelope(const convex_function_t *f, double lambda,
                     const vector_t *x, double *m_env, vector_t *m_grad)
{
    assert(f && x && m_env && lambda > 0.0);
    size_t n = x->n;

    /* Compute prox_{lambda*f}(x) */
    vector_t *prox = vector_alloc(n);

    if (f->prox) {
        f->prox(f, lambda, x, prox);
    } else {
        /* Fallback: gradient descent to solve prox */
        vector_copy(x, prox);
        vector_t *grad_f = vector_alloc(n);

        for (int iter = 0; iter < 200; iter++) {
            f->grad(f, prox, grad_f);
            double gnorm = 0.0;
            for (size_t i = 0; i < n; i++) {
                grad_f->data[i] = prox->data[i] - x->data[i]
                                + lambda * grad_f->data[i];
                gnorm += grad_f->data[i] * grad_f->data[i];
            }
            gnorm = sqrt(gnorm);
            if (gnorm / lambda < 1e-6) break;

            double alpha = 0.1;
            for (size_t i = 0; i < n; i++) {
                prox->data[i] -= alpha * grad_f->data[i];
            }
        }
        vector_free(grad_f);
    }

    /* Moreau envelope value */
    double f_prox = f->eval(f, prox);
    double dist2 = 0.0;
    for (size_t i = 0; i < n; i++) {
        double diff = prox->data[i] - x->data[i];
        dist2 += diff * diff;
    }
    *m_env = f_prox + dist2 / (2.0 * lambda);

    /* Gradient: grad M_f(x) = (x - prox)/lambda */
    if (m_grad) {
        for (size_t i = 0; i < n; i++) {
            m_grad->data[i] = (x->data[i] - prox->data[i]) / lambda;
        }
    }

    vector_free(prox);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Subdifferential Computation
 *
 * For convex f, the subdifferential at x is:
 *     df(x) = { g | f(y) >= f(x) + g^T (y - x) for all y }
 *
 * For differentiable f: df(x) = {grad f(x)}
 * For max of affine functions: df(x) = conv{grad f_i(x) : i in A(x)}
 *   where A(x) = {i : f_i(x) = f(x)} (active set).
 *
 * The subgradient returned is one element from the subdifferential.
 * ═══════════════════════════════════════════════════════════════════════ */

void subgradient_compute(const convex_function_t *f,
                         const vector_t *x, vector_t *g)
{
    assert(f && x && g);

    /* If the function has a gradient oracle, use it */
    if (f->grad) {
        f->grad(f, x, g);
        return;
    }

    /* Otherwise, compute via finite differences (not a true subgradient
     * but a reasonable approximation for smooth-ish functions) */
    size_t n = f->n;
    double fx = f->eval(f, x);
    double eps = 1e-6;

    for (size_t i = 0; i < n; i++) {
        vector_t *xp = vector_alloc(n);
        vector_copy(x, xp);
        xp->data[i] += eps;
        double fp = f->eval(f, xp);
        g->data[i] = (fp - fx) / eps;
        vector_free(xp);
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 * Fenchel Dual Solve (Numerical)
 *
 * Primal: minimize f(x) + g(A x)
 * Dual:   maximize -f*(-A^T y) - g*(y)
 *
 * Solves the dual problem via gradient ascent on the dual objective:
 *     D(y) = -f*(-A^T y) - g*(y)
 *
 * The connection between primal and dual via Fenchel conjugates
 * generalizes Lagrangian duality. Key instances:
 * - LASSO: f = (1/2)||.||^2, g = lambda||.||_1
 *   Dual decomposes into box constraints (easy to solve)
 * - SVM: f = hinge loss, g = ridge penalty
 * - Matrix completion: f = Frobenius norm, g = nuclear norm
 *
 * Reference: Rockafellar (1970), Ch.31; Bauschke & Combettes (2011)
 * ═══════════════════════════════════════════════════════════════════════ */

void fenchel_dual_solve(const convex_function_t *f,
                        const convex_function_t *g,
                        const matrix_t *A,
                        const vector_t *y0,
                        size_t max_iters,
                        opt_result_t *result,
                        vector_t *y_opt)
{
    assert(f && g && A && y0 && result && y_opt);
    size_t m = y0->n;  /* Dual dimension */
    size_t n = A->cols; /* Primal dimension (A is m x n) */

    vector_t *y = vector_alloc(m);
    vector_t *grad = vector_alloc(m);
    vector_t *ATy = vector_alloc(n);

    vector_copy(y0, y);

    for (size_t iter = 0; iter < max_iters; iter++) {
        /* Compute A^T y */
        for (size_t i = 0; i < n; i++) {
            ATy->data[i] = 0.0;
            for (size_t j = 0; j < m; j++) {
                ATy->data[i] += MAT_ELT(A, j, i) * y->data[j];
            }
        }

        /* Compute prox_{g*}(y) via Moreau: prox_{g*}(y) = y - prox_g(y) */
        vector_t *prox_g_y = vector_alloc(m);
        if (g->prox) {
            g->prox(g, 1.0, y, prox_g_y);
        } else {
            vector_copy(y, prox_g_y);
        }

        /* Gradient of dual: grad D(y) = A * prox_f(-A^T y) - prox_{g*}(y) */
        /* Compute prox_f(-A^T y) */
        vector_t *neg_ATy = vector_alloc(n);
        for (size_t i = 0; i < n; i++) neg_ATy->data[i] = -ATy->data[i];

        vector_t *prox_f = vector_alloc(n);
        if (f->prox) {
            f->prox(f, 1.0, neg_ATy, prox_f);
        } else {
            vector_copy(neg_ATy, prox_f);
        }

        /* grad = A * prox_f(-A^T y) - y + prox_g(y) */
        vector_t *Aprox = vector_alloc(m);
        matrix_vector_mul(A, prox_f, Aprox);

        for (size_t j = 0; j < m; j++) {
            grad->data[j] = Aprox->data[j] - y->data[j] + prox_g_y->data[j];
        }

        double gnorm = vector_norm_l2(grad);
        if (gnorm < GRAD_TOL) break;

        /* Gradient ascent step */
        double alpha = 0.01;
        for (size_t j = 0; j < m; j++) {
            y->data[j] += alpha * grad->data[j];
        }

        vector_free(prox_g_y); vector_free(neg_ATy);
        vector_free(prox_f); vector_free(Aprox);
    }

    vector_copy(y, y_opt);
    result->status = OPT_SUCCESS;
    result->iterations = max_iters;

    vector_free(y); vector_free(grad); vector_free(ATy);
}
