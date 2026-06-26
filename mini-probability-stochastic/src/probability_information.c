/**
 * @file probability_information.c
 * @brief Information-theoretic measures: Shannon entropy, differential
 *        entropy, Kullback-Leibler divergence, Jensen-Shannon divergence,
 *        mutual information.
 *
 * Knowledge: L5 engineering methods (information theory in probability).
 * Reference: Cover & Thomas, "Elements of Information Theory", 2nd ed.
 */
#include "probability_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E  2.71828182845904523536
#endif

static double sq(double x) { return x * x; }

/*
 * prob_entropy_disc — Shannon entropy H(X) = -sum p_i * log2(p_i) (bits).
 * H_max = log2(n) for uniform distribution. H_min = 0 for deterministic.
 * Knowledge: entropy is the fundamental limit of lossless compression
 * (Shannon, 1948). It quantifies uncertainty: the minimum expected
 * number of binary questions to identify X.
 * O(n) in number of mass points. Independent knowledge point.
 */
prob_real_t prob_entropy_disc(const prob_dist_disc_t *dist) {
    if (!dist || dist->n == 0) return 0.0;
    double h = 0.0;
    size_t i;
    for (i = 0; i < dist->n; i++) {
        double p = dist->points[i].pmass;
        if (p > 0.0) h -= p * log2(p);
    }
    return h;
}

/*
 * prob_differential_entropy — h(X) = -integral f(x) * ln(f(x)) dx (nats).
 * Supports Gaussian, exponential, uniform, gamma with analytic formulas.
 * Knowledge: differential entropy extends Shannon entropy to continuous
 * distributions but lacks some desirable properties (can be negative,
 * not invariant under nonlinear transformations — unlike discrete entropy).
 * Gaussian: h = 0.5 * ln(2*pi*e*sigma^2), maximum among all continuous
 * distributions with fixed variance (max-entropy principle).
 * Exponential: h = 1 - ln(lambda).
 * Uniform: h = ln(b - a).
 * Gamma: h = alpha + ln(beta*Gamma(alpha)) + (1-alpha)*psi(alpha),
 * where psi is the digamma function.
 * O(1). Independent knowledge point: continuous entropy.
 */
prob_real_t prob_differential_entropy(const prob_dist_cont_t *dist) {
    if (!dist) return NAN;
    switch (dist->type) {
    case PROB_DIST_GAUSSIAN:
        return 0.5 * log(2.0 * M_PI * M_E * dist->var);
    case PROB_DIST_EXPONENTIAL:
        return 1.0 - log(dist->param1);
    case PROB_DIST_UNIFORM:
        return log(dist->param2 - dist->param1);
    case PROB_DIST_GAMMA: {
        double a = dist->param1, b = dist->param2;
        /* Digamma approx: psi(x) ~ ln(x) - 1/(2x) - 1/(12x^2) */
        double psi_a = log(a) - 0.5 / a - 1.0 / (12.0 * a * a);
        return a + log(b * tgamma(a)) + (1.0 - a) * psi_a;
    }
    default:
        return NAN;
    }
}

/*
 * prob_kl_divergence — D_KL(P || Q) = sum_i p_i * log2(p_i / q_i) (bits).
 * Properties: D_KL >= 0 (Gibbs inequality), = 0 iff P = Q almost everywhere.
 * Non-symmetric: D_KL(P||Q) != D_KL(Q||P) in general.
 * Returns INFINITY if Q has zero mass where P has positive mass
 * (absolute continuity violation).
 * Knowledge: KL divergence is the information lost when Q is used to
 * approximate P. Not a metric but a premetric (f-divergence).
 * O(n^2) naive support matching. Independent knowledge point.
 */
prob_real_t prob_kl_divergence(const prob_dist_disc_t *p,
                               const prob_dist_disc_t *q) {
    if (!p || !q) return NAN;
    double kl = 0.0;
    size_t i, j;
    for (i = 0; i < p->n; i++) {
        double p_i = p->points[i].pmass;
        if (p_i <= 0.0) continue;
        double q_i = 0.0;
        for (j = 0; j < q->n; j++) {
            if (fabs(q->points[j].value - p->points[i].value) < 1e-12) {
                q_i = q->points[j].pmass; break;
            }
        }
        if (q_i <= 0.0) return INFINITY;
        kl += p_i * log2(p_i / q_i);
    }
    return kl;
}

/*
 * prob_js_divergence — Jensen-Shannon divergence.
 * JSD(P||Q) = 0.5 * D_KL(P||M) + 0.5 * D_KL(Q||M)
 * where M = 0.5 * (P + Q) is the midpoint distribution.
 * Bounded in [0, 1] when using base-2 logarithm.
 * sqrt(JSD) satisfies the triangle inequality (Endres & Schindelin, 2003),
 * making it a true metric on the space of probability distributions.
 * Assumes P and Q are defined on identical support points.
 * Knowledge: JSD is the canonical symmetric, smoothed version of KL
 * divergence. Used extensively in GANs and NLP (e.g., BLEU scores).
 * O(n). Independent knowledge point: symmetric information divergence.
 */
prob_real_t prob_js_divergence(const prob_dist_disc_t *p,
                               const prob_dist_disc_t *q) {
    if (!p || !q || p->n != q->n) return NAN;
    size_t n = p->n;
    prob_dpoint_t *m_pts = (prob_dpoint_t *)malloc(
        n * sizeof(prob_dpoint_t));
    if (!m_pts) return NAN;
    size_t i;
    for (i = 0; i < n; i++) {
        m_pts[i].value = p->points[i].value;
        m_pts[i].pmass = 0.5 * (p->points[i].pmass + q->points[i].pmass);
    }
    prob_dist_disc_t m_dist;
    m_dist.points = m_pts; m_dist.n = n;
    m_dist.mean = 0.0; m_dist.var = 0.0;
    double jsd = 0.5 * prob_kl_divergence(p, &m_dist)
               + 0.5 * prob_kl_divergence(q, &m_dist);
    free(m_pts);
    return jsd;
}

/*
 * prob_mutual_information — I(X;Y) = sum_{x,y} P(x,y) * log2(P(x,y) / (P(x)*P(y))).
 * Properties: I(X;Y) >= 0, = 0 iff X and Y are independent.
 * I(X;Y) = H(X) - H(X|Y) = H(Y) - H(Y|X) = H(X) + H(Y) - H(X,Y).
 * Knowledge: mutual information measures the amount of information shared
 * between X and Y. It is symmetric: I(X;Y) = I(Y;X).
 * O(nx * ny). Independent knowledge point: dependence quantification
 * beyond linear correlation.
 */
prob_real_t prob_mutual_information(const prob_real_t *joint,
                                    const prob_real_t *px, size_t nx,
                                    const prob_real_t *py, size_t ny) {
    if (!joint || !px || !py || nx == 0 || ny == 0) return NAN;
    double mi = 0.0;
    size_t i, j;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            double pxy = joint[i * ny + j];
            if (pxy <= 0.0) continue;
            double px_i = px[i], py_j = py[j];
            if (px_i <= 0.0 || py_j <= 0.0) continue;
            mi += pxy * log2(pxy / (px_i * py_j));
        }
    }
    return mi;
}
