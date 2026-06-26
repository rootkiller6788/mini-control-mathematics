/**
 * @file test_all.c
 * @brief Comprehensive test suite for mini-probability-stochastic.
 * Tests: distributions, statistics, estimation, information theory,
 * stochastic processes, Monte Carlo methods, Kalman/Particle filtering,
 * MCMC, bootstrap, and absorption analysis.
 *
 * Each assertion verifies a specific mathematical property or
 * theoretical relationship.  Written in C99 (no lambdas).
 */
#include "probability_core.h"
#include "stochastic_process.h"
#include "monte_carlo.h"
#include "estimation_filtering.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E  2.71828182845904523536
#endif

#define TOL 1e-10
#define TOL_WEAK 1e-6
#define TOL_STAT 1e-2

static int passed = 0, failed = 0;

#define CHECK(cond, msg) do { \
    if (cond) { passed++; } \
    else { failed++; printf("  FAIL: %s\n", msg); } \
} while(0)

#define CHECK_FLT(got, exp, tol, msg) do { \
    double _g = (got), _e = (exp); \
    if (fabs(_g - _e) <= (tol)) { passed++; } \
    else { failed++; printf("  FAIL: %s (got %.8f exp %.8f)\n", msg, _g, _e); } \
} while(0)

/* Named callback functions for MC/Kalman/IS tests (no C11 lambdas in C99) */
static double cb_sq(double x)        { return x * x; }
static double cb_id(double x)        { return x; }
static double cb_nlogpdf(double x)   { return -0.5 * x * x; }
static double cb_gauss_pdf(double x) { return exp(-0.5*x*x)/sqrt(2.0*M_PI); }
static double cb_cauchy_pdf(double x){ return 1.0/(M_PI*(1.0+x*x)); }
static double cb_cauchy_samp(uint64_t *s) {
    double u = prob_lcg64_uniform(s);
    return tan(M_PI*(u-0.5));
}
static prob_real_t cb_gibbs_s1(prob_real_t x2, uint64_t *s) { (void)x2; return prob_box_muller(s); }
static prob_real_t cb_gibbs_s2(prob_real_t x1, uint64_t *s) { (void)x1; return prob_box_muller(s); }
static prob_real_t cb_gauss_sample(uint64_t *s) { return prob_box_muller(s); }
static double cb_nd_id(const double *x, size_t d) { (void)d; return x[0]*x[0]; }

int main(void) {
    printf("=== mini-probability-stochastic Test Suite ===\n\n");

    /* ------ L1: Discrete Distribution ------ */
    printf("L1: Discrete Distributions\n");
    {
        prob_dpoint_t pts[] = {{1.0, 0.2}, {2.0, 0.3}, {3.0, 0.5}};
        prob_dist_disc_t *d = prob_dist_disc_create(pts, 3);
        CHECK(d != NULL, "create discrete dist");
        CHECK(d->n == 3, "n = 3");
        CHECK_FLT(d->mean, 2.3, TOL, "mean = 2.3");
        CHECK_FLT(d->var, 0.61, TOL_WEAK, "var = 0.61");
        CHECK_FLT(prob_pmf(d, 2.0), 0.3, TOL, "PMF(2)=0.3");
        CHECK_FLT(prob_pmf(d, 99.0), 0.0, TOL, "PMF(99)=0");
        CHECK_FLT(prob_cdf_disc(d, 2.0), 0.5, TOL, "CDF(2)=0.5");
        CHECK_FLT(prob_cdf_disc(d, 3.0), 1.0, TOL, "CDF(3)=1.0");
        CHECK_FLT(prob_mean_disc(d), 2.3, TOL, "mean()=2.3");
        CHECK_FLT(prob_var_disc(d), 0.61, TOL_WEAK, "var()=0.61");
        CHECK_FLT(prob_std_disc(d), sqrt(0.61), TOL_WEAK, "std()");
        prob_dist_disc_free(d);
    }

    /* ------ L1: Continuous Distribution Properties ------ */
    printf("L1: Continuous Distributions\n");
    {
        prob_dist_cont_t g;
        prob_dist_cont_init(&g, PROB_DIST_GAUSSIAN, 0.0, 1.0);
        CHECK_FLT(g.mean, 0.0, TOL, "N(0,1) mean=0");
        CHECK_FLT(g.var, 1.0, TOL, "N(0,1) var=1");
        CHECK_FLT(g.skew, 0.0, TOL, "N(0,1) skew=0");
        CHECK_FLT(g.kurt, 0.0, TOL, "N(0,1) kurt=0");

        prob_t c0 = prob_cdf_cont(&g, 0.0);
        CHECK_FLT(c0, 0.5, TOL_STAT, "Phi(0)=0.5");
        prob_t c2 = prob_cdf_cont(&g, 2.0);
        CHECK(c2 > 0.97 && c2 < 0.98, "Phi(2)~0.977");

        prob_dist_cont_t e;
        prob_dist_cont_init(&e, PROB_DIST_EXPONENTIAL, 1.0, 0.0);
        CHECK_FLT(e.mean, 1.0, TOL, "Exp(1) mean=1");
        CHECK_FLT(e.skew, 2.0, TOL, "Exp(1) skew=2");

        prob_dist_cont_t u;
        prob_dist_cont_init(&u, PROB_DIST_UNIFORM, 0.0, 1.0);
        CHECK_FLT(u.mean, 0.5, TOL, "U(0,1) mean=0.5");
        CHECK_FLT(u.kurt, -1.2, TOL, "U(0,1) kurt=-1.2");
    }

    /* ------ L2: Inverse CDF ------ */
    printf("L2: Inverse CDF\n");
    {
        prob_dist_cont_t g;
        prob_dist_cont_init(&g, PROB_DIST_GAUSSIAN, 0.0, 1.0);
        double q50 = prob_inverse_cdf(&g, 0.5);
        CHECK_FLT(q50, 0.0, TOL_WEAK, "N(0,1) median=0");
        double q975 = prob_inverse_cdf(&g, 0.975);
        CHECK(q975 > 1.9 && q975 < 2.0, "N(0,1) q(0.975)~1.96");

        prob_dist_cont_t e;
        prob_dist_cont_init(&e, PROB_DIST_EXPONENTIAL, 1.0, 0.0);
        double med = prob_inverse_cdf(&e, 0.5);
        CHECK_FLT(med, log(2.0), TOL, "Exp(1) median=ln(2)");
    }

    /* ------ L2: Bayes Theorem ------ */
    printf("L2: Bayes Theorem\n");
    {
        prob_t post = prob_bayes(0.01, 0.95, 0.01*0.95 + 0.99*0.05);
        CHECK(post > 0.15 && post < 0.17, "Bayes: posterior ~ 0.161");
    }

    /* ------ L3: Covariance and Correlation ------ */
    printf("L3: Covariance / Correlation\n");
    {
        double x[] = {1,2,3,4,5}, y[] = {2,4,6,8,10};
        double cov = prob_covariance(x, y, 5);
        CHECK(cov > 0.0, "cov > 0 for y=2x");
        double corr = prob_correlation(x, y, 5);
        CHECK_FLT(corr, 1.0, TOL_WEAK, "corr=1 for y=2x");

        double z[] = {5,4,3,2,1};
        double corr_neg = prob_correlation(x, z, 5);
        CHECK_FLT(corr_neg, -1.0, TOL_WEAK, "corr=-1 for z=6-x");

        /* Covariance matrix */
        double data[] = {1,2,3,4,5, 2,4,6,8,10, 3,6,9,12,15};
        double cmat[9];
        prob_covariance_matrix(data, 3, 5, cmat);
        CHECK(cmat[0] > 0.0, "cov[0,0] > 0 (variance)");
        CHECK(cmat[4] > 0.0, "cov[1,1] > 0");
    }

    /* ------ L3: Skewness and Kurtosis ------ */
    printf("L3: Skewness / Kurtosis\n");
    {
        prob_dpoint_t sym[] = {{-1,0.25},{0,0.5},{1,0.25}};
        prob_dist_disc_t *d = prob_dist_disc_create(sym, 3);
        double sk = prob_skewness_disc(d);
        CHECK_FLT(sk, 0.0, TOL_WEAK, "symmetric dist skew=0");
        prob_dist_disc_free(d);
    }

    /* ------ L4: Total Probability ------ */
    printf("L4: Total Probability\n");
    {
        prob_t part[] = {0.2, 0.3, 0.5};
        double d = prob_total_probability_check(part, 3);
        CHECK_FLT(d, 0.0, TOL_WEAK, "total prob sums to 1");
    }

    /* ------ L4: LLN and CLT ------ */
    printf("L4: LLN / CLT\n");
    {
        double s[] = {0.5,0.6,0.4,0.55,0.45,0.52,0.48,0.51};
        double err = prob_lln_check(s, 8, 0.5);
        CHECK(fabs(err) < 0.05, "LLN error small");

        double z[] = {0.1,-0.2,0.3,-0.1,0.0,0.4,-0.3,0.2};
        double ks = prob_clt_ks_test(z, 8);
        CHECK(ks >= 0.0 && ks <= 1.0, "KS in [0,1]");
    }

    /* ------ L5: Entropy ------ */
    printf("L5: Entropy\n");
    {
        prob_dpoint_t fair[] = {{0,0.5},{1,0.5}};
        prob_dist_disc_t *d = prob_dist_disc_create(fair, 2);
        double h = prob_entropy_disc(d);
        CHECK_FLT(h, 1.0, TOL, "fair coin H=1 bit");
        prob_dist_disc_free(d);

        prob_dpoint_t det[] = {{0,1.0}};
        prob_dist_disc_t *dd = prob_dist_disc_create(det, 1);
        CHECK_FLT(prob_entropy_disc(dd), 0.0, TOL, "deterministic H=0");
        prob_dist_disc_free(dd);
    }

    /* ------ L5: KL and JS Divergence ------ */
    printf("L5: KL / JS Divergence\n");
    {
        prob_dpoint_t p[] = {{0,0.5},{1,0.5}}, q_pts[] = {{0,0.1},{1,0.9}};
        prob_dist_disc_t *dp = prob_dist_disc_create(p, 2);
        prob_dist_disc_t *dq = prob_dist_disc_create(q_pts, 2);
        double kl = prob_kl_divergence(dp, dq);
        CHECK(kl > 0.0, "KL divergence > 0 for different dists");
        double jsd = prob_js_divergence(dp, dq);
        CHECK(jsd > 0.0 && jsd <= 1.0, "JSD in (0,1]");
        prob_dist_disc_free(dp); prob_dist_disc_free(dq);
    }

    /* ------ L5: Mutual Information ------ */
    printf("L5: Mutual Information\n");
    {
        double joint[] = {0.25, 0.25, 0.25, 0.25}; /* independent */
        double px[] = {0.5, 0.5}, py[] = {0.5, 0.5};
        double mi = prob_mutual_information(joint, px, 2, py, 2);
        CHECK_FLT(mi, 0.0, TOL_WEAK, "MI=0 for independent");

        double joint2[] = {0.5, 0.0, 0.0, 0.5}; /* perfectly dependent */
        double mi2 = prob_mutual_information(joint2, px, 2, py, 2);
        CHECK_FLT(mi2, 1.0, TOL_WEAK, "MI=1 bit for perfect dependence");
    }

    /* ------ L5: Differential Entropy ------ */
    printf("L5: Differential Entropy\n");
    {
        prob_dist_cont_t g;
        prob_dist_cont_init(&g, PROB_DIST_GAUSSIAN, 0.0, 1.0);
        double h = prob_differential_entropy(&g);
        double expected = 0.5 * log(2.0 * M_PI * M_E);
        CHECK_FLT(h, expected, TOL_WEAK, "Gaussian diff entropy");
    }

    /* ------ L5: MLE ------ */
    printf("L5: Maximum Likelihood\n");
    {
        double data[] = {2.0,2.1,1.9,2.05,1.95};
        double mu, s2;
        prob_mle_gaussian(data, 5, &mu, &s2);
        CHECK_FLT(mu, 2.0, TOL_WEAK, "Gaussian MLE mu~2.0");

        double edata[] = {1.0, 2.0, 1.5};
        double lam = prob_mle_exponential(edata, 3);
        CHECK(lam > 0.5 && lam < 1.0, "Exp MLE lambda~0.667");

        int32_t pdata[] = {2, 3, 1, 4, 2};
        double plam = prob_mle_poisson(pdata, 5);
        CHECK_FLT(plam, 2.4, TOL, "Poisson MLE=2.4");
    }

    /* ------ L6: KS Test ------ */
    printf("L6: KS Test\n");
    {
        prob_dist_cont_t g;
        prob_dist_cont_init(&g, PROB_DIST_GAUSSIAN, 0.0, 1.0);
        double data[] = {0.1, -0.2, 0.3, 0.5, -0.1};
        double d = prob_ks_test(data, 5, &g);
        CHECK(d >= 0.0 && d <= 1.0, "KS in [0,1]");
    }

    /* ------ L6: Empirical Distribution ------ */
    printf("L6: Empirical Distribution\n");
    {
        prob_empirical_t *e = prob_empirical_create(10);
        CHECK(e != NULL, "create empirical");
        int i; for (i = 1; i <= 10; i++) prob_empirical_insert(e, (double)i);
        CHECK(e->n == 10, "n=10");
        CHECK_FLT(e->mean, 5.5, TOL, "mean=5.5");
        double q50 = prob_empirical_quantile(e, 0.5);
        CHECK(q50 >= 5.0 && q50 <= 6.0, "median~5.5");
        prob_empirical_free(e);
    }

    /* ------ L7: Bootstrap ------ */
    printf("L7: Bootstrap CI\n");
    {
        double data[] = {1,2,3,4,5,6,7,8,9,10};
        double lo, hi;
        prob_bootstrap_ci(data, 10, 1000, 0.95, &lo, &hi);
        CHECK(lo < hi, "bootstrap: lo < hi");
        CHECK(lo < 5.5 && hi > 5.5, "bootstrap CI contains mean");
    }

    /* ------ L1/L2: Markov Chain ------ */
    printf("L1/L2: Markov Chain\n");
    {
        prob_t P[] = {0.5,0.5, 0.3,0.7};
        prob_t pi0[] = {1.0,0.0};
        markov_chain_dtmc_t *mc = markov_dtmc_create(P, pi0, 2);
        CHECK(mc != NULL, "DTMC create");

        prob_t pk[2];
        markov_dtmc_distribution(mc, pk, 1);
        CHECK_FLT(pk[0], 0.5, TOL, "pi_1[0]=0.5");

        prob_t stat[2];
        int r = markov_dtmc_stationary(mc, stat);
        CHECK(r == 0, "stationary solvable");
        CHECK_FLT(stat[0], 0.375, TOL_WEAK, "stat[0]=0.375");
        CHECK_FLT(stat[1], 0.625, TOL_WEAK, "stat[1]=0.625");

        /* Detailed balance */
        double db = markov_dtmc_detailed_balance(mc, stat);
        CHECK_FLT(db, 0.0, TOL_WEAK, "detailed balance holds");

        /* MFPT */
        double mfpt = markov_dtmc_mfpt(mc, 0, 1);
        CHECK(mfpt > 0.0 && isfinite(mfpt), "MFPT(0->1) finite");

        /* Simulate */
        size_t states[10];
        markov_dtmc_simulate(mc, 9, states, 42);
        CHECK(states[0] == 0, "start at state 0");

        /* N-step */
        prob_t P2[4];
        markov_dtmc_nstep(mc, P2, 2);
        CHECK_FLT(P2[0*2+0], 0.4, TOL_WEAK, "P^2[0,0]=0.4");

        markov_dtmc_free(mc);
    }

    /* ------ L1/L2: CTMC ------ */
    printf("L1/L2: CTMC\n");
    {
        prob_real_t Q[] = {-1.0,1.0, 2.0,-2.0};
        markov_chain_ctmc_t *mc = markov_ctmc_create(Q, 2);
        CHECK(mc != NULL, "CTMC create");

        prob_t stat[2];
        int r = markov_ctmc_stationary(mc, stat);
        CHECK(r == 0, "CTMC stationary solvable");
        CHECK_FLT(stat[0], 2.0/3.0, TOL_WEAK, "stat[0]=2/3");
        CHECK_FLT(stat[1], 1.0/3.0, TOL_WEAK, "stat[1]=1/3");

        prob_t Pt[4];
        markov_ctmc_transition(mc, Pt, 1.0);
        CHECK(Pt[0] >= 0.0 && Pt[0] <= 1.0, "P(t) valid");

        markov_ctmc_free(mc);
    }

    /* ------ L1/L2: Poisson Process ------ */
    printf("L1/L2: Poisson Process\n");
    {
        double events[100];
        size_t n = poisson_process_simulate(1.0, 10.0, events, 100, 42);
        CHECK(n > 0, "Poisson: events generated");
        CHECK(n < 100, "Poisson: bounded");
    }

    /* ------ L1/L2: Brownian Motion ------ */
    printf("L1/L2: Brownian Motion\n");
    {
        double W[101];
        brownian_motion_simulate(100, 1.0, W, 42);
        CHECK_FLT(W[0], 0.0, TOL, "W(0)=0");
        CHECK(fabs(W[100]) < 10.0, "W(1) bounded");

        double S[101];
        geometric_brownian_motion(100.0, 0.05, 0.2, 100, 1.0, S, 42);
        CHECK(S[0] == 100.0, "GBM: S(0)=100");

        double X[101];
        ornstein_uhlenbeck_simulate(0.0, 2.0, 0.0, 0.3, 100, 5.0, X, 42);
        CHECK_FLT(X[0], 0.0, TOL, "OU: X(0)=0");
    }

    /* ------ L5: Monte Carlo Integration ------ */
    printf("L5: Monte Carlo Integration\n");
    {
        double se;
        double est = mc_integrate_1d(cb_sq, 0.0, 1.0, 10000, 42, &se);
        CHECK_FLT(est, 1.0/3.0, TOL_STAT, "MC x^2=1/3");
        CHECK(se > 0.0, "MC SE positive");
    }

    /* ------ L5: Metropolis-Hastings ------ */
    printf("L5: Metropolis-Hastings\n");
    {
        double s[5000]; uint64_t seed = 42;
        double acc = mc_metropolis_hastings(cb_nlogpdf, 6000, 1000, 0.0, 1.0, s, &seed);
        CHECK(acc > 0.0 && acc <= 1.0, "MH acc rate valid");
        double m = 0.0; int i; for(i=0;i<5000;i++)m+=s[i]; m/=5000.0;
        CHECK(fabs(m) < 0.15, "MH: sample mean~0");
    }

    /* ------ L5: Importance Sampling ------ */
    printf("L5: Importance Sampling\n");
    {
        /* Estimate E[X^2] for N(0,1) using Cauchy proposal */
        double se;
        double est = mc_importance_sampling(
            cb_sq, cb_gauss_pdf, cb_cauchy_pdf, cb_cauchy_samp,
            5000, 42, &se);
        CHECK(fabs(est - 1.0) < 0.3, "IS: E[X^2]~1 for N(0,1)");
    }

    /* ------ L5: Antithetic Variates ------ */
    printf("L5: Antithetic Variates\n");
    {
        double se;
        double est = mc_antithetic_variates(cb_sq, 0.0, 1.0, 1000, 42, &se);
        CHECK_FLT(est, 1.0/3.0, TOL_STAT, "Antithetic: x^2=1/3");
    }

    /* ------ L5: Control Variates ------ */
    printf("L5: Control Variates\n");
    {
        double se, copt;
        double est = mc_control_variates(cb_sq, cb_id, 0.5, 0.0, 1.0, 1000, 42, &se, &copt);
        CHECK_FLT(est, 1.0/3.0, TOL_STAT, "Control: x^2=1/3");
    }

    /* ------ L5: Gibbs Sampler ------ */
    printf("L5: Gibbs Sampler\n");
    {
        double out1[1000], out2[1000]; uint64_t seed = 42;
        mc_gibbs_sampler_2d(cb_gibbs_s1, cb_gibbs_s2,
            2000, 1000, 0.0, 0.0, out1, out2, &seed);
        double m = 0.0; int i; for(i=0;i<1000;i++)m+=out1[i]; m/=1000.0;
        CHECK(fabs(m) < 0.15, "Gibbs: mean~0");
    }

    /* ------ L5: QMC ------ */
    printf("L5: QMC (Halton)\n");
    {
        double pts[32*2];
        qmc_halton_sequence(2, 32, pts, 1);
        CHECK(pts[0] > 0.0 && pts[0] < 1.0, "Halton point in [0,1]");
    }

    /* ------ L5: Kalman Filter ------ */
    printf("L5: Kalman Filter\n");
    {
        double F[]={1.0},H[]={1.0},Q[]={0.01},R[]={0.1};
        kalman_model_t m={1,0,1,F,NULL,H,Q,R};
        double xi[]={0.0},Pi[]={1.0},obs[]={1,2,3},xe[3],Pe[9];
        kalman_filter_run(&m,obs,NULL,3,xi,Pi,xe,Pe);
        CHECK(xe[0]>0.0,"KF: x_est[0]>0");
        CHECK(xe[2]>2.0,"KF: x_est[2]>2");
    }

    /* ------ L5: Kalman Smoother ------ */
    printf("L5: Kalman Smoother\n");
    {
        double F[]={1.0},H[]={1.0},Q[]={0.01},R[]={0.1};
        kalman_model_t m={1,0,1,F,NULL,H,Q,R};
        double xi[]={0.0},Pi[]={1.0},obs[]={1,2,3},xf[3],Pf[9],xp[3],Pp[9],xs[3],Ps[9];
        /* Run filter manually */
        double x[1]={0.0},P[1]={1.0},xp1[1],Pp1[1];
        int t;
        for(t=0;t<3;t++){
            kalman_predict(&m,x,P,NULL,xp1,Pp1);
            xf[t]=xp1[0];Pf[t]=Pp1[0];xp[t]=xp1[0];Pp[t]=Pp1[0];
            kalman_update(&m,xp1,Pp1,&obs[t],x,P,NULL);
        }
        kalman_smoother(&m,xf,Pf,xp,Pp,3,xs,Ps);
        CHECK(!isnan(xs[0]),"smoother: finite result");
    }

    /* ------ L6: Absorption ------ */
    printf("L6: Absorption Analysis\n");
    {
        prob_t P[]={1,0,0.5,0.5};
        int abs[]={1,0};
        prob_t aprob[4];
        int r=markov_absorption_prob(P,2,abs,aprob);
        CHECK(r==0,"absorption: solvable");
        CHECK_FLT(aprob[0],1.0,TOL,"abs prob from abs state=1");
        CHECK_FLT(aprob[1],1.0,TOL,"abs prob from trans state=1");

        double etime[2];
        r=markov_absorption_time(P,2,abs,etime);
        CHECK(r==0,"absorption time: solvable");
        CHECK_FLT(etime[0],0.0,TOL,"abs time from abs=0");
        CHECK_FLT(etime[1],2.0,TOL,"abs time from trans=2");
    }

    /* ------ L8: MCMC Diagnostics ------ */
    printf("L8: MCMC Diagnostics\n");
    {
        double chains[4*100];
        int j,i;
        for(j=0;j<4;j++) for(i=0;i<100;i++)
            chains[j*100+i] = (double)(j-2) * 0.1;
        double rhat = mcmc_gelman_rubin(chains, 4, 100);
        CHECK(rhat > 1.0, "R-hat > 1 when chains differ");

        double ess = mcmc_effective_sample_size(chains, 100, 20);
        CHECK(ess > 0.0, "ESS positive");
    }

    /* ------ L6: Time Series ACF ------ */
    printf("L6: Time Series ACF\n");
    {
        double x[]={1,2,3,4,5,6,7,8,9,10};
        double acf[6];
        time_series_autocorr(x,10,5,acf);
        CHECK_FLT(acf[0],1.0,TOL,"ACF[0]=1");
        CHECK(acf[1]>0.0,"ACF[1]>0 for trend");
    }

    /* ------ L6: Periodogram ------ */
    printf("L6: Periodogram\n");
    {
        double sig[]={0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1};
        double psd[9];
        periodogram_psd(sig,16,psd);
        CHECK(psd[0]>=0.0,"PSD non-negative");
    }

    /* ------ L2: Distribution Sampling ------ */
    printf("L2: Distribution Sampling\n");
    {
        uint64_t s = 42;
        double u = prob_sample_uniform(0.0, 1.0, &s);
        CHECK(u>0.0&&u<1.0,"sample uniform in (0,1)");
        double g = prob_sample_gaussian(5.0, 2.0, &s);
        CHECK(isfinite(g),"sample gaussian finite");
        double e = prob_sample_exponential(1.0, &s);
        CHECK(e>0.0,"sample exponential>0");
        int b = prob_sample_bernoulli(0.5, &s);
        CHECK(b==0||b==1,"sample Bernoulli 0/1");
        double gam = prob_sample_gamma(2.0, 1.0, &s);
        CHECK(gam>0.0,"sample gamma>0");
        double bet = prob_sample_beta(2.0, 5.0, &s);
        CHECK(bet>0.0&&bet<1.0,"sample beta in (0,1)");
    }

    /* ------ L5: AIC / BIC ------ */
    printf("L5: Model Selection (AIC/BIC)\n");
    {
        double aic = info_criterion_aic(-100.0, 3);
        double bic = info_criterion_bic(-100.0, 3, 50);
        CHECK(bic > aic, "BIC > AIC (stronger penalty)");
    }

    printf("\n=== Results: %d passed, %d failed ===\n", passed, failed);
    return failed ? 1 : 0;
}
