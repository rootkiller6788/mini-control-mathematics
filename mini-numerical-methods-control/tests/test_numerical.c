#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include "numerical_core.h"
#include "integration.h"
#include "linear_systems.h"
#include "eigen_methods.h"
#include "optimization.h"
#include "interpolation.h"
#include "control_solvers.h"
#include "nonlinear_methods.h"

#define TOL 1e-10
#define SOFT_TOL 1e-6
static int tests_run = 0, tests_passed = 0, tests_failed = 0;
#define TEST(n) do { tests_run++; printf("  TEST %s... ", #n); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(m) do { tests_failed++; printf("FAIL: %s\n", m); return; } while(0)
#define ASSERT_TRUE(c, m) do { if (!(c)) { FAIL(m); return; } } while(0)
#define ASSERT_NEAR(a, b, tol, m) do { if (fabs((a)-(b)) > (tol)) { \
    printf("FAIL: %s\n", m); tests_failed++; return; } } while(0)

static void test_vector_create(void) {
    TEST(vector_create);
    Vector *v = vector_create(5);
    ASSERT_TRUE(v != NULL, "create failed");
    vector_free(v); PASS();
}
static void test_vector_arithmetic(void) {
    TEST(vector_arithmetic);
    double a[]={1,2,3}, b[]={4,5,6};
    Vector *u=vector_create_from(3,a), *v=vector_create_from(3,b);
    Vector *w=vector_add(u,v);
    ASSERT_NEAR(w->data[0],5.0,TOL,"add[0]");
    ASSERT_NEAR(vector_dot(u,v),32.0,TOL,"dot");
    ASSERT_NEAR(vector_norm_l2(u),sqrt(14.0),TOL,"norm_l2");
    ASSERT_NEAR(vector_norm_l1(u),6.0,TOL,"norm_l1");
    vector_free(u); vector_free(v); vector_free(w); PASS();
}
static void test_vector_cross(void) {
    TEST(vector_cross);
    double a[]={1,0,0}, b[]={0,1,0};
    Vector *u=vector_create_from(3,a), *v=vector_create_from(3,b);
    Vector *w=vector_cross(u,v);
    ASSERT_NEAR(w->data[2],1.0,TOL,"cross z");
    vector_free(u); vector_free(v); vector_free(w); PASS();
}
static void test_vector_axpy(void) {
    TEST(vector_axpy);
    double xd[]={1,2,3}, yd[]={10,20,30};
    Vector *x=vector_create_from(3,xd), *y=vector_create_from(3,yd);
    vector_axpy(2.0,x,y);
    ASSERT_NEAR(y->data[0],12.0,TOL,"axpy[0]");
    vector_free(x); vector_free(y); PASS();
}
static void test_matrix_create(void) {
    TEST(matrix_create);
    Matrix *A=matrix_create(3,4);
    ASSERT_TRUE(A!=NULL && A->rows==3 && A->cols==4,"dims");
    matrix_free(A); PASS();
}
static void test_matrix_multiply(void) {
    TEST(matrix_multiply);
    double Ad[]={1,2,3,4,5,6}, Bd[]={1,2,3,4,5,6};
    Matrix *A=matrix_create_from(3,2,Ad), *B=matrix_create_from(2,3,Bd);
    Matrix *C=matrix_multiply(A,B);
    ASSERT_NEAR(C->data[0],9.0,TOL,"C[0,0]");
    ASSERT_NEAR(C->data[5],33.0,TOL,"C[1,2]");
    matrix_free(A); matrix_free(B); matrix_free(C); PASS();
}
static void test_matrix_transpose(void) {
    TEST(matrix_transpose);
    double Ad[]={1,2,3,4,5,6};
    Matrix *A=matrix_create_from(2,3,Ad);
    Matrix *At=matrix_transpose(A);
    ASSERT_TRUE(At->rows==3 && At->cols==2,"dims");
    matrix_free(A); matrix_free(At); PASS();
}
static void test_matrix_determinant(void) {
    TEST(matrix_determinant);
    double Ad[]={2,1,0,0,1,0,0,0,3};
    Matrix *A=matrix_create_from(3,3,Ad);
    ASSERT_NEAR(matrix_determinant(A),6.0,SOFT_TOL,"det");
    matrix_free(A); PASS();
}
static void test_matrix_identity(void) {
    TEST(matrix_identity);
    Matrix *I=matrix_create_identity(4);
    ASSERT_NEAR(matrix_trace(I),4.0,TOL,"trace");
    matrix_free(I); PASS();
}
static void test_gaussian_elimination(void) {
    TEST(gaussian_elimination);
    double Ad[]={2,1,0,1,3,1,0,1,2}, bd[]={3,5,3};
    Matrix *A=matrix_create_from(3,3,Ad);
    Vector *b=vector_create_from(3,bd), *x=vector_create(3);
    ASSERT_TRUE(gaussian_elimination_solve(A,b,x)==0,"solve");
    ASSERT_NEAR(x->data[0],1.0,SOFT_TOL,"x0=1");
    matrix_free(A); vector_free(b); vector_free(x); PASS();
}
static void test_lu_decomposition(void) {
    TEST(lu_decomposition);
    double Ad[]={4,3,2,1};
    Matrix *A=matrix_create_from(2,2,Ad);
    LUResult *lu=lu_decompose(A);
    ASSERT_TRUE(lu!=NULL&&!lu->singular,"LU");
    ASSERT_NEAR(lu_determinant(lu),-2.0,SOFT_TOL,"det");
    lu_free(lu); matrix_free(A); PASS();
}
static void test_cholesky(void) {
    TEST(cholesky);
    double Ad[]={4,2,2,3};
    Matrix *A=matrix_create_from(2,2,Ad);
    CholeskyResult *chol=cholesky_decompose(A);
    ASSERT_TRUE(chol!=NULL,"Cholesky");
    cholesky_free(chol); matrix_free(A); PASS();
}
static void test_qr_decomposition(void) {
    TEST(qr_decomposition);
    double Ad[]={3,2,1,2};
    Matrix *A=matrix_create_from(2,2,Ad);
    QRResult *qr=qr_decompose_householder(A);
    ASSERT_TRUE(qr!=NULL,"QR");
    Matrix *Qt=matrix_transpose(qr->Q);
    Matrix *QtQ=matrix_multiply(Qt,qr->Q);
    ASSERT_NEAR(QtQ->data[0],1.0,SOFT_TOL,"QtQ");
    matrix_free(Qt); matrix_free(QtQ);
    qr_free(qr); matrix_free(A); PASS();
}
static void test_power_iteration(void) {
    TEST(power_iteration);
    double Ad[]={2,-1,-1,2};
    Matrix *A=matrix_create_from(2,2,Ad);
    Vector *v=vector_create(2); v->data[0]=1; v->data[1]=-0.3;
    double lambda;
    ASSERT_TRUE(power_iteration(A,&lambda,v,100,1e-10)==0,"power iter");
    ASSERT_NEAR(lambda,3.0,SOFT_TOL,"lambda=3");
    matrix_free(A); vector_free(v); PASS();
}
static void test_qr_algorithm(void) {
    TEST(qr_algorithm);
    double Ad[]={2,-1,-1,2};
    Matrix *A=matrix_create_from(2,2,Ad);
    EigenResult *eig=qr_algorithm_eigenvalues(A);
    ASSERT_TRUE(eig!=NULL,"QR algorithm");
    double e1=eig->real_parts[0], e2=eig->real_parts[1];
    ASSERT_TRUE((fabs(e1-3.0)<SOFT_TOL&&fabs(e2-1.0)<SOFT_TOL)||(fabs(e1-1.0)<SOFT_TOL&&fabs(e2-3.0)<SOFT_TOL),"eigvals");
    eigen_result_free(eig); matrix_free(A); PASS();
}
static void test_spectral_radius(void) {
    TEST(spectral_radius);
    double Ad[]={-2,1,0,-3};
    Matrix *A=matrix_create_from(2,2,Ad);
    EigenResult *eig=qr_algorithm_eigenvalues(A);
    SpectralRadius *sr=spectral_radius_compute(eig);
    ASSERT_TRUE(sr->is_hurwitz,"Hurwitz");
    spectral_radius_free(sr); eigen_result_free(eig); matrix_free(A); PASS();
}
static void harmonic_rhs(double t, const double *s, double *d, void *p, size_t n) {
    (void)t;(void)p;(void)n; d[0]=s[1]; d[1]=-s[0];
}
static void test_rk4_single_step(void) {
    TEST(rk4_single_step);
    ODESystem sys={harmonic_rhs,NULL,2};
    double x[]={1,0}, xn[2];
    ASSERT_TRUE(rk4_step(&sys,0,x,0.1,xn)==0,"RK4");
    PASS();
}
static void test_rk4_integration(void) {
    TEST(rk4_integration);
    ODESystem sys={harmonic_rhs,NULL,2};
    double x0[]={1,0};
    ODESolution *sol=solve_ode_fixed_step(&sys,0,2*3.14159265,x0,0.01);
    ASSERT_TRUE(sol!=NULL,"integration");
    size_t last=sol->n_steps-1;
    ASSERT_NEAR(sol->states[last*2],1.0,0.02,"x(2pi)~1");
    odesolution_free(sol); PASS();
}
static double fq(double x,void*p){(void)p; return x-3;}
static double fpq(double x,void*p){(void)p; return 1.0;}
static void test_bisection(void) {
    TEST(bisection); int c;
    ASSERT_NEAR(bisection_method(fq,NULL,0,5,1e-8,100,&c),3.0,SOFT_TOL,"root");
    ASSERT_TRUE(c,"conv"); PASS();
}
static void test_newton_raphson(void) {
    TEST(newton_raphson); int c;
    ASSERT_NEAR(newton_raphson(fq,fpq,NULL,10,1e-8,50,&c),3.0,SOFT_TOL,"root");
    ASSERT_TRUE(c,"conv"); PASS();
}
static double gfv(const double*x,size_t n,void*p){(void)n;(void)p; return (x[0]-3)*(x[0]-3);}
static void gfg(const double*x,size_t n,void*p,double*g){(void)n;(void)p; g[0]=2*(x[0]-3);}
static void test_gradient_descent(void) {
    TEST(gradient_descent);
    double x0[]={0};
    OptResult*r=gradient_descent_backtracking(gfv,gfg,NULL,1,x0,0.1,1e-6,100);
    ASSERT_TRUE(r!=NULL&&r->converged,"GD");
    ASSERT_NEAR(r->x_opt[0],3.0,SOFT_TOL,"xopt");
    opt_result_free(r); PASS();
}
static void test_spline_interpolation(void) {
    TEST(spline_interpolation);
    double xs[]={0,1,2,3}, ys[]={0,1,4,9};
    InterpData d={xs,ys,4};
    Spline *s=spline_cubic_natural(&d);
    ASSERT_TRUE(s!=NULL,"spline");
    ASSERT_NEAR(spline_evaluate(s,1.5),2.25,0.5,"s(1.5)");
    spline_free(s); PASS();
}
static void test_controllability(void) {
    TEST(controllability);
    double Ad[]={0,1,0,0}, Bd[]={0,1};
    Matrix *A=matrix_create_from(2,2,Ad), *B=matrix_create_from(2,1,Bd);
    ASSERT_TRUE(is_controllable(A,B,1e-6),"controllable");
    matrix_free(A); matrix_free(B); PASS();
}
static void test_lqr_design(void) {
    TEST(lqr_design);
    double Ad[]={0,1,0,0}, Bd[]={0,1}, Qd[]={1,0,0,1}, Rd[]={1};
    Matrix *A=matrix_create_from(2,2,Ad), *B=matrix_create_from(2,1,Bd);
    Matrix *Q=matrix_create_from(2,2,Qd), *R=matrix_create_from(1,1,Rd);
    LQRResult *lqr=lqr_design(A,B,Q,R);
    ASSERT_TRUE(lqr!=NULL,"LQR");
    matrix_free(A); matrix_free(B); matrix_free(Q); matrix_free(R);
    lqr_result_free(lqr); PASS();
}
static void test_svd(void) {
    TEST(svd);
    double Ad[]={3,1,1,3};
    Matrix *A=matrix_create_from(2,2,Ad);
    SVDResult *svd=svd_decompose(A);
    ASSERT_TRUE(svd!=NULL,"SVD");
    ASSERT_NEAR(svd->S[0],4.0,0.3,"sigma1~4");
    Matrix *Ap=svd_pseudo_inverse(svd);
    ASSERT_TRUE(Ap!=NULL,"pseudo");
    matrix_free(Ap); svd_free(svd); matrix_free(A); PASS();
}
static void nlr_eval(const double*x,double*F,void*p,size_t n){
    (void)p;(void)n;F[0]=10*(x[1]-x[0]*x[0]);F[1]=1-x[0];
}
static void nlr_jac(const double*x,Matrix*J,void*p,size_t n){
    (void)p;(void)n;J->data[0]=-20*x[0];J->data[1]=10;J->data[2]=-1;J->data[3]=0;
}
static void test_newton_kantorovich(void) {
    TEST(newton_kantorovich);
    NLEquation eq={nlr_eval,nlr_jac,NULL,2};
    double x0[]={0.5,0.5};
    NLSolution*sol=newton_kantorovich(&eq,x0,1e-8,50);
    ASSERT_TRUE(sol!=NULL&&sol->converged,"NK");
    ASSERT_NEAR(sol->x[0],1.0,SOFT_TOL,"x*=1");
    nl_solution_free(sol); PASS();
}
static void test_condition_number(void) {
    TEST(condition_number);
    double Ad[]={1,2,3,4};
    Matrix *A=matrix_create_from(2,2,Ad);
    ConditionNumber*cn=condition_number_estimate(A);
    ASSERT_TRUE(cn!=NULL,"cond");
    condition_number_free(cn); matrix_free(A); PASS();
}
static void test_residual_norm(void) {
    TEST(residual_norm);
    double Ad[]={2,1,1,3},xd[]={1,1},bd[]={3,4};
    Matrix *A=matrix_create_from(2,2,Ad);
    Vector *x=vector_create_from(2,xd),*b=vector_create_from(2,bd);
    ASSERT_NEAR(residual_norm(A,x,b),0.0,SOFT_TOL,"residual");
    matrix_free(A); vector_free(x); vector_free(b); PASS();
}
int main(void) {
    printf("\n=== mini-numerical-methods-control Test Suite ===\n\n");
    test_vector_create(); test_vector_arithmetic(); test_vector_cross();
    test_vector_axpy(); test_matrix_create(); test_matrix_multiply();
    test_matrix_transpose(); test_matrix_determinant(); test_matrix_identity();
    test_condition_number(); test_residual_norm();
    test_gaussian_elimination(); test_lu_decomposition(); test_cholesky();
    test_qr_decomposition(); test_power_iteration(); test_qr_algorithm();
    test_spectral_radius(); test_rk4_single_step(); test_rk4_integration();
    test_bisection(); test_newton_raphson(); test_gradient_descent();
    test_spline_interpolation(); test_controllability(); test_lqr_design();
    test_svd(); test_newton_kantorovich();
    printf("\n=== Results: %d/%d passed, %d failed ===\n",tests_passed,tests_run,tests_failed);
    return tests_failed?1:0;
}
