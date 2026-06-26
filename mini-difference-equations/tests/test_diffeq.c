#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <complex.h>
#include "diffeq.h"
#include "ztransform.h"
#include "stability.h"
#include "statespace.h"
#include "signal.h"
#define TOL 1e-10
static double fcos_h(double x){(void)x;return cos(x);}
static double sqr_h(double x){return x*x;}
int main(void){
 int passed=0;
 printf("=== mini-difference-equations Test Suite ===\n\n");
 printf("  TEST: forward_difference... ");
 {DiscreteSignal*s=discrete_signal_alloc(5,0);s->values[0]=1;s->values[1]=2;s->values[2]=4;s->values[3]=7;s->values[4]=11;
  double d=forward_difference(s,2);assert(fabs(d-3.0)<TOL);discrete_signal_free(s);free(s);}
 printf("PASS\n");passed++;
 printf("  TEST: backward_difference... ");
 {DiscreteSignal*s=discrete_signal_alloc(5,0);s->values[0]=1;s->values[1]=3;s->values[2]=6;s->values[3]=10;s->values[4]=15;
  double d=backward_difference(s,3);assert(fabs(d-4.0)<TOL);discrete_signal_free(s);free(s);}
 printf("PASS\n");passed++;
 printf("  TEST: classify_difference_eq... ");
 {double a[]={0.5,-0.3},b[]={0.0,0.0};
  assert(classify_difference_eq(2,a,b)==DEQ_TYPE_HOMOGENEOUS);}
 printf("PASS\n");passed++;
 printf("  TEST: solve_iterative_lindiffeq... ");
 {LinearDiffEq*eq=linear_diffeq_alloc(1,0);double a[]={0.0};double b[]={1.0};
  linear_diffeq_set_coeffs(eq,a,b);double x[]={0,1,2,3,4};
  double iy[]={0.0},y[5];solve_iterative_lindiffeq(eq,x,iy,5,y);
  assert(fabs(y[0])<TOL);assert(fabs(y[1]-1.0)<TOL);assert(fabs(y[2]-2.0)<TOL);
  linear_diffeq_free(eq);free(eq);}
 printf("PASS\n");passed++;
 printf("  TEST: compute_impulse_response... ");
 {LinearDiffEq*eq=linear_diffeq_alloc(1,1);double a[]={0.0};double b[]={0.0,1.0};
  linear_diffeq_set_coeffs(eq,a,b);ImpulseResponse ir;
  ir.length=10;ir.h=(double*)calloc(10,sizeof(double));
  compute_impulse_response(eq,10,&ir);assert(fabs(ir.h[1]-1.0)<TOL);
  free(ir.h);linear_diffeq_free(eq);free(eq);}
 printf("PASS\n");passed++;
 printf("  TEST: check_existence_uniqueness... ");
 {LinearDiffEq*eq=linear_diffeq_alloc(2,0);double a[]={0.5,-0.3};
  linear_diffeq_set_coeffs(eq,a,NULL);double ic[]={1.0,0.5};
  ExistenceResult er=check_existence_uniqueness(eq,ic,2);
  assert(er.exists==1&&er.is_unique==1);linear_diffeq_free(eq);free(eq);}
 printf("PASS\n");passed++;
 printf("  TEST: time_constant_discrete... ");
 {CharPoly cp;cp.num_roots=1;cp.roots=(double complex*)calloc(1,sizeof(double complex));
  cp.multiplicity=(int*)calloc(1,sizeof(int));cp.roots[0]=0.5+0.0*I;
  double tau=time_constant_discrete(&cp);assert(tau>0.0&&tau<INFINITY);char_poly_free(&cp);}
 printf("PASS\n");passed++;
 printf("  TEST: jury_stability_test... ");
 {double coeffs[]={1.0,-0.5,0.06};assert(jury_stability_test(2,coeffs)==1);}
 printf("PASS\n");passed++;
 printf("  TEST: bilinear_stability_preservation... ");
 {double complex s=-2.0+1.0*I;assert(bilinear_stability_preservation(s,0.1,1e-8)==1);}
 printf("PASS\n");passed++;
 printf("  TEST: ztransform_at... ");
 {ZSequence*seq=zsequence_alloc(5);seq->x[0]=1;seq->x[1]=2;seq->x[2]=3;seq->x[3]=4;seq->x[4]=5;
  double complex z=2.0+0.0*I,Xz=ztransform_at(seq,z);
  assert(fabs(creal(Xz)-3.5625)<1e-8);zsequence_free(seq);free(seq);}
 printf("PASS\n");passed++;
 printf("  TEST: inverse_ztransform_long_division... ");
 {ZTransferFunc*tf=ztransfer_alloc(0,1);tf->den_coeffs[0]=-0.5;tf->num_coeffs[0]=1.0;
  double y[10];inverse_ztransform_long_division(tf,10,y);
  assert(fabs(y[0]-1.0)<TOL);assert(fabs(y[2]-0.25)<TOL);
  ztransfer_free(tf);free(tf);}
 printf("PASS\n");passed++;
 printf("  TEST: convolve_discrete... ");
 {double x[]={1,2,3},h[]={1,1},y[4];int Ny=convolve_discrete(x,3,h,2,y);
  assert(Ny==4&&fabs(y[1]-3.0)<TOL);}
 printf("PASS\n");passed++;
 printf("  TEST: discrete_lyapunov_solve... ");
 {double A[]={0.5};double Qd=1.0,Pd=0.0;double*Qa[1]={&Qd};double*Pa[1]={&Pd};
  assert(discrete_lyapunov_solve(A,1,(const double*const*)Qa,Pa)==0);
  assert(fabs(Pa[0][0]-4.0/3.0)<TOL);}
 printf("PASS\n");passed++;
 printf("  TEST: fibonacci_closed_form... ");
 {double F10;fibonacci_closed_form(10,&F10);assert(fabs(F10-55.0)<TOL);}
 printf("PASS\n");passed++;
 printf("  TEST: signal_sinusoid... ");
 {double sig[100];signal_sinusoid(sig,100,1.0,0.1*M_PI,0.0);
  assert(fabs(sig[5]-cos(0.5*M_PI))<TOL);}
 printf("PASS\n");passed++;
 printf("  TEST: signal_energy... ");
 {double sig[]={1,2,3};assert(fabs(signal_energy(sig,3)-14.0)<TOL);}
 printf("PASS\n");passed++;
 printf("  TEST: find_fixed_points... ");
 {FixedPoint fps[5];int n=find_fixed_points(fcos_h,0.0,2.0,1e-6,fps,5);
  assert(n>0&&fps[0].type==FIXED_STABLE);}
 printf("PASS\n");passed++;
 printf("  TEST: nonlinear_map_iterate... ");
 {double traj[10];nonlinear_map_iterate(sqr_h,2.0,10,traj);assert(fabs(traj[2]-16.0)<TOL);}
 printf("PASS\n");passed++;
 printf("  TEST: malthusian_growth... ");
 {double pop[10];malthusian_growth(100.0,0.02,10,pop);assert(fabs(pop[1]-102.0)<TOL);}
 printf("PASS\n");passed++;
 printf("  TEST: logistic_map_iterate... ");
 {double traj[20];logistic_map_iterate(2.0,0.2,20,traj);assert(fabs(traj[0]-0.2)<TOL);}
 printf("PASS\n");passed++;
 printf("  TEST: moving_average_filter... ");
 {double x[]={1,2,3,4,5,6,7,8,9,10},y[10];
  moving_average_filter(x,10,3,y);assert(fabs(y[2]-2.0)<TOL);}
 printf("PASS\n");passed++;
 printf("\n=== All %d tests passed ===\n",passed);
 return 0;
}
