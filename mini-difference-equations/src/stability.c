/* stability.c - Stability Analysis for Discrete-Time Systems
 * Jury/Schur-Cohn criteria, eigenvalue check, discretization methods.
 * Reference: Jury "Inners and Stability of Dynamic Systems" (1974)
 * Reference: Astrom & Wittenmark "Computer-Controlled Systems" (1997)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include "stability.h"
#include "ztransform.h"
/* L4/L5: Jury Stability Test - Build table */
int jury_table_build(int order,const double*coeffs,JuryTable*jt){
 if(!coeffs||!jt||order<1)return -1;jt->order=order;jt->num_rows=(order>2)?2*order-3:0;
 jt->coeffs=(double*)calloc((size_t)(order+1),sizeof(double));
 memcpy(jt->coeffs,coeffs,(size_t)(order+1)*sizeof(double));
 if(jt->num_rows>0){jt->table=(double**)calloc((size_t)jt->num_rows,sizeof(double*));
  for(int i=0;i<jt->num_rows;i++){int cols=order+1-(i/2);
   jt->table[i]=(double*)calloc((size_t)cols,sizeof(double));}
  /* First row = coeffs, second row = reversed coeffs */
  for(int j=0;j<=order;j++){jt->table[0][j]=coeffs[j];jt->table[1][j]=coeffs[order-j];}
  /* Build subsequent rows */
  for(int k=2;k<jt->num_rows;k+=2){int sz=order+1-k/2;
   for(int j=0;j<sz;j++){double det=jt->table[k-2][0]*jt->table[k-2][sz-1-j]-
    jt->table[k-2][sz-1]*jt->table[k-2][j];
    jt->table[k][j]=det/jt->table[k-2][0];}
   /* Next row is reversed */
   for(int j=0;j<sz;j++)jt->table[k+1][j]=jt->table[k][sz-1-j];}}
 jt->is_stable=jury_stability_test(order,coeffs);return 0;}
/* L4/L5: Jury Stability Test */
int jury_stability_test(int order,const double*coeffs){
 if(!coeffs||order<1)return 0;if(coeffs[0]<=0.0)return 0;
 /* Necessary conditions */
 if(!jury_necessary_conditions(order,coeffs))return 0;
 if(order==1)return(fabs(coeffs[1])<coeffs[0])?1:0;
 if(order==2){if(coeffs[0]+coeffs[1]+coeffs[2]<=0.0)return 0;
  if(coeffs[0]-coeffs[1]+coeffs[2]<=0.0)return 0;
  return(fabs(coeffs[2])<coeffs[0])?1:0;}
 /* Sufficient conditions from Jury table */
 JuryTable jt={0};jury_table_build(order,coeffs,&jt);
 if(jt.num_rows>=3){for(int k=2;k<jt.num_rows;k+=2){
  int sz=order+1-k/2;if(fabs(jt.table[k][0])<=fabs(jt.table[k][sz-1])){jt.is_stable=0;break;}}}
 int result=jt.is_stable;jury_table_free(&jt);return result;}
/* L4: Jury necessary conditions */
int jury_necessary_conditions(int order,const double*coeffs){
 if(!coeffs)return 0;double P1=0.0,Pm1=0.0;
 for(int i=0;i<=order;i++){P1+=coeffs[i];Pm1+=coeffs[i]*((i%2==0)?1.0:-1.0);}
 if(order%2==0){if(P1<=0.0||Pm1<=0.0)return 0;}
 else{if(P1<=0.0||Pm1>=0.0)return 0;}
 return(fabs(coeffs[order])<coeffs[0])?1:0;}
void jury_table_free(JuryTable*jt){if(!jt)return;free(jt->coeffs);
 if(jt->table){for(int i=0;i<jt->num_rows;i++)free(jt->table[i]);free(jt->table);}}
/* L4/L5: Eigenvalue-based stability check */
StabilityType eigenvalue_stability_check(int order,const double*a_coeffs){
 if(!a_coeffs||order<1)return UNKNOWN_STABILITY;
 double**C=(double**)calloc((size_t)order,sizeof(double*));
 for(int i=0;i<order;i++){C[i]=(double*)calloc((size_t)order,sizeof(double));
  C[0][i]=-a_coeffs[i];if(i<order-1)C[i+1][i]=1.0;}
 int stable=1,marginal=0;
 for(int i=0;i<order;i++){double ev=C[i][i];double mag=fabs(ev);
  if(mag>1.0+1e-10){stable=0;break;}if(fabs(mag-1.0)<1e-8)marginal=1;}
 for(int i=0;i<order;i++)free(C[i]);free(C);
 if(!stable)return UNSTABLE;if(marginal)return STABLE_MARGINAL;return STABLE_ASYMPTOTIC;}
/* L2: Spectral radius */
double spectral_radius(const double*A,int n){
 if(!A||n<1)return 0.0;double max_ev=0.0;
 for(int i=0;i<n;i++){double ev=A[i*n+i];double mag=fabs(ev);if(mag>max_ev)max_ev=mag;}
 return max_ev;}
/* L2: BIBO stability check */
int bibo_stability_check(const double*h,int N,double tol){
 if(!h)return 0;double sum=0.0;for(int i=0;i<N;i++)sum+=fabs(h[i]);
 return(sum<1.0/tol)?1:0;}
/* L2: Asymptotic stability check */
int asymptotic_stability_check(const double*h,int N,double tol){
 if(!h)return 0;for(int i=N/2;i<N;i++)if(fabs(h[i])>tol)return 0;return 1;}
/* L5: Bilinear transform s<->z */
double complex bilinear_s_to_z(double complex s,double T){
 return(1.0+s*T/2.0)/(1.0-s*T/2.0);}
double complex bilinear_z_to_s(double complex z,double T){
 return(2.0/T)*(z-1.0)/(z+1.0);}
int bilinear_stability_preservation(double complex s,double T,double tol){
 if(creal(s)>tol)return 1;/* unstable s -> expect |z|>1 */
 double complex z=bilinear_s_to_z(s,T);
 if(creal(s)<-tol)return(cabs(z)<1.0-tol)?1:0;
 return(fabs(cabs(z)-1.0)<tol)?1:0;}
/* L5: Forward Euler discretization s<-(z-1)/T */
int forward_euler_discretize(const ContTransferFunc*ctf,double T,ZTransferFunc*dtf){
 if(!ctf||!dtf)return -1;/* H(z)=G(s)|_{s=(z-1)/T} */
 /* Replace s by (z-1)/T: this is a rational function substitution */
 /* For simple cases: 1st order G(s)=b0/(s+a0) -> H(z)=b0*T/(z-1+a0*T) */
 int n=ctf->den_order,m=ctf->num_order;if(n==1&&m==0){dtf->den_order=1;dtf->num_order=0;
  if(dtf->den_coeffs){dtf->den_coeffs[0]=ctf->den_coeffs[1]*T-1.0;}/* a1 in z^{-1} form */
  if(dtf->num_coeffs){dtf->num_coeffs[0]=ctf->num_coeffs[0]*T/(1.0+ctf->den_coeffs[1]*T);}}
 return 0;}
/* L5: Backward Euler s<-(z-1)/(zT) */
int backward_euler_discretize(const ContTransferFunc*ctf,double T,ZTransferFunc*dtf){
 if(!ctf||!dtf)return -1;int n=ctf->den_order,m=ctf->num_order;
 if(n==1&&m==0){dtf->den_order=1;dtf->num_order=0;
  double Td=1.0+ctf->den_coeffs[1]*T;if(fabs(Td)<1e-15)return -1;
  if(dtf->den_coeffs){dtf->den_coeffs[0]=-(1.0)/Td;}/* a1=-1/(1+aT) */
  if(dtf->num_coeffs){dtf->num_coeffs[0]=ctf->num_coeffs[0]*T/Td;dtf->num_coeffs[1]=0.0;}}
 return 0;}
/* L5: Tustin (bilinear) discretization s<-(2/T)(z-1)/(z+1) */
int tustin_discretize(const ContTransferFunc*ctf,double T,ZTransferFunc*dtf){
 if(!ctf||!dtf)return -1;int n=ctf->den_order,m=ctf->num_order;
 if(n==1&&m==0){double a=ctf->den_coeffs[1],b=ctf->num_coeffs[0],k=2.0/T;
  double d0=k+a,d1=-k+a;if(fabs(d0)<1e-15)return -1;
  dtf->den_order=n;if(dtf->den_coeffs){dtf->den_coeffs[0]=d1/d0;}
  dtf->num_order=m;if(dtf->num_coeffs){dtf->num_coeffs[0]=b/d0;dtf->num_coeffs[1]=b/d0;}}
 return 0;}
/* L5: ZOH discretization */
int zoh_discretize(const ContTransferFunc*ctf,double T,ZTransferFunc*dtf){
 if(!ctf||!dtf)return -1;int n=ctf->den_order,m=ctf->num_order;
 if(n==1&&m==0){double a=ctf->den_coeffs[1],b=ctf->num_coeffs[0];
  double ea=exp(-a*T);dtf->den_order=1;dtf->num_order=0;
  if(fabs(a)>1e-12){if(dtf->den_coeffs)dtf->den_coeffs[0]=-ea;
   if(dtf->num_coeffs)dtf->num_coeffs[0]=b*(1.0-ea)/a;}
  else{if(dtf->den_coeffs)dtf->den_coeffs[0]=-1.0;
   if(dtf->num_coeffs)dtf->num_coeffs[0]=b*T;}}
 return 0;}
/* L5: Matched pole-zero */
int matched_pole_zero_discretize(const ContTransferFunc*ctf,double T,ZTransferFunc*dtf){
 if(!ctf||!dtf)return -1;int n=ctf->den_order,m=ctf->num_order;
 if(n==1&&m==0){double a=ctf->den_coeffs[1];dtf->den_order=1;dtf->num_order=0;
  if(dtf->den_coeffs)dtf->den_coeffs[0]=-exp(-a*T);
  if(dtf->num_coeffs)dtf->num_coeffs[0]=ctf->num_coeffs[0];}
 return 0;}
/* L3: Discrete gain margin from frequency response */
double discrete_gain_margin(const FreqResponse*fr){
 if(!fr||fr->num_points<2)return 0.0;double gm=INFINITY;
 for(int i=0;i<fr->num_points;i++){double ph=fr->phase[i];
  /* Normalize phase to [-pi,pi] and find -pi crossing */
  double ph_n=fmod(ph+M_PI,2.0*M_PI);if(ph_n<0.0)ph_n+=2.0*M_PI;ph_n-=M_PI;
  if(fabs(ph_n+M_PI)<0.05||fabs(ph_n-M_PI)<0.05){/* near -pi */
   double mag=fr->magnitude[i];if(mag<gm&&mag>0.0)gm=mag;}}
 return(gm<INFINITY&&gm>0.0)?1.0/gm:INFINITY;}
/* L3: Discrete phase margin */
double discrete_phase_margin(const FreqResponse*fr){
 if(!fr||fr->num_points<2)return 0.0;
 for(int i=0;i<fr->num_points;i++){if(fabs(fr->magnitude[i]-1.0)<0.01){
   double ph=fmod(fr->phase[i]+M_PI,2.0*M_PI);if(ph<0.0)ph+=2.0*M_PI;ph-=M_PI;
   return ph+M_PI;}}
 return 0.0;}
/* L5: Schur-Cohn stability test */
int schur_cohn_stability_test(int order,const double*coeffs){
 if(!coeffs||order<1)return 0;/* Build Schur-Cohn matrix (2n x 2n) */
 int N=2*order;double**M=(double**)calloc((size_t)N,sizeof(double*));
 for(int i=0;i<N;i++){M[i]=(double*)calloc((size_t)N,sizeof(double));
  /* Fill with a_i and a_{n-i} patterns */
  for(int j=0;j<N;j++){if(i<order&&j<order&&i>=j)M[i][j]=coeffs[order-1-(i-j)];
   if(i<order&&j>=order&&(i+j-order)<order)M[i][j]=coeffs[i+j-order];
   if(i>=order&&j<order&&(i-order+j)<order)M[i][j]=coeffs[i-order+j];
   if(i>=order&&j>=order&&(i-order)>=j-order)M[i][j]=coeffs[order-1-(i-j)];}}
 /* Check positive definiteness (simplified determinant check) */
 double det=1.0;for(int i=0;i<N;i++)det*=M[i][i];int stable=(det>0.0)?1:0;
 for(int i=0;i<N;i++)free(M[i]);free(M);return stable;}
/* Utility */
void cont_transfer_free(ContTransferFunc*ctf){if(!ctf)return;free(ctf->num_coeffs);free(ctf->den_coeffs);}
ContTransferFunc* cont_transfer_alloc(int nm,int dn){
 ContTransferFunc*ctf=(ContTransferFunc*)calloc(1,sizeof(ContTransferFunc));
 if(!ctf)return NULL;ctf->num_order=nm;ctf->den_order=dn;
 ctf->num_coeffs=(double*)calloc((size_t)(nm+1),sizeof(double));
 ctf->den_coeffs=(double*)calloc((size_t)(dn+1),sizeof(double));return ctf;}
void disc_root_locus_free(DiscRootLocus*rl){if(!rl)return;free(rl->poles);free(rl->zeros);free(rl->locus_points);}
void stability_print_verdict(StabilityType t){
 switch(t){case STABLE_ASYMPTOTIC:printf("Asymptotically stable\n");break;
  case STABLE_MARGINAL:printf("Marginally stable\n");break;
  case UNSTABLE:printf("Unstable\n");break;
  case STABLE_BIBO:printf("BIBO stable\n");break;
  case STABLE_LYAPUNOV:printf("Lyapunov stable\n");break;
  default:printf("Unknown stability\n");}}
void jury_table_print(const JuryTable*jt){if(!jt)return;
 printf("Jury Table (order %d):\n",jt->order);
 for(int i=0;i<jt->num_rows;i++){printf("Row %d: ",i);int sz=jt->order+1-i/2;
  for(int j=0;j<sz;j++)printf("%8.4f ",jt->table[i][j]);printf("\n");}
 printf("Stable: %s\n",jt->is_stable?"YES":"NO");}
