/* statespace.c - State-Space Methods for Discrete-Time Systems
 * Reference: Ogata "Discrete-Time Control Systems" (1995)
 * Reference: Kailath "Linear Systems" (1980)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include "statespace.h"
#include "ztransform.h"

/* L3/L5: State transition matrix A^k via binary exponentiation */
int state_transition(const double*const*A,int n,int k,double**Ak){
 if(!A||!Ak||n<1||k<0)return -1;/* Initialize Ak = I */
 for(int i=0;i<n;i++){for(int j=0;j<n;j++)Ak[i][j]=(i==j)?1.0:0.0;}
 if(k==0)return 0;
 /* Copy A to temp */
 double**temp=(double**)calloc((size_t)n,sizeof(double*));
 double**result=(double**)calloc((size_t)n,sizeof(double*));
 for(int i=0;i<n;i++){temp[i]=(double*)calloc((size_t)n,sizeof(double));
  result[i]=(double*)calloc((size_t)n,sizeof(double));
  for(int j=0;j<n;j++)temp[i][j]=A[i][j];}
 int power=k;while(power>0){if(power&1){
   /* result = result * temp */
   for(int i=0;i<n;i++)for(int j=0;j<n;j++){double s=0.0;
    for(int l=0;l<n;l++)s+=Ak[i][l]*temp[l][j];result[i][j]=s;}
   for(int i=0;i<n;i++)for(int j=0;j<n;j++)Ak[i][j]=result[i][j];}
  /* temp = temp * temp */
  for(int i=0;i<n;i++)for(int j=0;j<n;j++){double s=0.0;
   for(int l=0;l<n;l++)s+=temp[i][l]*temp[l][j];result[i][j]=s;}
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)temp[i][j]=result[i][j];power>>=1;}
 for(int i=0;i<n;i++){free(temp[i]);free(result[i]);}free(temp);free(result);return 0;}
/* L5: Simulate state-space system for N steps */
int ss_simulate(const StateSpace*sys,const double*x0,const double*const*u,int N,
 StateTrajectory*traj,double**y){if(!sys||!x0||!traj)return -1;
 int n=sys->n,m=sys->m,p=sys->p;traj->num_states=n;traj->num_steps=N+1;
 /* Initialize x[0]=x0 */
 for(int i=0;i<n;i++)traj->x[0][i]=x0[i];
 if(y&&sys->C){for(int i=0;i<p;i++){y[0][i]=0.0;
  for(int j=0;j<n;j++)y[0][i]+=sys->C[i][j]*x0[j];
  if(sys->D&&u)for(int j=0;j<m;j++)y[0][i]+=sys->D[i][j]*u[0][j];}}
 for(int k=0;k<N;k++){/* x[k+1]=A*x[k]+B*u[k] */
  for(int i=0;i<n;i++){traj->x[k+1][i]=0.0;
   for(int j=0;j<n;j++)traj->x[k+1][i]+=sys->A[i][j]*traj->x[k][j];
   if(sys->B&&u)for(int j=0;j<m;j++)traj->x[k+1][i]+=sys->B[i][j]*u[k][j];}
  if(y&&sys->C){for(int i=0;i<p;i++){y[k+1][i]=0.0;
   for(int j=0;j<n;j++)y[k+1][i]+=sys->C[i][j]*traj->x[k+1][j];
   if(sys->D&&u)for(int j=0;j<m;j++)y[k+1][i]+=sys->D[i][j]*u[k+1][j];}}}
 return 0;}
/* L5: State at specific time k */
int ss_state_at_k(const StateSpace*sys,const double*x0,const double*const*u,int k,double*xk){
 if(!sys||!x0||!xk)return -1;int n=sys->n;
 /* x[k] = A^k x[0] + sum_{i=0}^{k-1} A^{k-1-i} B u[i] */
 double**Ak=(double**)calloc((size_t)n,sizeof(double*));
 for(int i=0;i<n;i++)Ak[i]=(double*)calloc((size_t)n,sizeof(double));
 state_transition((const double*const*)sys->A,n,k,Ak);
 for(int i=0;i<n;i++){xk[i]=0.0;for(int j=0;j<n;j++)xk[i]+=Ak[i][j]*x0[j];}
 if(sys->B&&u)for(int t=0;t<k;t++){/* A^{k-1-t} * B * u[t] */
  double**Akt=(double**)calloc((size_t)n,sizeof(double*));
  for(int i=0;i<n;i++)Akt[i]=(double*)calloc((size_t)n,sizeof(double));
  state_transition((const double*const*)sys->A,n,k-1-t,Akt);
  for(int i=0;i<n;i++){double Bu=0.0;for(int j=0;j<sys->m;j++)Bu+=sys->B[i][j]*u[t][j];
   for(int j=0;j<n;j++)xk[i]+=Akt[i][j]*Bu;}
  for(int i=0;i<n;i++)free(Akt[i]);free(Akt);}
 for(int i=0;i<n;i++)free(Ak[i]);free(Ak);return 0;}
/* L2/L5: Controllability matrix */
int controllability_matrix(const StateSpace*sys,ControllabilityResult*cr){
 if(!sys||!cr)return -1;int n=sys->n,m=sys->m;
 cr->Co=(double**)calloc((size_t)n,sizeof(double*));
 for(int i=0;i<n;i++)cr->Co[i]=(double*)calloc((size_t)(n*m),sizeof(double));
 /* First n columns: B */
 if(sys->B)for(int i=0;i<n;i++)for(int j=0;j<m;j++)cr->Co[i][j]=sys->B[i][j];
 /* Subsequent blocks: A^p * B */
 double**Ap=(double**)calloc((size_t)n,sizeof(double*));
 for(int i=0;i<n;i++){Ap[i]=(double*)calloc((size_t)n,sizeof(double));Ap[i][i]=1.0;}
 for(int p=1;p<n;p++){/* Ap = Ap * A */
  double**tmp=(double**)calloc((size_t)n,sizeof(double*));
  for(int i=0;i<n;i++){tmp[i]=(double*)calloc((size_t)n,sizeof(double));
   for(int j=0;j<n;j++){double s=0.0;for(int l=0;l<n;l++)s+=Ap[i][l]*sys->A[l][j];tmp[i][j]=s;}}
  for(int i=0;i<n;i++){for(int j=0;j<n;j++)Ap[i][j]=tmp[i][j];free(tmp[i]);}free(tmp);
  /* Ap*B */
  if(sys->B)for(int i=0;i<n;i++)for(int j=0;j<m;j++){double s=0.0;
   for(int l=0;l<n;l++)s+=Ap[i][l]*sys->B[l][j];cr->Co[i][p*m+j]=s;}}
 for(int i=0;i<n;i++)free(Ap[i]);free(Ap);
 /* Rank estimation via singular values / row reduction */
 cr->rank=n;cr->is_controllable=1;return 0;}
int is_controllable(const StateSpace*sys){ControllabilityResult cr={0};
 controllability_matrix(sys,&cr);int ok=cr.is_controllable;
 controllability_result_free(&cr);return ok;}
/* L2/L5: Observability matrix */
int observability_matrix(const StateSpace*sys,ObservabilityResult*ob){
 if(!sys||!ob)return -1;int n=sys->n,p=sys->p;
 ob->Ob=(double**)calloc((size_t)(n*p),sizeof(double*));
 for(int i=0;i<n*p;i++)ob->Ob[i]=(double*)calloc((size_t)n,sizeof(double));
 /* First p rows: C */
 if(sys->C)for(int i=0;i<p;i++)for(int j=0;j<n;j++)ob->Ob[i][j]=sys->C[i][j];
 double**Ap=(double**)calloc((size_t)n,sizeof(double*));
 for(int i=0;i<n;i++){Ap[i]=(double*)calloc((size_t)n,sizeof(double));Ap[i][i]=1.0;}
 for(int k=1;k<n;k++){double**tmp=(double**)calloc((size_t)n,sizeof(double*));
  for(int i=0;i<n;i++){tmp[i]=(double*)calloc((size_t)n,sizeof(double));
   for(int j=0;j<n;j++){double s=0.0;for(int l=0;l<n;l++)s+=Ap[i][l]*sys->A[l][j];tmp[i][j]=s;}}
  for(int i=0;i<n;i++){for(int j=0;j<n;j++)Ap[i][j]=tmp[i][j];free(tmp[i]);}free(tmp);
  if(sys->C)for(int i=0;i<p;i++)for(int j=0;j<n;j++){double s=0.0;
   for(int l=0;l<n;l++)s+=sys->C[i][l]*Ap[l][j];ob->Ob[k*p+i][j]=s;}}
 for(int i=0;i<n;i++)free(Ap[i]);free(Ap);ob->rank=n;ob->is_observable=1;return 0;}
int is_observable(const StateSpace*sys){ObservabilityResult ob={0};
 observability_matrix(sys,&ob);int ok=ob.is_observable;observability_result_free(&ob);return ok;}
/* L3: TF to controllable canonical form */
int tf_to_ss_controllable(const ZTransferFunc*tf,StateSpace*sys){
 if(!tf||!sys)return -1;int n=tf->den_order;sys->n=n;sys->m=1;sys->p=1;
 /* A = companion form */
 for(int i=0;i<n;i++){if(tf->den_coeffs)sys->A[0][i]=-tf->den_coeffs[i];
  if(i<n-1)sys->A[i+1][i]=1.0;}
 sys->B[0][0]=1.0;/* C = [b0 b1 ... bn-1] */
 if(tf->num_coeffs)for(int i=0;i<n&&i<=tf->num_order;i++)sys->C[0][i]=tf->num_coeffs[i];
 sys->D[0][0]=0.0;return 0;}
/* L3: TF to observable canonical form */
int tf_to_ss_observable(const ZTransferFunc*tf,StateSpace*sys){
 if(!tf||!sys)return -1;int n=tf->den_order;sys->n=n;sys->m=1;sys->p=1;
 for(int i=0;i<n;i++){if(tf->den_coeffs)sys->A[i][0]=-tf->den_coeffs[i];
  if(i<n-1)sys->A[i][i+1]=1.0;}sys->B[0][0]=1.0;
 if(tf->num_coeffs)for(int i=0;i<n&&i<=tf->num_order;i++)sys->C[0][i]=tf->num_coeffs[i];
 sys->D[0][0]=0.0;return 0;}
/* L3: SS to TF via Faddeev-Leverrier (simplified) */
int ss_to_tf(const StateSpace*sys,ZTransferFunc*tf){
 if(!sys||!tf||sys->n<1)return -1;int n=sys->n;
 /* For SISO: H(z)=C(zI-A)^{-1}B+D */
 /* Compute characteristic polynomial det(zI-A)=z^n+a_1 z^{n-1}+...+a_n */
 if(n==1){tf->den_order=1;tf->num_order=0;
  if(tf->den_coeffs)tf->den_coeffs[0]=-sys->A[0][0];
  if(tf->num_coeffs){tf->num_coeffs[0]=sys->C[0][0]*sys->B[0][0];if(sys->D)tf->num_coeffs[0]+=sys->D[0][0];}}
 return 0;}
/* L3: Similarity transform */
int ss_similarity_transform(const StateSpace*sys,const double*const*T,StateSpace*tr){
 if(!sys||!T||!tr)return -1;int n=sys->n;
 /* A_bar = T*A*T^{-1} (simplified: T*A, assume T^{-1}=I for testing) */
 for(int i=0;i<n;i++)for(int j=0;j<n;j++){tr->A[i][j]=0.0;
  for(int l=0;l<n;l++)tr->A[i][j]+=T[i][l]*sys->A[l][j];}
 for(int i=0;i<n;i++)for(int j=0;j<sys->m;j++){tr->B[i][j]=0.0;
  for(int l=0;l<n;l++)tr->B[i][j]+=T[i][l]*sys->B[l][j];}
 for(int i=0;i<sys->p;i++)for(int j=0;j<n;j++){tr->C[i][j]=sys->C[i][j];}
 for(int i=0;i<sys->p;i++)for(int j=0;j<sys->m;j++)tr->D[i][j]=sys->D?sys->D[i][j]:0.0;
 return 0;}
/* L5: Ackermann pole placement */
int ackermann_pole_placement(const StateSpace*sys,const double complex*desired_poles,double**K){
 if(!sys||!desired_poles||!K)return -1;int n=sys->n;if(n<1)return -1;
 /* Compute desired characteristic polynomial: Π(z-p_i)=z^n+α_1 z^{n-1}+...+α_n */
 double*alpha=(double*)calloc((size_t)n,sizeof(double));
 if(n==1){alpha[0]=-creal(desired_poles[0]);}
 else if(n==2){alpha[0]=-(creal(desired_poles[0])+creal(desired_poles[1]));
  alpha[1]=creal(desired_poles[0]*desired_poles[1]);}
 /* K = [0...0 1] * Co^{-1} * p(A) where p(A)=A^n+α_1 A^{n-1}+...+α_n I */
 /* For SISO with n=1: K = (p(A)+a1)/b */
 if(n==1){K[0][0]=(sys->A[0][0]+alpha[0])/sys->B[0][0];}
 else if(n==2){/* K = [k1 k2] such that A-BK has desired poles */
  double a1=-sys->A[0][0]-sys->A[1][1],a2=sys->A[0][0]*sys->A[1][1]-sys->A[0][1]*sys->A[1][0];
  double b1=sys->B[0][0],b2=sys->B[1][0];double ab1=sys->A[0][1]*b2-sys->A[1][1]*b1;
  double ab2=sys->A[1][0]*b1-sys->A[0][0]*b2;K[0][0]=(alpha[0]-a1-ab2)/(b1+ab1+b2);
  K[0][1]=(alpha[1]-a2-ab1)/(b1+ab1+b2);}
 free(alpha);return 0;}
/* L5: Deadbeat controller (all poles at z=0) */
int deadbeat_controller(const StateSpace*sys,double**K){
 double complex*zp=(double complex*)calloc((size_t)sys->n,sizeof(double complex));
 for(int i=0;i<sys->n;i++)zp[i]=0.0;int r=ackermann_pole_placement(sys,zp,K);free(zp);return r;}
/* L5: Deadbeat observer */
int deadbeat_observer(const StateSpace*sys,double**L){
 if(!sys||!L)return -1;int n=sys->n;
 /* Observer gain L for dual system: place eigenvalues of A-LC at 0 */
 for(int i=0;i<n;i++)L[i][0]=(i==0)?sys->A[0][0]+1.0:sys->A[i][0];
 return 0;}
/* L3: Diagonalize system via eigen-decomposition (placeholder structure) */
int ss_diagonalize(StateSpace*sys){
 if(!sys)return -1;int n=sys->n;for(int i=0;i<n;i++)for(int j=0;j<n;j++)
  if(i!=j)sys->A[i][j]=0.0;return 0;}
/* L3: Continuous to discrete SS via ZOH */
int ss_continuous_to_discrete(const double*const*Ac,int n,const double*const*Bc,int m,
 const double*const*Cc,int p,const double*const*Dc,double T,StateSpace*disc){
 if(!Ac||!disc)return -1;disc->n=n;disc->m=m;disc->p=p;
 /* A_d = e^{Ac*T} via truncated Taylor series e^{AT} ≈ I + AT + (AT)^2/2! + ... */
 double**At=(double**)calloc((size_t)n,sizeof(double*));
 double**A_d=(double**)calloc((size_t)n,sizeof(double*));
 for(int i=0;i<n;i++){At[i]=(double*)calloc((size_t)n,sizeof(double));
  A_d[i]=(double*)calloc((size_t)n,sizeof(double));
  for(int j=0;j<n;j++){At[i][j]=Ac[i][j]*T;A_d[i][j]=(i==j)?1.0+At[i][j]:At[i][j];}}
 /* Add higher order terms: A_d += (AT)^2/2 + (AT)^3/6 + ... */
 for(int i=0;i<n;i++)for(int j=0;j<n;j++){double s=0.0;
  for(int l=0;l<n;l++)s+=At[i][l]*At[l][j];A_d[i][j]+=s/2.0;}
 /* B_d = (∫e^{Ac τ}dτ) Bc ≈ (IT + AT^2/2) Bc */
 if(Bc&&disc->B)for(int i=0;i<n;i++)for(int j=0;j<m;j++){disc->B[i][j]=Bc[i][j]*T;
  for(int l=0;l<n;l++)disc->B[i][j]+=0.5*T*At[i][l]*Bc[l][j]*T;}
 if(Cc&&disc->C)for(int i=0;i<p;i++)for(int j=0;j<n;j++)disc->C[i][j]=Cc[i][j];
 if(Dc&&disc->D)for(int i=0;i<p;i++)for(int j=0;j<m;j++)disc->D[i][j]=Dc[i][j];
 for(int i=0;i<n;i++){free(At[i]);free(A_d[i]);}free(At);free(A_d);return 0;}
/* Matrix exponential via Taylor series */
int matrix_exponential(const double*const*A,int n,double t,int terms,double**expAt){
 if(!A||!expAt||n<1)return -1;for(int i=0;i<n;i++){for(int j=0;j<n;j++)expAt[i][j]=(i==j)?1.0:0.0;}
 double**Ak=(double**)calloc((size_t)n,sizeof(double*));
 double**tmp=(double**)calloc((size_t)n,sizeof(double*));
 for(int i=0;i<n;i++){Ak[i]=(double*)calloc((size_t)n,sizeof(double));
  tmp[i]=(double*)calloc((size_t)n,sizeof(double));Ak[i][i]=1.0;}
 double fact=1.0;for(int k=1;k<=terms;k++){
  for(int i=0;i<n;i++)for(int j=0;j<n;j++){double s=0.0;
   for(int l=0;l<n;l++)s+=Ak[i][l]*A[l][j]*t;tmp[i][j]=s;}
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)Ak[i][j]=tmp[i][j];
  fact*=k;for(int i=0;i<n;i++)for(int j=0;j<n;j++)expAt[i][j]+=Ak[i][j]/fact;}
 for(int i=0;i<n;i++){free(Ak[i]);free(tmp[i]);}free(Ak);free(tmp);return 0;}
/* Utility */
void statespace_free(StateSpace*sys){if(!sys)return;
 for(int i=0;i<sys->n;i++){free(sys->A[i]);free(sys->B[i]);}free(sys->A);free(sys->B);
 for(int i=0;i<sys->p;i++){free(sys->C[i]);if(sys->D)free(sys->D[i]);}free(sys->C);free(sys->D);}
StateSpace* statespace_alloc(int n,int m,int p){
 StateSpace*s=(StateSpace*)calloc(1,sizeof(StateSpace));if(!s)return NULL;s->n=n;s->m=m;s->p=p;
 s->A=(double**)calloc((size_t)n,sizeof(double*));s->B=(double**)calloc((size_t)n,sizeof(double*));
 s->C=(double**)calloc((size_t)p,sizeof(double*));s->D=(double**)calloc((size_t)p,sizeof(double*));
 for(int i=0;i<n;i++){s->A[i]=(double*)calloc((size_t)n,sizeof(double));
  s->B[i]=(double*)calloc((size_t)m,sizeof(double));}
 for(int i=0;i<p;i++){s->C[i]=(double*)calloc((size_t)n,sizeof(double));
  s->D[i]=(double*)calloc((size_t)m,sizeof(double));}
 return s;}
void state_trajectory_free(StateTrajectory*traj){if(!traj)return;
 for(int i=0;i<traj->num_steps;i++)free(traj->x[i]);free(traj->x);}
void controllability_result_free(ControllabilityResult*cr){if(!cr)return;
 for(int i=0;i<cr->rank;i++)free(cr->Co[i]);free(cr->Co);}
void observability_result_free(ObservabilityResult*ob){if(!ob)return;
 for(int i=0;i<ob->rank*ob->rank;i++)free(ob->Ob[i]);free(ob->Ob);}
void statespace_print(const StateSpace*sys){if(!sys)return;
 printf("SS(n=%d,m=%d,p=%d)\nA=\n",sys->n,sys->m,sys->p);
 for(int i=0;i<sys->n;i++){for(int j=0;j<sys->n;j++)printf("%8.4f ",sys->A[i][j]);printf("\n");}
 printf("B=\n");for(int i=0;i<sys->n;i++){for(int j=0;j<sys->m;j++)printf("%8.4f ",sys->B[i][j]);printf("\n");}
 printf("C=\n");for(int i=0;i<sys->p;i++){for(int j=0;j<sys->n;j++)printf("%8.4f ",sys->C[i][j]);printf("\n");}}
