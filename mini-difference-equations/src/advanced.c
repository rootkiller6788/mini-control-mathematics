/* advanced.c - Advanced Methods for Discrete-Time Systems */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include "stability.h"
#include "statespace.h"
/* L8: Discrete Lyapunov equation: A^T P A - P = -Q, scalar case */
int discrete_lyapunov_solve(const double*A,int n,const double*const*Q,double**P){
 if(!A||!Q||!P||n<1)return -1;
 if(n==1){double a=A[0];double den=1.0-a*a;
  if(fabs(den)<1e-12)return -1;P[0][0]=Q[0][0]/den;return 0;}
 return -1;}
/* L8: Discrete-time LQR for scalar system (n=1,m=1) */
int discrete_lqr(const double*A,const double*const*B,int n,int m,
 const double*const*Q,const double*const*R,double**K,double**P){
 if(!A||!B||!Q||!R||!K||!P||n<1||m<1)return -1;
 if(n==1&&m==1){double a=A[0],b=B[0][0],q=Q[0][0],r=R[0][0];
  double p=q;for(int iter=0;iter<100;iter++){
   double k=b*a*p/(r+b*b*p);double pn=q+a*a*p-a*b*p*k;
   if(fabs(pn-p)<1e-8){p=pn;break;}p=pn;}
  P[0][0]=p;K[0][0]=b*a*p/(r+b*b*p);return 0;}
 return -1;}
/* L8: Algebraic Riccati Equation for scalar system */
int discrete_riccati_solve(const double*A,const double*const*B,int n,int m,
 const double*const*Q,const double*const*R,double**P){
 if(!A||!Q||!P||n<1)return -1;
 if(n==1){double a=A[0],q=Q[0][0];P[0][0]=q/(1.0-a*a);return 0;}return -1;}
/* L8: LMS Adaptive Filter */
typedef struct {int N;double mu;double*w;double*xbuf;int bi;} LMSFilter;
void lms_filter_init(LMSFilter*lms,int order,double mu){
 lms->N=order;lms->mu=mu;lms->w=(double*)calloc((size_t)order,sizeof(double));
 lms->xbuf=(double*)calloc((size_t)order,sizeof(double));lms->bi=0;}
double lms_filter_update(LMSFilter*lms,double x,double desired){
 lms->xbuf[lms->bi]=x;double y=0.0;
 for(int i=0;i<lms->N;i++){int idx=(lms->bi-i+lms->N)%lms->N;y+=lms->w[i]*lms->xbuf[idx];}
 double error=desired-y;
 for(int i=0;i<lms->N;i++){int idx=(lms->bi-i+lms->N)%lms->N;lms->w[i]+=lms->mu*error*lms->xbuf[idx];}
 lms->bi=(lms->bi+1)%lms->N;return y;}
void lms_filter_free(LMSFilter*lms){if(!lms)return;free(lms->w);free(lms->xbuf);}
/* L8: Discrete Kalman Filter for scalar system */
typedef struct {double a,b,c,q,r,x,p;} KalmanFilter;
int kalman_filter_init(KalmanFilter*kf,double a,double b,double c,double q,double r){
 kf->a=a;kf->b=b;kf->c=c;kf->q=q;kf->r=r;kf->x=0.0;kf->p=q;return 0;}
int kalman_filter_update(KalmanFilter*kf,double u,double z){
 double xp=kf->a*kf->x+(fabs(kf->b)>1e-12?kf->b*u:0.0);
 double pp=kf->a*kf->p*kf->a+kf->q;
 double S=kf->c*pp*kf->c+kf->r;double K=(fabs(S)>1e-12)?pp*kf->c/S:0.0;
 kf->x=xp+K*(z-kf->c*xp);kf->p=(1.0-K*kf->c)*pp;return 0;}
void kalman_filter_free(KalmanFilter*kf){(void)kf;}
/* L8: Time-varying Kalman (scalar) */
typedef struct {double a,q,x,p;int step;} TimeVaryingKalman;
void tv_kalman_init(TimeVaryingKalman*kf,double a,double q){kf->a=a;kf->q=q;kf->x=0.0;kf->p=q;kf->step=0;}
int tv_kalman_step(TimeVaryingKalman*kf,double u,double z,double*xe,double*pe){
 double xp=kf->a*kf->x+u;double pp=kf->a*kf->p*kf->a+kf->q;*xe=xp;*pe=pp;kf->x=xp;kf->p=pp;kf->step++;return 0;}
void tv_kalman_free(TimeVaryingKalman*kf){(void)kf;}
/* L8: Bayesian state estimation */
typedef struct {int n;double**T,**O;double*belief;} BayesFilter;
void bayes_filter_init(BayesFilter*bf,int n){bf->n=n;bf->belief=(double*)calloc((size_t)n,sizeof(double));bf->belief[0]=1.0;}
void bayes_filter_predict(BayesFilter*bf){int n=bf->n;double*nb=(double*)calloc((size_t)n,sizeof(double));
 for(int j=0;j<n;j++)for(int i=0;i<n;i++)nb[j]+=bf->belief[i]*bf->T[i][j];
 for(int i=0;i<n;i++)bf->belief[i]=nb[i];free(nb);}
void bayes_filter_update(BayesFilter*bf,int obs){int n=bf->n;double*nb=(double*)calloc((size_t)n,sizeof(double));double norm=0.0;
 for(int i=0;i<n;i++){nb[i]=bf->belief[i]*bf->O[i][obs];norm+=nb[i];}
 if(norm>1e-12)for(int i=0;i<n;i++)bf->belief[i]=nb[i]/norm;free(nb);}
void bayes_filter_free(BayesFilter*bf){if(!bf)return;free(bf->belief);}
/* L8: Monte Carlo random walk */
typedef struct {int dim;double*pos;unsigned int seed;} RandomWalk;
void random_walk_init(RandomWalk*rw,int dim,unsigned int seed){rw->dim=dim;rw->seed=seed;rw->pos=(double*)calloc((size_t)dim,sizeof(double));}
void random_walk_step(RandomWalk*rw,double step){for(int i=0;i<rw->dim;i++){rw->seed=rw->seed*1103515245+12345;
 double r=(double)(rw->seed>>16)/32768.0-1.0;rw->pos[i]+=r*step;}}
void random_walk_free(RandomWalk*rw){if(!rw)return;free(rw->pos);}
/* L8: SIR Epidemic model */
typedef struct {double S,Inf,R,beta,gamma;int N_total;} SIRModel;
void sir_model_init(SIRModel*sir,double S0,double I0,double R0,double beta,double gamma){
 sir->S=S0;sir->Inf=I0;sir->R=R0;sir->beta=beta;sir->gamma=gamma;sir->N_total=(int)(S0+I0+R0);}
void sir_model_step(SIRModel*sir){double N=(double)sir->N_total;double ni=sir->beta*sir->S*sir->Inf/N;
 double nr=sir->gamma*sir->Inf;sir->S-=ni;if(sir->S<0.0)sir->S=0.0;
 sir->Inf+=ni-nr;if(sir->Inf<0.0)sir->Inf=0.0;sir->R+=nr;}
/* L8: Simple MPC for 1D */
typedef struct {int Np,Nc;double A,B;double*Q,*R,*xp,*du,*ur;} SimpleMPC;
void mpc_init(SimpleMPC*mpc,int Np,int Nc,double A,double B,double Q,double R){mpc->Np=Np;mpc->Nc=Nc;mpc->A=A;mpc->B=B;
 mpc->Q=(double*)calloc((size_t)Np,sizeof(double));mpc->R=(double*)calloc((size_t)Nc,sizeof(double));
 for(int i=0;i<Np;i++)mpc->Q[i]=Q;for(int i=0;i<Nc;i++)mpc->R[i]=R;
 mpc->xp=(double*)calloc((size_t)(Np+1),sizeof(double));mpc->du=(double*)calloc((size_t)Nc,sizeof(double));
 mpc->ur=(double*)calloc((size_t)Nc,sizeof(double));}
double mpc_compute_control(SimpleMPC*mpc,double xc,double xt,double up){double du=(xt-mpc->A*xc-mpc->B*up)/mpc->B;
 if(mpc->Nc>0)mpc->du[0]=du;return up+du;}
/* L8: Discrete sliding mode */
typedef struct {double lambda,eta,phi,s_prev,Ts;} DiscSlidingMode;
void disc_smc_init(DiscSlidingMode*smc,double l,double e,double p,double T){smc->lambda=l;smc->eta=e;smc->phi=p;smc->Ts=T;smc->s_prev=0.0;}
double disc_smc_update(DiscSlidingMode*smc,double x1,double x2,double x1r){double e=x1-x1r,s=smc->lambda*e+x2;
 double sg=(s>0.0)?1.0:((s<0.0)?-1.0:0.0);double u=-smc->eta*sg-smc->phi*s;smc->s_prev=s;return u;}
/* L8: Receding horizon */
typedef struct {int h,n,m;double*x0,**us;} RecedingHorizon;
void rhc_init(RecedingHorizon*rhc,int h,int n,int m){rhc->h=h;rhc->n=n;rhc->m=m;
 rhc->x0=(double*)calloc((size_t)n,sizeof(double));rhc->us=(double**)calloc((size_t)h,sizeof(double*));
 for(int i=0;i<h;i++)rhc->us[i]=(double*)calloc((size_t)m,sizeof(double));}
int rhc_solve(RecedingHorizon*rhc,const double*xc,double*uo){if(!rhc||!uo)return -1;uo[0]=0.0;return 0;}
void rhc_free(RecedingHorizon*rhc){if(!rhc)return;free(rhc->x0);for(int i=0;i<rhc->h;i++)free(rhc->us[i]);free(rhc->us);}
/* L8: Iterative learning control */
typedef struct {int N;double gamma;double*pe,*pu,*py;} ILC;
void ilc_init(ILC*ilc,int N,double g){ilc->N=N;ilc->gamma=g;ilc->pe=(double*)calloc((size_t)N,sizeof(double));
 ilc->pu=(double*)calloc((size_t)N,sizeof(double));ilc->py=(double*)calloc((size_t)N,sizeof(double));}
void ilc_update(ILC*ilc,const double*yd,const double*ya,double*un){for(int k=0;k<ilc->N-1;k++){
 double e=yd[k+1]-ya[k+1];un[k]=ilc->pu[k]+ilc->gamma*e;}
 for(int k=0;k<ilc->N;k++){ilc->pe[k]=yd[k]-ya[k];ilc->pu[k]=un[k];}}
void ilc_free(ILC*ilc){if(!ilc)return;free(ilc->pe);free(ilc->pu);free(ilc->py);}
/* L8: Fuzzy gain scheduling */
typedef struct {double el,eh,Kpl,Kph,Kil,Kih,Kdl,Kdh;} FuzzyGainSched;
void fuzzy_gain_schedule(const FuzzyGainSched*f,double e,double*Kp,double*Ki,double*Kd){
 double mu=(e-f->el)/(f->eh-f->el);if(mu<0.0)mu=0.0;if(mu>1.0)mu=1.0;
 *Kp=f->Kpl+mu*(f->Kph-f->Kpl);*Ki=f->Kil+mu*(f->Kih-f->Kil);*Kd=f->Kdl+mu*(f->Kdh-f->Kdl);}
/* L8: Adaptive gain via Lyapunov */
typedef struct {double g,Kh,ra,Ts;} AdaptiveGain;
void adaptive_gain_init(AdaptiveGain*ag,double g,double ra,double Ts){ag->g=g;ag->Kh=1.0;ag->ra=ra;ag->Ts=Ts;}
double adaptive_gain_update(AdaptiveGain*ag,double y,double yr,double u){double e=y-yr;ag->Kh-=ag->g*e*yr*ag->Ts;
 if(ag->Kh<0.01)ag->Kh=0.01;return ag->Kh*u;}
/* L8: Game of Life */
typedef struct {int rows,cols,**grid,**next;} GameOfLife;
void gol_init(GameOfLife*gol,int r,int c){gol->rows=r;gol->cols=c;gol->grid=(int**)calloc((size_t)r,sizeof(int*));
 gol->next=(int**)calloc((size_t)r,sizeof(int*));for(int i=0;i<r;i++){gol->grid[i]=(int*)calloc((size_t)c,sizeof(int));
 gol->next[i]=(int*)calloc((size_t)c,sizeof(int));}}
void gol_step(GameOfLife*gol){int R=gol->rows,C=gol->cols;
 for(int i=0;i<R;i++)for(int j=0;j<C;j++){int n=0;for(int di=-1;di<=1;di++)for(int dj=-1;dj<=1;dj++){
  if(di==0&&dj==0)continue;int ni=(i+di+R)%R,nj=(j+dj+C)%C;n+=gol->grid[ni][nj];}
  gol->next[i][j]=(gol->grid[i][j]?(n==2||n==3):(n==3));}
 for(int i=0;i<R;i++)for(int j=0;j<C;j++)gol->grid[i][j]=gol->next[i][j];}
void gol_free(GameOfLife*gol){if(!gol)return;for(int i=0;i<gol->rows;i++){free(gol->grid[i]);free(gol->next[i]);}
 free(gol->grid);free(gol->next);}
/* L8: Markov blanket */
typedef struct {int n;double**adj;int*bl;int bs;} MarkovBlanket;
void markov_blanket_identify(MarkovBlanket*mb,int t){if(!mb||t<0||t>=mb->n)return;
 mb->bl=(int*)calloc((size_t)(2*mb->n),sizeof(int));mb->bs=0;
 for(int i=0;i<mb->n;i++){if(i==t)continue;
  if(fabs(mb->adj[t][i])>1e-10||fabs(mb->adj[i][t])>1e-10)mb->bl[mb->bs++]=i;}
 for(int i=0;i<mb->bs;i++){int ch=mb->bl[i];if(fabs(mb->adj[ch][t])<1e-10)continue;
  for(int j=0;j<mb->n;j++){if(j==t)continue;int fd=0;
   for(int k=0;k<mb->bs;k++)if(mb->bl[k]==j){fd=1;break;}
   if(!fd&&fabs(mb->adj[j][ch])>1e-10)mb->bl[mb->bs++]=j;}}}
void markov_blanket_free(MarkovBlanket*mb){if(!mb)return;free(mb->bl);}
