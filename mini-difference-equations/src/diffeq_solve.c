/* diffeq_solve.c - Solution Methods for Difference Equations
 * Reference: Elaydi "An Introduction to Difference Equations" (2005)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include "diffeq.h"
/* L1: Forward/Backward/Shift operators */
double forward_difference(const DiscreteSignal*s,int n){
 if(!s||!s->values)return 0.0;int i=n-s->start_index;
 if(i<0||i+1>=s->length)return 0.0;return s->values[i+1]-s->values[i];}
void forward_difference_signal(const DiscreteSignal*in,DiscreteSignal*out){
 if(!in||!out||!in->values||!out->values)return;int N=in->length;
 if(N<2){if(out->length>0)out->values[0]=0.0;return;}
 for(int i=0;i<N-1&&i<out->length;i++)out->values[i]=in->values[i+1]-in->values[i];
 if(N-1<out->length)out->values[N-1]=0.0;}
double backward_difference(const DiscreteSignal*s,int n){
 if(!s||!s->values)return 0.0;int i=n-s->start_index;
 if(i-1<0||i>=s->length)return 0.0;return s->values[i]-s->values[i-1];}
void backward_difference_signal(const DiscreteSignal*in,DiscreteSignal*out){
 if(!in||!out||!in->values||!out->values)return;int N=in->length;
 if(out->length>0)out->values[0]=0.0;
 for(int i=1;i<N&&i<out->length;i++)out->values[i]=in->values[i]-in->values[i-1];}
double shift_operator(const DiscreteSignal*s,int n,int k){
 if(!s||!s->values)return 0.0;int i=n+k-s->start_index;
 if(i<0||i>=s->length)return 0.0;return s->values[i];}
void shift_operator_signal(const DiscreteSignal*in,DiscreteSignal*out,int k){
 if(!in||!out||!in->values||!out->values)return;int N=in->length;
 for(int i=0;i<N&&i<out->length;i++){int si=i+k;out->values[i]=(si>=0&&si<N)?in->values[si]:0.0;}}
DiffEqType classify_difference_eq(int order,const double*a,const double*b){
 if(order<1)return DEQ_TYPE_HOMOGENEOUS;int h=0;
 if(b)for(int i=0;i<order;i++)if(fabs(b[i])>1e-15){h=1;break;}
 return h?DEQ_TYPE_NONHOMOGENEOUS:DEQ_TYPE_HOMOGENEOUS;}
void build_characteristic_polynomial(const LinearDiffEq*eq,CharPoly*cp){
 if(!eq||!cp)return;cp->num_roots=eq->order;
 if(cp->num_roots>0){cp->roots=(double complex*)calloc((size_t)cp->num_roots,sizeof(double complex));
 cp->multiplicity=(int*)calloc((size_t)cp->num_roots,sizeof(int));}}
/* L5: Iterative solution */
int solve_iterative_lindiffeq(const LinearDiffEq*eq,const double*x,const double*iy,
 int N,double*y){if(!eq||!y)return -1;int k=eq->order,m=eq->input_order;
 for(int i=0;i<k&&i<N;i++)y[i]=(iy&&i<k)?iy[i]:0.0;
 for(int n=k;n<N;n++){double s=0.0;
  for(int i=1;i<=k&&i<=n;i++)if(eq->a_coeffs)s-=eq->a_coeffs[i-1]*y[n-i];
  if(eq->b_coeffs&&x)for(int j=0;j<=m&&j<=n;j++)s+=eq->b_coeffs[j]*x[n-j];y[n]=s;}return 0;}
/* L5: Characteristic Root Method - Horner + Newton + Deflation */
static void horner_c(const double*cf,int d,double zr,double zi,double*Pr,double*Pi,double*dPr,double*dPi){
 *Pr=cf[0];*Pi=0.0;*dPr=0.0;*dPi=0.0;
 for(int i=1;i<=d;i++){double tr=*dPr*zr-*dPi*zi+*Pr,ti=*dPr*zi+*dPi*zr+*Pi;*dPr=tr;*dPi=ti;
  double t2r=*Pr*zr-*Pi*zi+cf[i],t2i=*Pr*zi+*Pi*zr;*Pr=t2r;*Pi=t2i;}}
static int newton_root(const double*cf,int d,double*zr,double*zi){
 for(int it=0;it<200;it++){double Pr,Pi,dPr,dPi;horner_c(cf,d,*zr,*zi,&Pr,&Pi,&dPr,&dPi);
  double dn=dPr*dPr+dPi*dPi;if(dn<1e-15)break;
  double dr=-(Pr*dPr+Pi*dPi)/dn,di=-(Pi*dPr-Pr*dPi)/dn;*zr+=dr;*zi+=di;if(dr*dr+di*di<1e-14)break;}return 0;}
static void deflate_poly(double*cf,int d,double rr,double ri){
 double*nc=(double*)calloc((size_t)d,sizeof(double));nc[0]=cf[0];
 for(int i=1;i<d;i++)nc[i]=cf[i]+nc[i-1]*rr;memcpy(cf,nc,(size_t)d*sizeof(double));free(nc);}
int solve_characteristic_roots(const LinearDiffEq*eq,const double*iy,int N,double*y){
 if(!eq||!y||eq->order<1)return -1;int k=eq->order;
 double*rr=(double*)calloc((size_t)k,sizeof(double));
 double*ri=(double*)calloc((size_t)k,sizeof(double));
 if(k==1){rr[0]=eq->a_coeffs?-eq->a_coeffs[0]:0.0;ri[0]=0.0;}
 else if(k==2){double a1=eq->a_coeffs?eq->a_coeffs[0]:0.0,a2=eq->a_coeffs?eq->a_coeffs[1]:0.0,
  disc=a1*a1-4.0*a2;if(disc>=0.0){double sd=sqrt(disc);rr[0]=(-a1+sd)/2.0;rr[1]=(-a1-sd)/2.0;ri[0]=ri[1]=0.0;}
  else{rr[0]=rr[1]=-a1/2.0;ri[0]=sqrt(-disc)/2.0;ri[1]=-ri[0];}}
 else{double*pl=(double*)calloc((size_t)(k+1),sizeof(double));pl[0]=1.0;
  if(eq->a_coeffs)for(int i=0;i<k;i++)pl[i+1]=eq->a_coeffs[i];int cd=k;
  for(int ri0=0;ri0<k;ri0++){double zr=0.4+0.3*ri0,zi=0.2;newton_root(pl,cd,&zr,&zi);rr[ri0]=zr;ri[ri0]=zi;
   if(ri0<k-1){deflate_poly(pl,cd,zr,0.0);cd--;}}free(pl);}
 double**V=(double**)calloc((size_t)k,sizeof(double*));
 for(int i=0;i<k;i++)V[i]=(double*)calloc((size_t)k,sizeof(double));
 for(int i=0;i<k;i++)for(int j=0;j<k;j++){double m=sqrt(rr[j]*rr[j]+ri[j]*ri[j]),a=atan2(ri[j],rr[j]);
  V[i][j]=pow(m,i)*cos(i*a);}
 double*rhs=(double*)calloc((size_t)k,sizeof(double));
 if(iy)for(int i=0;i<k;i++)rhs[i]=iy[i];
 for(int col=0;col<k;col++){int pv=col;double mv=fabs(V[col][col]);
  for(int row=col+1;row<k;row++)if(fabs(V[row][col])>mv){mv=fabs(V[row][col]);pv=row;}
  if(pv!=col){for(int j=0;j<k;j++){double t=V[col][j];V[col][j]=V[pv][j];V[pv][j]=t;}
   double t=rhs[col];rhs[col]=rhs[pv];rhs[pv]=t;}
  if(fabs(V[col][col])>1e-12)for(int row=col+1;row<k;row++){double f=V[row][col]/V[col][col];
   for(int j=col;j<k;j++)V[row][j]-=f*V[col][j];rhs[row]-=f*rhs[col];}}
 double*c=(double*)calloc((size_t)k,sizeof(double));
 for(int i=k-1;i>=0;i--){double s=rhs[i];for(int j=i+1;j<k;j++)s-=V[i][j]*c[j];
  c[i]=(fabs(V[i][i])>1e-12)?s/V[i][i]:0.0;}
 for(int n=0;n<N;n++){double yn=0.0;for(int j=0;j<k;j++){
   double m=sqrt(rr[j]*rr[j]+ri[j]*ri[j]),a=atan2(ri[j],rr[j]);yn+=c[j]*pow(m,n)*cos(n*a);}y[n]=yn;}
 for(int i=0;i<k;i++)free(V[i]);free(V);free(c);free(rhs);free(rr);free(ri);return 0;}
/* L5: Undetermined Coefficients */
int solve_undetermined_coeffs(const LinearDiffEq*eq,int rt,double p1,double p2,
 const double*iy,int N,double*y){
 if(!eq||!y)return -1;int k=eq->order;solve_characteristic_roots(eq,iy,N,y);
 switch(rt){
 case 0:{double dn=1.0;if(eq->a_coeffs)for(int i=0;i<k;i++)dn+=eq->a_coeffs[i];
  if(fabs(dn)>1e-12){double Ap=p1/dn;for(int n=0;n<N;n++)y[n]+=Ap;}
  else{double K=p1/(double)k;for(int n=0;n<N;n++)y[n]+=K*n;}break;}
 case 1:{double r=p1,C=p2,Pr=pow(r,k);
  if(eq->a_coeffs)for(int i=0;i<k;i++)Pr+=eq->a_coeffs[i]*pow(r,k-1-i);
  if(fabs(Pr)>1e-12){double Ap=C/Pr;for(int n=0;n<N;n++)y[n]+=Ap*pow(r,n);}break;}
 case 3:{double w=p1,C=p2,Pc=cos(w*k),Ps=sin(w*k);
  if(eq->a_coeffs)for(int i=0;i<k;i++){Pc+=eq->a_coeffs[i]*cos(w*(k-1-i));Ps+=eq->a_coeffs[i]*sin(w*(k-1-i));}
  double det=Pc*Pc+Ps*Ps;if(fabs(det)>1e-12){double Ap=C*Pc/det,Bp=C*Ps/det;
   for(int n=0;n<N;n++)y[n]+=Ap*cos(w*n)+Bp*sin(w*n);}break;}
 }return 0;}
/* L2: Impulse & Step Response */
int compute_impulse_response(const LinearDiffEq*eq,int N,ImpulseResponse*ir){
 if(!eq||!ir||!ir->h)return -1;double*d=(double*)calloc((size_t)N,sizeof(double));
 if(N>0)d[0]=1.0;double*iz=(double*)calloc((size_t)eq->order,sizeof(double));
 int ret=solve_iterative_lindiffeq(eq,d,iz,N,ir->h);ir->is_iir=1;double ts=0.0;int ts0=N/2;
 for(int i=ts0;i<N;i++)ts+=fabs(ir->h[i]);if(ts<1e-10)ir->is_iir=0;free(d);free(iz);return ret;}
int compute_step_response(const ImpulseResponse*ir,StepResponse*sr){
 if(!ir||!sr||!ir->h||!sr->s)return -1;int N=ir->length;if(N>sr->length)N=sr->length;
 double cs=0.0;for(int i=0;i<N;i++){cs+=ir->h[i];sr->s[i]=cs;}sr->steady_state=cs;return 0;}
int compute_step_response_direct(const LinearDiffEq*eq,int N,StepResponse*sr){
 if(!eq||!sr)return -1;double*si=(double*)calloc((size_t)N,sizeof(double));
 for(int i=0;i<N;i++)si[i]=1.0;double*iz=(double*)calloc((size_t)eq->order,sizeof(double));
 int ret=solve_iterative_lindiffeq(eq,si,iz,N,sr->s);free(si);free(iz);return ret;}
/* L4: Superposition & Linearity */
int verify_superposition(const LinearDiffEq*eq,const DiscreteSignal*x1,const DiscreteSignal*x2,
 double alpha,double beta,int N,double tol){
 if(!eq||!x1||!x2)return 0;double*y1=(double*)calloc((size_t)N,sizeof(double));
 double*y2=(double*)calloc((size_t)N,sizeof(double)),*yd=(double*)calloc((size_t)N,sizeof(double));
 double*ys=(double*)calloc((size_t)N,sizeof(double)),*xc=(double*)calloc((size_t)N,sizeof(double));
 double*iz=(double*)calloc((size_t)eq->order,sizeof(double));
 solve_iterative_lindiffeq(eq,x1->values,iz,N,y1);
 solve_iterative_lindiffeq(eq,x2->values,iz,N,y2);
 int N1=x1->length,N2=x2->length;
 for(int i=0;i<N;i++){double v1=(i<N1)?x1->values[i]:0.0,v2=(i<N2)?x2->values[i]:0.0;
  xc[i]=alpha*v1+beta*v2;}solve_iterative_lindiffeq(eq,xc,iz,N,yd);
 for(int i=0;i<N;i++)ys[i]=alpha*y1[i]+beta*y2[i];int ok=1;
 for(int i=0;i<N;i++)if(fabs(yd[i]-ys[i])>tol){ok=0;break;}
 free(y1);free(y2);free(yd);free(ys);free(xc);free(iz);return ok;}
int verify_linearity(const LinearDiffEq*eq,int N,double tol){
 if(!eq)return 0;double*x1=(double*)calloc((size_t)N,sizeof(double));
 double*x2=(double*)calloc((size_t)N,sizeof(double));
 for(int i=0;i<N;i++){x1[i]=sin(0.3*i);x2[i]=cos(0.7*i);}
 double*y1=(double*)calloc((size_t)N,sizeof(double));
 double*y2=(double*)calloc((size_t)N,sizeof(double));
 double*ysm=(double*)calloc((size_t)N,sizeof(double));
 double*xsm=(double*)calloc((size_t)N,sizeof(double));
 double*ysc=(double*)calloc((size_t)N,sizeof(double));
 double*xsc=(double*)calloc((size_t)N,sizeof(double));
 double*iz=(double*)calloc((size_t)eq->order,sizeof(double));
 solve_iterative_lindiffeq(eq,x1,iz,N,y1);solve_iterative_lindiffeq(eq,x2,iz,N,y2);
 for(int i=0;i<N;i++)xsm[i]=x1[i]+x2[i];solve_iterative_lindiffeq(eq,xsm,iz,N,ysm);
 int ao=1;for(int i=0;i<N;i++)if(fabs(ysm[i]-(y1[i]+y2[i]))>tol){ao=0;break;}
 for(int i=0;i<N;i++)xsc[i]=2.5*x1[i];solve_iterative_lindiffeq(eq,xsc,iz,N,ysc);
 int ho=1;for(int i=0;i<N;i++)if(fabs(ysc[i]-2.5*y1[i])>tol){ho=0;break;}
 free(x1);free(x2);free(y1);free(y2);free(ysm);free(xsm);free(ysc);free(xsc);free(iz);return ao&&ho;}
/* L4: Existence & Uniqueness */
ExistenceResult check_existence_uniqueness(const LinearDiffEq*eq,
 const double*ic,int nc){ExistenceResult r={0};
 if(!eq){snprintf(r.message,256,"Null eq");return r;}r.order=eq->order;
 if(eq->order<1){snprintf(r.message,256,"Degenerate order=%d",eq->order);return r;}
 if(nc<eq->order){r.exists=1;r.is_unique=0;
  snprintf(r.message,256,"Insufficient IC: %d of %d",nc,eq->order);return r;}
 if(eq->a_coeffs)for(int i=0;i<eq->order;i++)if(!isfinite(eq->a_coeffs[i])){
  snprintf(r.message,256,"Non-finite a[%d]",i);return r;}
 r.exists=1;r.is_unique=1;r.num_conditions=nc;
 snprintf(r.message,256,"Unique solution exists (order %d)",eq->order);return r;}
/* L1: Fixed Points via Bisection */
int find_fixed_points(double(*f)(double),double xmin,double xmax,double tol,
 FixedPoint*fps,int mp){if(!f||!fps||mp<1)return 0;
 int nf=0,ns=1000;double dx=(xmax-xmin)/ns,h=1e-6;
 for(int i=0;i<ns&&nf<mp;i++){double a=xmin+i*dx,b=a+dx,ga=f(a)-a,gb=f(b)-b;
  if(ga*gb>0.0)continue;if(fabs(ga)<tol){fps[nf].value=a;
   fps[nf].multiplier=fabs((f(a+h)-f(a-h))/(2.0*h));
   fps[nf].type=classify_fixed_point(fps[nf].multiplier);nf++;continue;}
  double lo=a,hi=b;for(int it=0;it<60;it++){double mid=(lo+hi)/2.0,gm=f(mid)-mid;
   if(fabs(gm)<tol){fps[nf].value=mid;fps[nf].multiplier=fabs((f(mid+h)-f(mid-h))/(2.0*h));
    fps[nf].type=classify_fixed_point(fps[nf].multiplier);nf++;break;}
   if(ga*gm<0.0){hi=mid;gb=gm;}else{lo=mid;ga=gm;}}}return nf;}
FixedPointType classify_fixed_point(double m){
 if(m<1.0-1e-10)return FIXED_STABLE;if(m>1.0+1e-10)return FIXED_UNSTABLE;return FIXED_NEUTRAL;}
/* L3: Engineering Response Metrics */
void compute_response_metrics(const StepResponse*sr,double Ts,DiscreteResponseMetrics*m){
 if(!sr||!m)return;memset(m,0,sizeof(DiscreteResponseMetrics));int N=sr->length;if(N<2)return;
 double ss=sr->s[N-1];m->steady_state_error=(fabs(ss)>1e-10)?fabs(1.0-ss):0.0;
 double ssr=(fabs(ss)<1e-10)?1.0:ss;
 for(int i=0;i<N;i++)if(fabs(sr->s[i])>=0.632*fabs(ssr)){m->time_constant=i*Ts;break;}
 int r10=-1,r90=-1;for(int i=0;i<N;i++){
  if(r10<0&&fabs(sr->s[i])>=0.10*fabs(ssr))r10=i;
  if(r90<0&&fabs(sr->s[i])>=0.90*fabs(ssr)){r90=i;break;}}
 if(r10>=0&&r90>=0)m->rise_time=(r90-r10)*Ts;
 double peak=sr->s[0];int pi=0;for(int i=1;i<N;i++)if(sr->s[i]>peak){peak=sr->s[i];pi=i;}
 m->peak_time=pi*Ts;if(fabs(ss)>1e-10){double os=100.0*(peak-fabs(ss))/fabs(ss);m->overshoot_percent=os>0.0?os:0.0;}
 double hi=ssr*1.02,lo=ssr*0.98;int settle=N-1;
 for(int i=N-1;i>=0;i--){if(sr->s[i]>hi||sr->s[i]<lo){settle=((i+1)<N)?i+1:N-1;break;}if(i==0)settle=0;}
 m->settling_time=settle*Ts;
 if(m->overshoot_percent>0.01){double los=log(m->overshoot_percent/100.0);
  m->damping_ratio=-los/sqrt(M_PI*M_PI+los*los);}else m->damping_ratio=1.0;
 if(m->peak_time>1e-10&&m->damping_ratio<0.999)
  m->natural_frequency=M_PI/(m->peak_time*sqrt(1.0-m->damping_ratio*m->damping_ratio));}
/* L3: Discrete-domain metrics from poles */
double time_constant_discrete(const CharPoly*cp){
 if(!cp||cp->num_roots<1)return 0.0;double mm=0.0;
 for(int i=0;i<cp->num_roots;i++){double m=cabs(cp->roots[i]);if(m<1.0&&m>mm)mm=m;}
 if(mm<1e-15||mm>=1.0)return INFINITY;return -1.0/log(mm);}
double settling_time_discrete(const CharPoly*cp){double t=time_constant_discrete(cp);return isinf(t)?INFINITY:4.0*t;}
double damping_ratio_discrete(const CharPoly*cp){
 if(!cp||cp->num_roots<1)return 0.0;double bz=1.0;
 for(int i=0;i<cp->num_roots;i++){double r=cabs(cp->roots[i]),th=fabs(carg(cp->roots[i]));
  if(r>=1.0||r<1e-10||th<1e-10)continue;double lr=-log(r),z=lr/sqrt(lr*lr+th*th);if(z<bz)bz=z;}
 if(bz<0.0)bz=0.0;if(bz>1.0)bz=1.0;return bz;}
double natural_frequency_discrete(const CharPoly*cp){
 if(!cp||cp->num_roots<1)return 0.0;double bw=0.0;
 for(int i=0;i<cp->num_roots;i++){double r=cabs(cp->roots[i]),th=fabs(carg(cp->roots[i]));
  if(r>=1.0||th<1e-10)continue;double lr=log(r),wn=sqrt(lr*lr+th*th);if(wn>bw)bw=wn;}return bw;}
double overshoot_percent_discrete(const CharPoly*cp){double z=damping_ratio_discrete(cp);
 if(z>=1.0)return 0.0;return 100.0*exp(-M_PI*z/sqrt(1.0-z*z));}
/* L1: Nonlinear Map Iteration */
double nonlinear_map_eval(double(*f)(double),double yn){return f?f(yn):yn;}
int nonlinear_map_iterate(double(*f)(double),double y0,int N,double*traj){
 if(!f||!traj||N<1)return -1;traj[0]=y0;for(int n=0;n<N-1;n++)traj[n+1]=f(traj[n]);return 0;}
/* Utility: memory management */
void linear_diffeq_free(LinearDiffEq*eq){if(!eq)return;free(eq->a_coeffs);eq->a_coeffs=NULL;free(eq->b_coeffs);eq->b_coeffs=NULL;}
void discrete_signal_free(DiscreteSignal*s){if(!s)return;free(s->values);s->values=NULL;}
void char_poly_free(CharPoly*cp){if(!cp)return;free(cp->roots);cp->roots=NULL;free(cp->multiplicity);cp->multiplicity=NULL;}
void impulse_response_free(ImpulseResponse*ir){if(!ir)return;free(ir->h);ir->h=NULL;}
void step_response_free(StepResponse*sr){if(!sr)return;free(sr->s);sr->s=NULL;}
DiscreteSignal* discrete_signal_alloc(int length,int start_index){
 DiscreteSignal*s=(DiscreteSignal*)calloc(1,sizeof(DiscreteSignal));if(!s)return NULL;
 s->length=length;s->start_index=start_index;
 if(length>0){s->values=(double*)calloc((size_t)length,sizeof(double));if(!s->values){free(s);return NULL;}}return s;}
LinearDiffEq* linear_diffeq_alloc(int order,int input_order){
 LinearDiffEq*e=(LinearDiffEq*)calloc(1,sizeof(LinearDiffEq));if(!e)return NULL;
 e->order=order;e->input_order=input_order;
 if(order>0){e->a_coeffs=(double*)calloc((size_t)order,sizeof(double));if(!e->a_coeffs){free(e);return NULL;}}
 if(input_order>=0){e->b_coeffs=(double*)calloc((size_t)(input_order+1),sizeof(double));
  if(!e->b_coeffs){linear_diffeq_free(e);free(e);return NULL;}}e->type=DEQ_TYPE_LINEAR_CONSTANT;return e;}
int linear_diffeq_set_coeffs(LinearDiffEq*eq,const double*a,const double*b){
 if(!eq)return -1;if(a&&eq->a_coeffs)memcpy(eq->a_coeffs,a,(size_t)eq->order*sizeof(double));
 if(b&&eq->b_coeffs)memcpy(eq->b_coeffs,b,(size_t)(eq->input_order+1)*sizeof(double));return 0;}
void linear_diffeq_print(const LinearDiffEq*eq){
 if(!eq){printf("null\n");return;}printf("y[n]");
 if(eq->a_coeffs)for(int i=0;i<eq->order;i++)if(fabs(eq->a_coeffs[i])>1e-14)printf(" %+.4fy[n-%d]",eq->a_coeffs[i],i+1);
 printf(" =");
 if(eq->b_coeffs)for(int i=0;i<=eq->input_order;i++)if(fabs(eq->b_coeffs[i])>1e-14||i==0)printf(" %+.4fx[n-%d]",eq->b_coeffs[i],i);
 printf("\n");}
