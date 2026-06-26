/* applications.c - Real-World Applications of Difference Equations
 * L6: Digital filter design, Fibonacci, logistic map, ARMA, PID
 * L7: Kalman prediction, DC motor, exponential smoothing, resonator
 * Reference: Ogata "Discrete-Time Control Systems" (1995)
 * Reference: Oppenheim & Schafer "Discrete-Time Signal Processing" (2010)
 * Reference: Box & Jenkins "Time Series Analysis" (2015)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include "diffeq.h"
#include "signal.h"
#include "ztransform.h"
#include "stability.h"
#include "statespace.h"
/* L6: Fibonacci closed form (Binet formula) via characteristic equation */
void fibonacci_closed_form(int n,double*Fn){
 /* F[n] satisfies: F[n+2]=F[n+1]+F[n], F[0]=0,F[1]=1
  * Char eq: λ²-λ-1=0 → λ=(1±√5)/2 */
 double phi=(1.0+sqrt(5.0))/2.0,psi=(1.0-sqrt(5.0))/2.0;
 *Fn=(pow(phi,n)-pow(psi,n))/sqrt(5.0);}
/* L6: Fibonacci iterative (verification) */
int fibonacci_iterative(int n,double*values,int count){
 if(!values||count<2)return -1;values[0]=0.0;if(count>1)values[1]=1.0;
 for(int i=2;i<count;i++)values[i]=values[i-1]+values[i-2];
 double fn;fibonacci_closed_form(n,&fn);return 0;}
/* L6: Logistic map iteration (population dynamics) */
int logistic_map_iterate(double r,double x0,int N,double*trajectory){
 /* x[n+1] = r * x[n] * (1 - x[n]), x in [0,1], r in [0,4]
  * Classic chaotic system at r>3.57 */
 if(!trajectory||N<1)return -1;trajectory[0]=x0;
 for(int n=0;n<N-1;n++){double xn=trajectory[n];
  trajectory[n+1]=r*xn*(1.0-xn);}return 0;
}
/* L6: Logistic map bifurcation diagram helper */
void logistic_bifurcation_data(double r_min,double r_max,int r_steps,int transient,int steady,double*results){
 double dr=(r_max-r_min)/r_steps;
 for(int i=0;i<=r_steps;i++){double r=r_min+i*dr;double x=0.5;
  for(int t=0;t<transient;t++)x=r*x*(1.0-x);
  for(int s=0;s<steady;s++){x=r*x*(1.0-x);results[i*steady+s]=x;}}
}
/* L6: Design IIR Butterworth lowpass via bilinear transform */
int design_iir_lowpass(double cutoff_freq,double sampling_freq,int order,IIRFilter*filter){
 if(!filter||order<1||order>8)return -1;filter->type=FILTER_LOWPASS;filter->order=order;
 filter->f_cutoff=cutoff_freq;filter->sampling_freq=sampling_freq;
 double wc=2.0*sampling_freq*tan(M_PI*cutoff_freq/sampling_freq);/* Pre-warp */
 /* Butterworth prototype poles in s-plane */
 /* For order 2: s^2 + sqrt(2)*wc*s + wc^2 */
 if(order==2){double a0=wc*wc,a1=sqrt(2.0)*wc;
  /* Bilinear: H(z)=H(s)|_{s=2/T*(z-1)/(z+1)} */
  double T=1.0/sampling_freq,k=2.0/T;
  double b0=a0,b1=0.0,b2=-a0;double den=k*k+a1*k+a0;
  filter->b_coeffs=(double*)calloc(3,sizeof(double));
  filter->a_coeffs=(double*)calloc(2,sizeof(double));
  filter->b_coeffs[0]=b0/den;filter->b_coeffs[1]=2.0*b0/den;filter->b_coeffs[2]=b0/den;
  filter->a_coeffs[0]=(2.0*k*k-2.0*a0)/den;filter->a_coeffs[1]=(k*k-a1*k+a0)/den;}
 return 0;}
/* L6: Design FIR filter via window method */
int design_fir_filter(FilterType type,int order,double fc_low,double fc_high,WindowType win,FIRFilter*fir){
 if(!fir||order<1)return -1;fir->type=type;fir->order=order;fir->window=win;
 fir->coefficients=(double*)calloc((size_t)(order+1),sizeof(double));
 int M=order;double*ideal_h=(double*)calloc((size_t)(M+1),sizeof(double));
 int mid=M/2;double wc_low=2.0*M_PI*fc_low,wc_high=2.0*M_PI*fc_high;
 /* Ideal impulse response */
 for(int n=0;n<=M;n++){int k=n-mid;
  if(type==FILTER_LOWPASS){ideal_h[n]=(k==0)?wc_low/M_PI:sin(wc_low*k)/(M_PI*k);}
  else if(type==FILTER_HIGHPASS){ideal_h[n]=(k==0)?1.0-wc_low/M_PI:-sin(wc_low*k)/(M_PI*k);}
  else if(type==FILTER_BANDPASS){ideal_h[n]=(k==0)?(wc_high-wc_low)/M_PI:
   (sin(wc_high*k)-sin(wc_low*k))/(M_PI*k);}
  else{ideal_h[n]=(k==0)?1.0-(wc_high-wc_low)/M_PI:
   (sin(wc_low*k)-sin(wc_high*k))/(M_PI*k);}}
 /* Window function */
 double*w=(double*)calloc((size_t)(M+1),sizeof(double));
 window_generate(w,M+1,win,5.0);
 for(int n=0;n<=M;n++)fir->coefficients[n]=ideal_h[n]*w[n];
 free(ideal_h);free(w);return 0;}
/* L6: ARMA(p,q) process simulation */
int arma_model_simulate(const double*ar_coeffs,int p,const double*ma_coeffs,int q,
 double noise_std,int N,double*y){
 if(!y||N<1)return -1;double*noise=(double*)calloc((size_t)N,sizeof(double));
 signal_white_noise(noise,N,0.0,noise_std);/* AR part: y[n]+Σa_i y[n-i]=noise[n]+Σb_j noise[n-j] */
 for(int n=0;n<N;n++){y[n]=noise[n];
  if(ar_coeffs)for(int i=1;i<=p&&i<=n;i++)y[n]-=ar_coeffs[i-1]*y[n-i];
  if(ma_coeffs)for(int j=0;j<=q&&j<=n;j++)y[n]+=ma_coeffs[j]*noise[n-j];}
 free(noise);return 0;}
/* L6: Economic cobweb model */
int cobweb_model(double a,double b,double initial_price,int N,double*prices,double*quantities){
 /* Demand: Qd[n]=a-b*P[n], Supply: Qs[n]=c+d*P[n-1] (c=0,d=1 for simplicity)
  * Equilibrium: P[n]=(a-Q[n-1])/b, Q[n]=P[n-1] */
 if(!prices||!quantities||N<1)return -1;prices[0]=initial_price;quantities[0]=prices[0];
 for(int n=1;n<N;n++){quantities[n]=prices[n-1];prices[n]=(a-quantities[n])/b;}return 0;}
/* L6: Moving average filter */
int moving_average_filter(const double*x,int N,int window,double*y){
 if(!x||!y||window<1)return -1;for(int n=0;n<N;n++){double s=0.0;int cnt=0;
  for(int k=0;k<window&&n-k>=0;k++){s+=x[n-k];cnt++;}y[n]=cnt>0?s/cnt:0.0;}return 0;}
/* L7: Discrete PID controller (Tustin discretization with anti-windup) */
void discrete_pid_init(DiscretePID*pid,double Kp,double Ki,double Kd,double Ts,double out_min,double out_max){
 pid->Kp=Kp;pid->Ki=Ki;pid->Kd=Kd;pid->Ts=Ts;pid->integral=0.0;pid->prev_error=0.0;
 pid->prev_output=0.0;pid->out_min=out_min;pid->out_max=out_max;pid->has_aw=1;}
double discrete_pid_update(DiscretePID*pid,double setpoint,double measurement){
 double error=setpoint-measurement;pid->integral+=pid->Ki*pid->Ts*error;
 /* Anti-windup: clamp integral */
 if(pid->has_aw){if(pid->integral>pid->out_max)pid->integral=pid->out_max;
  if(pid->integral<pid->out_min)pid->integral=pid->out_min;}
 double derivative=pid->Kd*(error-pid->prev_error)/pid->Ts;
 double output=pid->Kp*error+pid->integral+derivative;
 if(output>pid->out_max)output=pid->out_max;if(output<pid->out_min)output=pid->out_min;
 pid->prev_error=error;pid->prev_output=output;return output;}
/* L7: Kalman filter prediction step */
typedef struct {int n,m;double**A;double**B;double**Q;double*x_hat;double**P;} KalmanPredict;
void kalman_predict_init(KalmanPredict*kf,int n,int m){
 kf->n=n;kf->m=m;kf->x_hat=(double*)calloc((size_t)n,sizeof(double));
 kf->P=(double**)calloc((size_t)n,sizeof(double*));
 for(int i=0;i<n;i++)kf->P[i]=(double*)calloc((size_t)n,sizeof(double));}
int kalman_predict_step(KalmanPredict*kf,const double*u,double*x_pred,double**P_pred){
 if(!kf||!x_pred)return -1;int n=kf->n;
 /* x_pred = A*x_hat + B*u */
 for(int i=0;i<n;i++){x_pred[i]=0.0;
  for(int j=0;j<n;j++)x_pred[i]+=kf->A[i][j]*kf->x_hat[j];
  if(kf->B&&u)for(int j=0;j<kf->m;j++)x_pred[i]+=kf->B[i][j]*u[j];}
 /* P_pred = A*P*A^T + Q */
 for(int i=0;i<n;i++)for(int j=0;j<n;j++){P_pred[i][j]=kf->Q[i][j];
  for(int l=0;l<n;l++)for(int r=0;r<n;r++)P_pred[i][j]+=kf->A[i][l]*kf->P[l][r]*kf->A[j][r];}
 return 0;}
void kalman_predict_free(KalmanPredict*kf){if(!kf)return;free(kf->x_hat);
 for(int i=0;i<kf->n;i++)free(kf->P[i]);free(kf->P);}
/* L7: Exponential smoothing (1st-order IIR lowpass) */
double exponential_smoothing(double alpha,double x,double*prev_y){
 double y=alpha*x+(1.0-alpha)*(*prev_y);*prev_y=y;return y;}
/* L7: Digital resonator (2nd-order oscillator) */
double digital_resonator(double r,double omega0,double x,double*y1,double*y2){
 /* y[n]=2r*cos(ω0)*y[n-1]-r²*y[n-2]+x[n] */
 double y=2.0*r*cos(omega0)*(*y1)-r*r*(*y2)+x;*y2=*y1;*y1=y;return y;}
/* L7: DC motor discrete model (1st order: J dw/dt + B w = Kt i) */
void dc_motor_discrete_init(DCMotorDisc*motor,double J,double B,double Kt,double Ts){
 motor->J=J;motor->B=B;motor->Kt=Kt;motor->Ts=Ts;motor->prev_w=0.0;motor->prev_i=0.0;}
double dc_motor_discrete_update(DCMotorDisc*motor,double voltage,double R){
 /* w[n] = (Ts*Kt*i[n] + J*w[n-1])/(J+Ts*B) */
 double i=voltage/R;double w=(motor->Ts*motor->Kt*i+motor->J*motor->prev_w)/(motor->J+motor->Ts*motor->B);
 motor->prev_w=w;motor->prev_i=i;return w;}
/* L6: AR(1) autoregressive model estimation */
typedef struct {double phi,sigma2;int N;double*y;} AR1Model;
void ar1_model_fit(const double*data,int N,AR1Model*model){
 if(!data||!model||N<2)return;model->N=N;double sum_yy=0.0,sum_yy1=0.0;
 for(int i=1;i<N;i++){sum_yy+=data[i]*data[i];sum_yy1+=data[i]*data[i-1];}
 double sum_y1y1=0.0;for(int i=0;i<N-1;i++)sum_y1y1+=data[i]*data[i];
 model->phi=(fabs(sum_y1y1)>1e-12)?sum_yy1/sum_y1y1:0.0;
 double ss_res=0.0;for(int i=1;i<N;i++){double e=data[i]-model->phi*data[i-1];ss_res+=e*e;}
 model->sigma2=ss_res/(N-1);}
/* L6: AR(2) model via Yule-Walker equations */
typedef struct {double phi1,phi2,sigma2;} AR2Model;
void ar2_yule_walker(const double*data,int N,AR2Model*model){
 if(!data||!model||N<3)return;double r0=0.0,r1=0.0,r2=0.0;
 for(int i=0;i<N;i++)r0+=data[i]*data[i];r0/=N;
 for(int i=0;i<N-1;i++)r1+=data[i]*data[i+1];r1/=N;
 for(int i=0;i<N-2;i++)r2+=data[i]*data[i+2];r2/=N;
 /* [r0 r1; r1 r0] [phi1; phi2] = [r1; r2] */
 double det=r0*r0-r1*r1;if(fabs(det)<1e-12){model->phi1=model->phi2=0.0;model->sigma2=r0;return;}
 double r1_r0=r1/r0,r2_r0=r2/r0;model->phi1=(r1_r0*(1.0-r2_r0/r1_r0))/det;
 model->phi2=(r2_r0-r1_r0*model->phi1);
 double ss=0.0;for(int i=2;i<N;i++){double e=data[i]-model->phi1*data[i-1]-model->phi2*data[i-2];ss+=e*e;}
 model->sigma2=ss/(N-2);}
/* L7: Audio echo effect via difference equation (comb filter) */
typedef struct {int delay;double decay;double*buffer;int buf_idx;int buf_size;} AudioEcho;
void audio_echo_init(AudioEcho*echo,int delay,double decay){
 echo->delay=delay;echo->decay=decay;echo->buf_size=delay+1;
 echo->buffer=(double*)calloc((size_t)echo->buf_size,sizeof(double));echo->buf_idx=0;}
double audio_echo_process(AudioEcho*echo,double x){
 /* y[n]=x[n]+decay*y[n-delay] */
 double y_delayed=echo->buffer[(echo->buf_idx-echo->delay+echo->buf_size)%echo->buf_size];
 double y=x+echo->decay*y_delayed;echo->buffer[echo->buf_idx]=y;
 echo->buf_idx=(echo->buf_idx+1)%echo->buf_size;return y;}
void audio_echo_free(AudioEcho*echo){if(!echo)return;free(echo->buffer);}
/* L7: Reverb via Schroeder allpass comb filter */
typedef struct {int*delays;double*gains;int num_stages;double**buffers;int*buf_idxs;int*buf_sizes;} SchroederReverb;
void schroeder_reverb_init(SchroederReverb*rev,const int*delays,const double*gains,int num_stages){
 rev->num_stages=num_stages;rev->delays=(int*)calloc((size_t)num_stages,sizeof(int));
 rev->gains=(double*)calloc((size_t)num_stages,sizeof(double));
 rev->buffers=(double**)calloc((size_t)num_stages,sizeof(double*));
 rev->buf_idxs=(int*)calloc((size_t)num_stages,sizeof(int));
 rev->buf_sizes=(int*)calloc((size_t)num_stages,sizeof(int));
 for(int i=0;i<num_stages;i++){rev->delays[i]=delays[i];rev->gains[i]=gains[i];
  rev->buf_sizes[i]=delays[i]+1;rev->buffers[i]=(double*)calloc((size_t)rev->buf_sizes[i],sizeof(double));
  rev->buf_idxs[i]=0;}}
double schroeder_reverb_process(SchroederReverb*rev,double x){
 double y=x;for(int i=0;i<rev->num_stages;i++){int d=rev->delays[i];double g=rev->gains[i];
  int idx=(rev->buf_idxs[i]-d+rev->buf_sizes[i])%rev->buf_sizes[i];
  double delayed=rev->buffers[i][idx];double out=y+g*delayed;
  rev->buffers[i][rev->buf_idxs[i]]=out-delayed*g;rev->buf_idxs[i]=(rev->buf_idxs[i]+1)%rev->buf_sizes[i];
  y=out;}return y;}
void schroeder_reverb_free(SchroederReverb*rev){if(!rev)return;
 for(int i=0;i<rev->num_stages;i++)free(rev->buffers[i]);free(rev->buffers);free(rev->delays);free(rev->gains);
 free(rev->buf_idxs);free(rev->buf_sizes);}
/* L7: Tesla/SpaceX relevant - discrete battery SoC estimation via Coulomb counting */
double battery_soc_update(double prev_soc,double current,double dt,double capacity_ah){
 /* SoC[n+1]=SoC[n]-I*dt/(3600*capacity_Ah) Coulomb counting difference equation */
 double delta_soc=-current*dt/(3600.0*capacity_ah);double soc=prev_soc+delta_soc;
 if(soc>1.0)soc=1.0;if(soc<0.0)soc=0.0;return soc;}
/* L7: GPS position estimation via constant velocity model (discrete state-space) */
void gps_cv_model_predict(double*x,double*y,double*vx,double*vy,double dt){
 /* x[n+1]=x[n]+vx*dt, y[n+1]=y[n]+vy*dt, vx[n+1]=vx[n], vy[n+1]=vy[n] */
 *x+=*vx*dt;*y+=*vy*dt;/* velocity stays constant */}
/* L6: Compound interest / mortgage calculator via difference equation */
typedef struct {double principal;double rate;int periods;double payment;} LoanModel;
void loan_amortization_schedule(const LoanModel*loan,int N,double*balance,double*interest_paid,double*principal_paid){
 if(!loan||!balance||N<1)return;double b=loan->principal;double r=loan->rate;double pmt=loan->payment;
 for(int n=0;n<N;n++){double int_paid=b*r;double prin_paid=pmt-int_paid;
  if(prin_paid>b)prin_paid=b;b-=prin_paid;if(b<0.0)b=0.0;
  balance[n]=b;if(interest_paid)interest_paid[n]=int_paid;
  if(principal_paid)principal_paid[n]=prin_paid;}}
/* L6: Gambler's ruin probability via difference equation */
double gamblers_ruin_probability(double p,int initial,int N){
 /* P_k = probability of reaching N before 0, starting from k.
  * Difference eq: P_k = p*P_{k+1} + (1-p)*P_{k-1}, P_0=0, P_N=1.
  * Solution: P_k = (1-(q/p)^k)/(1-(q/p)^N) for p!=0.5, else P_k=k/N */
 if(fabs(p-0.5)<1e-12)return(double)initial/N;
 double q=1.0-p,r=q/p;if(fabs(r-1.0)<1e-12)return(double)initial/N;
 return(1.0-pow(r,initial))/(1.0-pow(r,N));}
/* L6: Malthusian population growth (geometric): P[n+1]=(1+r)*P[n] */
void malthusian_growth(double P0,double r,int N,double*population){
 if(!population||N<1)return;population[0]=P0;
 for(int n=1;n<N;n++)population[n]=population[n-1]*(1.0+r);}
/* L6: Verhulst (logistic) continuous discretized: P[n+1]=P[n]+r*P[n]*(1-P[n]/K) */
void verhulst_logistic_discrete(double P0,double r,double K,int N,double*population){
 if(!population||N<1)return;population[0]=P0;
 for(int n=0;n<N-1;n++){double Pn=population[n];
  population[n+1]=Pn+r*Pn*(1.0-Pn/K);}}
/* L7: Electrical circuit RC low-pass filter discretized (backward Euler) */
void rc_lowpass_init(RCLowPass*rc,double R,double C,double Ts){
 rc->R=R;rc->C=C;rc->Ts=Ts;rc->Vc_prev=0.0;}
double rc_lowpass_update(RCLowPass*rc,double Vin){
 /* Vc[n] = (Ts*Vin + RC*Vc[n-1])/(RC + Ts) where RC=R*C */
 double RC=rc->R*rc->C;double Vc=(rc->Ts*Vin+RC*rc->Vc_prev)/(RC+rc->Ts);
 rc->Vc_prev=Vc;return Vc;}
/* L7: Mechanical mass-spring-damper discrete model */
void msd_discrete_init(MassSpringDamper*msd,double m,double c,double k,double Ts){
 msd->m=m;msd->c=c;msd->k=k;msd->Ts=Ts;msd->x_prev=0.0;msd->v_prev=0.0;}
void msd_discrete_update(MassSpringDamper*msd,double F,double*x,double*v){
 /* m*a+c*v+k*x=F, discretized with backward Euler:
  * v[n]=(m*v[n-1]+Ts*(F-c*v[n-1]-k*x[n-1]))/(m+Ts*c+Ts^2*k)
  * x[n]=x[n-1]+Ts*v[n] */
 double m=msd->m,c=msd->c,k=msd->k,Ts=msd->Ts;
 double den=m+Ts*c+Ts*Ts*k;if(fabs(den)<1e-15)den=1e-15;
 double v_new=(m*msd->v_prev+Ts*(F-c*msd->v_prev-k*msd->x_prev))/den;
 double x_new=msd->x_prev+Ts*v_new;
 *v=v_new;*x=x_new;msd->x_prev=x_new;msd->v_prev=v_new;}
/* L7: Newton's law of cooling discrete model */
void newton_cooling_init(NewtonCooling*nc,double k,double T_amb,double Ts){
 nc->k=k;nc->T_amb=T_amb;nc->Ts=Ts;nc->T_prev=T_amb;}
double newton_cooling_update(NewtonCooling*nc,double T_current){
 /* T[n+1]=T[n]+k*(T_amb-T[n])*Ts */
 double T_next=nc->T_prev+nc->k*(nc->T_amb-nc->T_prev)*nc->Ts;
 nc->T_prev=T_next;return T_next;}
/* L7: Quadcopter altitude discrete model (1D) */
void quad_alt_disc_init(QuadAltDisc*q,double mass,double k_thrust,double Ts){
 q->mass=mass;q->g=9.81;q->k_thrust=k_thrust;q->Ts=Ts;q->z_prev=0.0;q->v_prev=0.0;}
void quad_alt_disc_update(QuadAltDisc*q,double throttle,double*z,double*v){
 /* m*a = k_thrust*throttle - m*g, discretized */
 double F=q->k_thrust*throttle-q->mass*q->g;double a=F/q->mass;
 double v_new=q->v_prev+a*q->Ts;double z_new=q->z_prev+q->v_prev*q->Ts+0.5*a*q->Ts*q->Ts;
 *v=v_new;*z=z_new;q->z_prev=z_new;q->v_prev=v_new;}
/* L7: Mars rover wheel slip detection via difference of encoder counts */
typedef struct {double wheel_radius;int prev_left_ticks,prev_right_ticks;double Ts;} RoverSlip;
typedef struct {double Kp,Ki,Ts;double integral,prev_freq_error;} LFCController;
typedef struct {double C_atm,C_ocean,C_land;double k_ao,k_oa,k_al,k_la;double Ts;} CarbonModel;
void rover_slip_init(RoverSlip*rs,double wheel_radius,double Ts){
 rs->wheel_radius=wheel_radius;rs->prev_left_ticks=0;rs->prev_right_ticks=0;rs->Ts=Ts;}
double rover_slip_detect(RoverSlip*rs,int left_ticks,int right_ticks,int ticks_per_rev){
 /* Slip ratio = (omega*r - v)/omega*r where omega from encoder, v from difference */
 int dl=left_ticks-rs->prev_left_ticks,dr=right_ticks-rs->prev_right_ticks;
 double omega=2.0*M_PI*(dl+dr)/(2.0*ticks_per_rev*rs->Ts);
 double v=rs->wheel_radius*2.0*M_PI*(dl-dr)/(2.0*ticks_per_rev*rs->Ts);
 rs->prev_left_ticks=left_ticks;rs->prev_right_ticks=right_ticks;
 double v_wheel=omega*rs->wheel_radius;if(fabs(v_wheel)<1e-6)return 0.0;
 return 1.0-v/v_wheel;}
/* L7: Smart grid load frequency control (discrete PI) */
void lfc_init(LFCController*lfc,double Kp,double Ki,double Ts){
 lfc->Kp=Kp;lfc->Ki=Ki;lfc->Ts=Ts;lfc->integral=0.0;lfc->prev_freq_error=0.0;}
double lfc_update(LFCController*lfc,double freq_error){
 lfc->integral+=lfc->Ki*lfc->Ts*freq_error;double u=lfc->Kp*freq_error+lfc->integral;
 lfc->prev_freq_error=freq_error;return u;}
/* L7: Climate carbon cycle box model (3-reservoir discrete) */
void carbon_model_init(CarbonModel*cm,double k_ao,double k_oa,double k_al,double k_la,double Ts){
 cm->C_atm=600.0;cm->C_ocean=38000.0;cm->C_land=2300.0;/* PgC initial */
 cm->k_ao=k_ao;cm->k_oa=k_oa;cm->k_al=k_al;cm->k_la=k_la;cm->Ts=Ts;}
void carbon_model_step(CarbonModel*cm,double emissions){
 /* dC_atm/dt = emissions + k_oa*C_ocean - k_ao*C_atm + k_la*C_land - k_al*C_atm
  * Discretized via forward Euler */
 double dC_atm=emissions+cm->k_oa*cm->C_ocean-cm->k_ao*cm->C_atm+cm->k_la*cm->C_land-cm->k_al*cm->C_atm;
 double dC_ocean=cm->k_ao*cm->C_atm-cm->k_oa*cm->C_ocean;
 double dC_land=cm->k_al*cm->C_atm-cm->k_la*cm->C_land;
 cm->C_atm+=dC_atm*cm->Ts;cm->C_ocean+=dC_ocean*cm->Ts;cm->C_land+=dC_land*cm->Ts;}
/* L7: Beer fermentation temperature control (discrete PI with constraints) */
double beer_fermentation_control(double T_target,double T_measured,double*integral,
 double Kp,double Ki,double Ts,double T_min,double T_max){
 double error=T_target-T_measured;*integral+=Ki*Ts*error;
 if(*integral>T_max)*integral=T_max;if(*integral<T_min)*integral=T_min;
 double u=Kp*error+*integral;if(u>T_max)u=T_max;if(u<T_min)u=T_min;return u;}
/* L7: 2020 pandemic SEIR discrete model */
void seir_model_init(SEIRModel*m,double S0,double E0,double I0,double R0,
 double beta,double sigma,double gamma){
 m->S=S0;m->E=E0;m->Inf=I0;m->R=R0;m->beta=beta;m->sigma=sigma;m->gamma=gamma;
 m->N=(int)(S0+E0+I0+R0);}
void seir_model_step(SEIRModel*m){
 double N=(double)m->N;double new_exposed=m->beta*m->S*m->Inf/N;
 double new_infectious=m->sigma*m->E;double new_recovered=m->gamma*m->Inf;
 m->S-=new_exposed;m->E+=new_exposed-new_infectious;m->Inf+=new_infectious-new_recovered;
 m->R+=new_recovered;}
/* L6: Forced harmonic oscillator discrete: my''+cy'+ky=F(t) */
void forced_oscillator_discrete(double m,double c,double k,const double*F,int N,double Ts,
 double*displacement,double*velocity){
 if(!F||!displacement||!velocity||N<1)return;displacement[0]=0.0;velocity[0]=0.0;
 for(int n=0;n<N-1;n++){double x=displacement[n],v=velocity[n],Fn=F[n];
  double a=(Fn-c*v-k*x)/m;velocity[n+1]=v+a*Ts;
  displacement[n+1]=x+v*Ts+0.5*a*Ts*Ts;}}
/* L6: Sampled-data system with ZOH (zero-order hold) */
void sampled_data_zoh_init(SampledDataZOH*sys,double A,double B,double C,double Ts){
 sys->A=A;sys->B=B;sys->C=C;sys->Ts=Ts;sys->x_prev=0.0;}
double sampled_data_zoh_update(SampledDataZOH*sys,double u){
 /* Exact discretization: x[n+1]=e^{AT}*x[n]+(e^{AT}-1)/A*B*u[n] */
 double eAT=exp(sys->A*sys->Ts);double Bd=(fabs(sys->A)>1e-12)?(eAT-1.0)/sys->A*sys->B:sys->B*sys->Ts;
 double x_new=eAT*sys->x_prev+Bd*u;double y=sys->C*x_new;sys->x_prev=x_new;return y;}
/* L6: Discretized predator-prey (Lotka-Volterra) */
void lotka_volterra_discrete(double alpha,double beta,double gamma,double delta,
 double*prey,double*predator,int N){
 /* dx/dt=alpha*x-beta*x*y, dy/dt=delta*x*y-gamma*y. Discrete: simple Euler */
 if(!prey||!predator||N<1)return;for(int n=0;n<N-1;n++){double x=prey[n],y=predator[n];
  prey[n+1]=x+(alpha*x-beta*x*y);predator[n+1]=y+(delta*x*y-gamma*y);
  if(prey[n+1]<0.0)prey[n+1]=0.0;if(predator[n+1]<0.0)predator[n+1]=0.0;}}
