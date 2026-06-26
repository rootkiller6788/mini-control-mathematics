/* signal.c - Signal generation and processing implementations */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "signal.h"
void signal_impulse(double*out,int N,int delay){for(int i=0;i<N;i++)out[i]=(i==delay)?1.0:0.0;}
void signal_step(double*out,int N,int start){for(int i=0;i<N;i++)out[i]=(i>=start)?1.0:0.0;}
void signal_ramp(double*out,int N){for(int i=0;i<N;i++)out[i]=(double)i;}
void signal_exponential(double*out,int N,double base){out[0]=1.0;for(int i=1;i<N;i++)out[i]=out[i-1]*base;}
void signal_complex_exponential(double*re,double*im,int N,double w){for(int i=0;i<N;i++){re[i]=cos(w*i);im[i]=sin(w*i);}}
void signal_sinusoid(double*out,int N,double A,double w,double ph){for(int i=0;i<N;i++)out[i]=A*cos(w*i+ph);}
void signal_white_noise(double*out,int N,double mean,double std){
 static unsigned int seed=12345;for(int i=0;i<N;i++){
  double u1=(double)((seed=seed*1103515245+12345)>>16)/32768.0;
  double u2=(double)((seed=seed*1103515245+12345)>>16)/32768.0;
  double r=sqrt(-2.0*log(u1+1e-15));out[i]=mean+std*r*cos(2.0*M_PI*u2);}}
void signal_chirp(double*out,int N,double A,double w0,double k){for(int i=0;i<N;i++)out[i]=A*cos((w0+0.5*k*i)*i);}
double signal_energy(const double*x,int N){double e=0.0;for(int i=0;i<N;i++)e+=x[i]*x[i];return e;}
double signal_power(const double*x,int N){return signal_energy(x,N)/N;}
double signal_rms(const double*x,int N){return sqrt(signal_power(x,N));}
double signal_peak(const double*x,int N,int*pi){double p=x[0];*pi=0;for(int i=1;i<N;i++)if(x[i]>p){p=x[i];*pi=i;}return p;}
double signal_mean(const double*x,int N){double s=0.0;for(int i=0;i<N;i++)s+=x[i];return s/N;}
double signal_variance(const double*x,int N){double m=signal_mean(x,N),v=0.0;for(int i=0;i<N;i++){double d=x[i]-m;v+=d*d;}return v/N;}
void window_generate(double*w,int N,WindowType t,double kb){for(int i=0;i<N;i++){
 double n=(double)i;switch(t){case WINDOW_RECTANGULAR:w[i]=1.0;break;
 case WINDOW_HANN:w[i]=0.5*(1.0-cos(2.0*M_PI*n/(N-1)));break;
 case WINDOW_HAMMING:w[i]=0.54-0.46*cos(2.0*M_PI*n/(N-1));break;
 case WINDOW_BLACKMAN:w[i]=0.42-0.5*cos(2.0*M_PI*n/(N-1))+0.08*cos(4.0*M_PI*n/(N-1));break;
 case WINDOW_BARTLETT:w[i]=1.0-fabs(2.0*n/(N-1)-1.0);break;
 case WINDOW_KAISER:{double a=kb*sqrt(1.0-pow(2.0*n/(N-1)-1.0,2.0));w[i]=sinh(a)/sinh(kb);break;}
 default:w[i]=1.0;}}}
void fir_filter_apply(const double*h,int M,const double*x,int N,double*y){
 for(int n=0;n<N;n++){y[n]=0.0;for(int k=0;k<=M&&k<=n;k++)y[n]+=h[k]*x[n-k];}}
void iir_filter_apply(const double*b,int nb,const double*a,int na,const double*x,int N,double*y){
 for(int n=0;n<N;n++){y[n]=0.0;for(int k=0;k<=nb&&k<=n;k++)y[n]+=b[k]*x[n-k];
  for(int k=1;k<=na&&k<=n;k++)y[n]-=a[k-1]*y[n-k];}}
void moving_avg_init(MovingAvgFilter*maf,int ws){maf->window_size=ws;maf->buffer=(double*)calloc((size_t)ws,sizeof(double));
 maf->buffer_index=0;maf->sum=0.0;}
double moving_avg_update(MovingAvgFilter*maf,double x){maf->sum+=x-maf->buffer[maf->buffer_index];
 maf->buffer[maf->buffer_index]=x;maf->buffer_index=(maf->buffer_index+1)%maf->window_size;
 return maf->sum/maf->window_size;}
void exp_moving_avg_init(ExpMovingAvg*ema,double alpha,double fs){ema->alpha=alpha;ema->sampling_rate=fs;ema->y_prev=0.0;}
double exp_moving_avg_update(ExpMovingAvg*ema,double x){double y=ema->alpha*x+(1.0-ema->alpha)*ema->y_prev;ema->y_prev=y;return y;}
void digital_resonator_init(DigitalResonator*res,double r,double w0){res->r=r;res->omega0=w0;res->y1=0.0;res->y2=0.0;}
double digital_resonator_update(DigitalResonator*res,double x){double y=2.0*res->r*cos(res->omega0)*res->y1-res->r*res->r*res->y2+x;
 res->y2=res->y1;res->y1=y;return y;}
void fir_filter_free(FIRFilter*f){if(!f)return;free(f->coefficients);}
void iir_filter_free(IIRFilter*f){if(!f)return;free(f->b_coeffs);free(f->a_coeffs);}
void moving_avg_free(MovingAvgFilter*maf){if(!maf)return;free(maf->buffer);}
