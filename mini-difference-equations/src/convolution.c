/* convolution.c - Discrete Convolution Sum and Correlation
 * Reference: Oppenheim & Schafer "Discrete-Time Signal Processing" (2010)
 * Reference: Proakis & Manolakis "Digital Signal Processing" (2013)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include "diffeq.h"
/* L2: Direct convolution y[n] = Σ x[k] h[n-k] */
int convolve_discrete(const double*x,int Nx,const double*h,int Nh,double*y){
 if(!x||!h||!y)return -1;int Ny=Nx+Nh-1;
 for(int n=0;n<Ny;n++){double s=0.0;for(int k=0;k<=n;k++){
  double xk=(k<Nx)?x[k]:0.0;double hnk=((n-k)<Nh&&(n-k)>=0)?h[n-k]:0.0;s+=xk*hnk;}
  y[n]=s;}return Ny;}
/* L5: FFT-based convolution (radix-2 Cooley-Tukey) */
static void fft_radix2(double complex*x,int N,int inv){
 /* Bit-reversal permutation */
 for(int i=1,j=0;i<N;i++){int bit=N>>1;for(;j&bit;bit>>=1)j^=bit;j^=bit;if(i<j){double complex t=x[i];x[i]=x[j];x[j]=t;}}
 for(int len=2;len<=N;len<<=1){double ang=2.0*M_PI/len*(inv?-1:1);
  double complex wlen=cos(ang)+sin(ang)*I;
  for(int i=0;i<N;i+=len){double complex w=1.0;
   for(int j=0;j<len/2;j++){double complex u=x[i+j],v=x[i+j+len/2]*w;
    x[i+j]=u+v;x[i+j+len/2]=u-v;w*=wlen;}}}
 if(inv)for(int i=0;i<N;i++)x[i]/=N;}
int convolve_fft(const double*x,int Nx,const double*h,int Nh,double*y){
 if(!x||!h||!y)return -1;int N=1;while(N<Nx+Nh-1)N<<=1;
 double complex*X=(double complex*)calloc((size_t)N,sizeof(double complex));
 double complex*H=(double complex*)calloc((size_t)N,sizeof(double complex));
 for(int i=0;i<Nx;i++)X[i]=x[i]+0.0*I;for(int i=0;i<Nh;i++)H[i]=h[i]+0.0*I;
 fft_radix2(X,N,0);fft_radix2(H,N,0);
 for(int i=0;i<N;i++)X[i]*=H[i];fft_radix2(X,N,1);
 int Ny=Nx+Nh-1;for(int i=0;i<Ny;i++)y[i]=creal(X[i]);free(X);free(H);return Ny;}
/* L5: Deconvolution via long division */
int deconvolve_discrete(const double*y,int Ny,const double*h,int Nh,double*x){
 if(!y||!h||!x)return -1;/* x = deconv(y,h) where y = x * h */
 /* Solve Toeplitz system using Levinson recursion (simplified one-step) */
 if(Nh==1&&fabs(h[0])>1e-15){double inv=1.0/h[0];
  for(int i=0;i<Ny;i++)x[i]=y[i]*inv;return Ny-Nh+1;}
 /* General: solve lower triangular system */
 int Nx=Ny-Nh+1;if(Nx<1)return 0;
 for(int n=0;n<Nx;n++){double s=y[n];
  for(int k=1;k<=n&&k<Nh;k++)s-=x[n-k]*h[k];
  if(fabs(h[0])>1e-15)x[n]=s/h[0];else x[n]=0.0;}return Nx;}
/* L3: Autocorrelation Rxx[k] = Σ x[n] x[n+k] */
void auto_correlation_discrete(const double*x,int N,double*Rxx,int max_lag){
 if(!x||!Rxx)return;for(int k=0;k<=max_lag&&k<N;k++){Rxx[k]=0.0;
  for(int n=0;n<N-k;n++)Rxx[k]+=x[n]*x[n+k];}}
/* L3: Cross-correlation Rxy[k] = Σ x[n] y[n+k] */
void cross_correlation_discrete(const double*x,const double*y,int N,double*Rxy,int max_lag){
 if(!x||!y||!Rxy)return;for(int k=-max_lag;k<=max_lag;k++){Rxy[k+max_lag]=0.0;
  for(int n=0;n<N;n++){int nk=n+k;if(nk>=0&&nk<N)Rxy[k+max_lag]+=x[n]*y[nk];}}}
/* L3: Circular convolution */
int circular_convolve_discrete(const double*x,const double*h,int N,double*y){
 if(!x||!h||!y)return -1;
 for(int n=0;n<N;n++){y[n]=0.0;for(int k=0;k<N;k++)y[n]+=x[k]*h[(n-k+N)%N];}return N;}
/* L4: Convolution properties verification */
int convolution_properties_verify(const double*x,const double*h1,const double*h2,int Nx,int Nh,double tol){
 if(!x||!h1||!h2)return 0;int Ny=Nx+Nh-1;
 double*y1=(double*)calloc((size_t)Ny,sizeof(double));
 double*y2=(double*)calloc((size_t)Ny,sizeof(double));
 double*y3=(double*)calloc((size_t)Ny,sizeof(double));
 double*y4=(double*)calloc((size_t)Ny,sizeof(double));
 double*y5=(double*)calloc((size_t)Ny,sizeof(double));
 /* Commutative: x*h1 = h1*x */
 convolve_discrete(x,Nx,h1,Nh,y1);convolve_discrete(h1,Nh,x,Nx,y2);
 int comm=1;for(int i=0;i<Ny;i++)if(fabs(y1[i]-y2[i])>tol){comm=0;break;}
 /* Associative: (x*h1)*h2 = x*(h1*h2) */
 int N12=Nx+Nh-1;double*h12=(double*)calloc((size_t)N12,sizeof(double));
 convolve_discrete(x,Nx,h1,Nh,y3);convolve_discrete(y3,Ny,h2,Nh,y4);
 convolve_discrete(h1,Nh,h2,Nh,h12);convolve_discrete(x,Nx,h12,N12,y5);
 int assoc=1;int Nmax=Ny+Nh-1;for(int i=0;i<Nmax;i++){double v4=(i<Ny+Nh-1)?y4[i]:0.0,
  v5=(i<Nx+N12-1)?y5[i]:0.0;if(fabs(v4-v5)>tol){assoc=0;break;}}
 /* Distributive: x*(h1+h2) = x*h1 + x*h2 */
 double*hsum=(double*)calloc((size_t)(Nh>Nh?Nh:Nh),sizeof(double));
 for(int i=0;i<Nh;i++)hsum[i]=h1[i]+h2[i];
 convolve_discrete(x,Nx,hsum,Nh,y3);convolve_discrete(x,Nx,h1,Nh,y4);
 convolve_discrete(x,Nx,h2,Nh,y5);int dist=1;
 for(int i=0;i<Ny;i++)if(fabs(y3[i]-(y4[i]+y5[i]))>tol){dist=0;break;}
 free(y1);free(y2);free(y3);free(y4);free(y5);free(h12);free(hsum);
 return comm&&assoc&&dist;}
/* L3: Toeplitz matrix from impulse response */
void toeplitz_from_impulse_response(const double*h,int Nh,int N,double**T){
 if(!h||!T)return;for(int i=0;i<N;i++)for(int j=0;j<N;j++){int idx=i-j;
  T[i][j]=(idx>=0&&idx<Nh)?h[idx]:0.0;}}
/* L3: Overlap-add convolution for long signals */
int overlap_add_convolution(const double*x,int Nx,const double*h,int Nh,int block_size,double*y){
 if(!x||!h||!y||block_size<1)return -1;int Ny=Nx+Nh-1;
 for(int i=0;i<Ny;i++)y[i]=0.0;/* Zero output */
 int L=block_size;int M=Nh;/* Filter length */
 for(int block_start=0;block_start<Nx;block_start+=L){
  int block_len=L;if(block_start+block_len>Nx)block_len=Nx-block_start;
  double*block_y=(double*)calloc((size_t)(block_len+M-1),sizeof(double));
  convolve_discrete(x+block_start,block_len,h,M,block_y);
  for(int i=0;i<block_len+M-1;i++)if(block_start+i<Ny)y[block_start+i]+=block_y[i];
  free(block_y);}
 return Ny;}
/* L3: Overlap-save convolution */
int overlap_save_convolution(const double*x,int Nx,const double*h,int Nh,int block_size,double*y){
 if(!x||!h||!y||block_size<Nh)return -1;int Ny=Nx+Nh-1;
 for(int i=0;i<Ny;i++)y[i]=0.0;int L=block_size-Nh+1;/* Valid output per block */
 for(int block_start=0;block_start<Nx;block_start+=L){
  int input_start=(block_start-Nh+1>=0)?block_start-Nh+1:0;
  int input_len=block_size;if(input_start+input_len>Nx)input_len=Nx-input_start;
  double*block_in=(double*)calloc((size_t)input_len,sizeof(double));
  for(int i=0;i<input_len;i++)block_in[i]=(input_start+i>=0&&input_start+i<Nx)?x[input_start+i]:0.0;
  double*block_y=(double*)calloc((size_t)(input_len+Nh-1),sizeof(double));
  convolve_discrete(block_in,input_len,h,Nh,block_y);/* Save only valid part */
  for(int i=0;i<L;i++){int out_idx=block_start+i;if(out_idx<Ny)y[out_idx]=block_y[Nh-1+i];}
  free(block_in);free(block_y);}
 return Ny;}
