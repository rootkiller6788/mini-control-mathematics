#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "signal.h"
int main(void) {
    printf("=== Digital Filter Design via Difference Equation ===\n\n");
    FIRFilter fir;design_fir_filter(FILTER_LOWPASS,20,0.25,0.0,WINDOW_HAMMING,&fir);
    printf("FIR Low-pass filter (order=%d, fc=0.25*fs, Hamming window):\n",fir.order);
    printf("Coefficients h[0..%d]:\n",fir.order);
    for(int i=0;i<=fir.order;i++)printf("  h[%2d] = %+.6f\n",i,fir.coefficients[i]);
    printf("\nFilter test: input = sin(0.1*pi*n) + sin(0.4*pi*n)\n");
    double x[50],y[50];
    for(int i=0;i<50;i++)x[i]=sin(0.1*M_PI*i)+0.5*sin(0.4*M_PI*i);
    fir_filter_apply(fir.coefficients,fir.order,x,50,y);
    printf("n   |  x[n]     |  y[n]\n");
    printf("----|-----------|----------\n");
    for(int i=0;i<15;i++)printf("%3d | %+.4f  | %+.4f\n",i,x[i],y[i]);
    printf("...\n");fir_filter_free(&fir);
    IIRFilter iir;design_iir_lowpass(100.0,1000.0,2,&iir);
    printf("\nIIR Butterworth Low-pass (order=2, fc=100Hz, fs=1kHz):\n");
    printf("b = [%.4f, %.4f, %.4f]\n",iir.b_coeffs[0],iir.b_coeffs[1],iir.b_coeffs[2]);
    printf("a = [1, %.4f, %.4f]\n",iir.a_coeffs[0],iir.a_coeffs[1]);
    iir_filter_free(&iir);
    return 0;
}
