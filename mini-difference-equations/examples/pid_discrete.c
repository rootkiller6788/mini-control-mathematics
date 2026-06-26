#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffeq.h"
#include "signal.h"
#include "stability.h"
int main(void) {
    printf("=== Discrete PID Controller & Engineering Applications ===\n\n");
    printf("1. Discrete PID step response:\n");
    DiscretePID pid;discrete_pid_init(&pid,2.0,0.5,0.1,0.01,-10.0,10.0);
    printf("   Kp=2.0, Ki=0.5, Kd=0.1, Ts=0.01\n");
    printf("   t   |  SP  |  PV  |  Output\n");
    double pv=0.0;
    for(int i=0;i<=50;i+=5){double sp=(i>=5)?1.0:0.0;
        double u=discrete_pid_update(&pid,sp,pv);
        pv+=0.1*(u-pv);/* simple 1st-order plant */
        printf("  %4.2f | %3.1f  | %5.3f | %7.4f\n",i*0.01,sp,pv,u);}
    printf("\n2. RC Low-pass Filter (R=1kOhm, C=100uF, Ts=1ms):\n");
    RCLowPass rc;rc_lowpass_init(&rc,1000.0,100e-6,0.001);
    printf("   t(ms) | Vin | Vout\n");
    double vc=0.0;
    for(int i=0;i<=50;i+=5){
        double vin=(i>=2)?5.0:0.0;vc=rc_lowpass_update(&rc,vin);
        printf("   %4d  | %3.0f  | %6.4f\n",i,vin,vc);}
    printf("\n3. DC Motor (J=0.01, B=0.1, Kt=0.05):\n");
    DCMotorDisc motor;dc_motor_discrete_init(&motor,0.01,0.1,0.05,0.001);
    printf("   t(ms) | Vin(V) |  w(rad/s)\n");
    double w=0.0;
    for(int i=0;i<=100;i+=10){w=dc_motor_discrete_update(&motor,12.0,2.0);
        printf("   %4d  |   12   | %8.4f\n",i,w);}
    printf("\n4. Newton Cooling (k=0.1, Tamb=25C, T0=100C, Ts=0.1s):\n");
    NewtonCooling nc;newton_cooling_init(&nc,0.1,25.0,0.1);
    printf("   t(s) |  T(C)\n");
    double T=100.0;
    for(int i=0;i<=100;i+=10){printf("   %3.0f  | %6.2f\n",i*0.1,T);T=newton_cooling_update(&nc,T);}
    printf("\n5. Quadcopter Altitude (m=1kg, k_thrust=20, Ts=0.01s, throttle=0.6):\n");
    QuadAltDisc quad;quad_alt_disc_init(&quad,1.0,20.0,0.01);
    double z=0.0,v=0.0;
    printf("   t(s) |  z(m)  |  v(m/s)\n");
    for(int i=0;i<=100;i+=10){if(i>=10)quad_alt_disc_update(&quad,0.6,&z,&v);
        printf("  %4.2f | %6.3f | %7.3f\n",i*0.01,z,v);}
    return 0;
}
