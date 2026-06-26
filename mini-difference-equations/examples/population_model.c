#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffeq.h"
int main(void) {
    printf("=== Population Dynamics via Difference Equations ===\n\n");
    printf("1. Malthusian Growth: P[n+1] = (1+r) * P[n]\n");
    double pop_m[50];malthusian_growth(100.0,0.05,50,pop_m);
    printf("   P0=100, r=5%%, n=0..10:\n   ");
    for(int i=0;i<=10;i++)printf("%.0f ",pop_m[i]);printf("\n\n");
    printf("2. Logistic Map: x[n+1] = r * x[n] * (1 - x[n])\n");
    double r_values[]={2.0,3.2,3.5,3.8};int nv=4;
    for(int ri=0;ri<nv;ri++){double r=r_values[ri];double traj[100];
        logistic_map_iterate(r,0.5,100,traj);
        printf("   r=%.1f: last 20 values:\n   ",r);
        for(int i=80;i<100;i++)printf("%.4f ",traj[i]);printf("\n");}
    printf("\n3. Verhulst Logistic: P[n+1]=P[n]+r*P[n]*(1-P[n]/K)\n");
    double pop_v[100];verhulst_logistic_discrete(10.0,0.8,1000.0,100,pop_v);
    printf("   P0=10, r=0.8, K=1000, n=0,10,20,...:\n   ");
    for(int i=0;i<=100;i+=10)printf("%.1f ",pop_v[i]);printf("\n\n");
    printf("4. Lotka-Volterra Predator-Prey:\n");
    double prey[200],pred[200];prey[0]=40;pred[0]=9;
    lotka_volterra_discrete(0.1,0.02,0.01,0.01,prey,pred,200);
    printf("   t  |  prey  |  predator\n");
    for(int i=0;i<=200;i+=20)printf("  %3d | %6.1f | %6.1f\n",i,prey[i],pred[i]);
    printf("\n5. SIR Epidemic Model:\n");
    SEIRModel m;seir_model_init(&m,9900,100,0,0,0.3,0.2,0.1);
    printf("   day |    S    |    E    |    I    |    R\n");
    for(int d=0;d<=50;d+=5){printf("   %3d | %7.0f | %7.0f | %7.0f | %7.0f\n",
        d,m.S,m.E,m.Inf,m.R);for(int s=0;s<5;s++)seir_model_step(&m);}
    return 0;
}
