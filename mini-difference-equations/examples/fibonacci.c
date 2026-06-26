#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffeq.h"
int main(void) {
    printf("=== Fibonacci via Difference Equation ===\n\n");
    double F_closed, F_closed_prev;
    for(int n=0;n<=20;n++){fibonacci_closed_form(n,&F_closed);
        printf("F[%2d] = %.0f",n,F_closed);
        if(n>=1){fibonacci_closed_form(n-1,&F_closed_prev);
            printf("  (ratio F[n]/F[n-1] = %.6f)",F_closed/F_closed_prev);}
        printf("\n");}
    printf("\nGolden ratio phi = (1+sqrt(5))/2 = %.10f\n",(1.0+sqrt(5.0))/2.0);
    printf("As n->infinity, F[n]/F[n-1] -> phi\n");
    printf("\nUsing iterative difference equation:\n");
    double values[21];fibonacci_iterative(20,values,21);
    printf("F[0..20]: ");for(int i=0;i<=20;i++)printf("%.0f ",values[i]);printf("\n");
    return 0;
}
