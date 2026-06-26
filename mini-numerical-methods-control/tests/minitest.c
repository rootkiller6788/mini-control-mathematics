#include <stdio.h>
#include <stdlib.h>
#include "../include/numerical_core.h"
#include "../include/control_solvers.h"

int main() {
    printf("Test 1: controllability_matrix\n");
    double Ad[]={0,1,0,0}, Bd[]={0,1};
    Matrix *A=matrix_create_from(2,2,Ad), *B=matrix_create_from(2,1,Bd);
    Matrix *Ctrb=controllability_matrix(A,B);
    if(!Ctrb){printf("FAIL NULL\n");return 1;}
    printf("  Ctrb rows=%zu cols=%zu\n",Ctrb->rows,Ctrb->cols);
    printf("  Ctrb[0,0]=%f Ctrb[0,1]=%f\n",Ctrb->data[0],Ctrb->data[1]);
    printf("  Ctrb[1,0]=%f Ctrb->data[1]=%f\n",Ctrb->data[2],Ctrb->data[3]);
    matrix_free(Ctrb);
    
    printf("Test 2: is_controllable\n");
    int c=is_controllable(A,B,1e-6);
    printf("  controllable=%d\n",c);
    
    printf("Test 3: LQR design\n");
    double Qd[]={1,0,0,1}, Rd[]={1};
    Matrix *Q=matrix_create_from(2,2,Qd), *R=matrix_create_from(1,1,Rd);
    printf("  calling lqr_design...\n"); fflush(stdout);
    LQRResult *lqr=lqr_design(A,B,Q,R);
    printf("  lqr=%p\n",(void*)lqr);
    if(lqr){
        printf("  K[0]=%f K[1]=%f\n",lqr->K->data[0],lqr->K->data[1]);
        lqr_result_free(lqr);
    }
    matrix_free(A); matrix_free(B); matrix_free(Q); matrix_free(R);
    printf("All done!\n");
    return 0;
}
