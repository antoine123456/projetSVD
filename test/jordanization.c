#include <stdio.h>
#include <stdlib.h>
#include "kernel.h"
#include "jordan.h"
void matA(double *A)
{
    A[0] = 1;
    A[1] = -1;
    A[2] = 1;
    A[3] = -1;
    A[4] = 0;
    A[5] = -2;
    A[6] = 0;
    A[7] = -3;
    A[8] = 2;
    A[9] = 3;
    A[10] = 0;
    A[11] = 3;
    A[12] = 1;
    A[13] = 5;
    A[14] = -1;
    A[15] = 6;
}
void jordanTest(){
      int size = 4;
    double *A = (double *)calloc(size * size, sizeof(double));
    matA(A);
    double V[4] = {2, 1, 1, 1};
    Jordan_Bidiag_t res = jordan(A, 4, 4, V);
    __repr__(res.P, 4, 4);
    printf("  ======================================\n");
    __repr__(res.J, 4, 4);
}