// #pragma once
// #include <math_mat.h>
// #include <QR.h>
typedef struct
{
    double *J;
    double *P;
} Jordan_Bidiag_t;

typedef struct
{
    int nEigVal;
    int *mult;
    int *from;
} DiagInfo_t;

// pour Proto //////
typedef struct
{
    double *values;
    double *vectors;
} eigen_t;
////////////////////
void freeJordan(Jordan_Bidiag_t j);
void freeInfo(DiagInfo_t d);

int listIndex(double *list, int *vus, double val, int n);

void to_square_mat(double *Mat, int m, int n);

DiagInfo_t getMult(double *eigVal, int n);

Jordan_Bidiag_t jordan(double *Mat, int m, int n, double *V);