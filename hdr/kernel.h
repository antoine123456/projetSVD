#include <stdbool.h>
typedef struct
{
    int i;
    int j;
    bool ok;
} pivot_t;

void echanger2(double *A, int j1, int j2, double *B, int nA, int nB);
void mult2(double *A, int j, double t, double *B, int nA, int nB);
void div2(double *A, int j, double t, double *B, int nA, int nB);
void add2(double *A, int j1, double t, int j2, double *B, int nA, int nB);
double gcd(double a, double b);
void div(double *A, int j, double t, int n);
double gcdcol(double *A, int j, int n);
void normaliser_colonne(double *A, int j, int n);
void pivoter2(double *A, int l, int k, double *B, int nA, int nB);
pivot_t chercher_pivot(double *A, int k, int n);
double *eye(int n);
double *imker(double *A, int n);
double *transpose(double *A, int n);
double *col(double *A, int j, int n);
bool est_col_nulle(double *A, int j, int n);
int *col_nulles(double *A, int n);
double *noyau(double *A, int n);