#include <stdbool.h>
typedef struct
{
    int i;
    int j;
    bool ok;
} pivot_t;

typedef struct
{
    double *A;
    double *B;
} imker_t;

typedef struct
{
    int *col;
    int size;
} col_t;

typedef struct
{
    double *vectors;
    int size;
} noyau_t;

void __repr__(double *A, int n, int m);
void defMat(double *A);
int euclidean_division(double a, double b);
double modulo(double x, double y);
void echanger2(double *A, int j1, int j2, double *B, int nA, int nB);
void mult2(double *A, int j, double t, double *B, int nA, int nB);
void div2(double *A, int j, double t, double *B, int nA, int nB);
void add2(double *A, int j1, double t, int j2, double *B, int nA, int nB);
double gcd(double a, double b);
void divide(double *A, int j, double t, int n);
double gcdcol(double *A, int j, int n);
void normaliser_colonne(double *A, int j, int n);
void pivoter2(double *A, int l, int k, double *B, int nA, int nB);
pivot_t chercher_pivot(double *A, int k, int n);
double *eye(int n);
imker_t imker(double *A, int n);
double *transpose(double *A, int n, int m);
double *col(double *A, int j, int n);
bool est_col_nulle(double *A, int j, int n);
col_t col_nulles(double *A, int n);
noyau_t noyau(double *A, int n);