#include "kernel.h"
#include "math_mat.h" // -> transpose
#include <stdio.h>
#include <stdlib.h>
int euclidean_division(double a, double b)
{
    if ((b == 0) | (a == 0))
        return 0;
    int sign = ((a < 0) ^ (b < 0)) ? -1 : 1;
    int q = 0;
    b = (b < 0) ? -b : b;
    double r = (a < 0) ? -a : a;
    if (b > r && sign < 0)
    {
        double tmp = r;
        r = b;
        b = tmp;
    }
    while (r >= b)
    {
        r = r - b;
        q += 1;
    }
    // printf("%f // %f = %d\n", a, r, sign * q);
    return sign * q;
}
double modulo(double x, double y)
{
    /*x modulo y*/
    if (x == 0)
        return 0;
    int sign = ((x < 0) ^ (y < 0)) ? -1 : 1;
    x = x - sign * y * abs(euclidean_division(x, y));
    return x;
}
void echanger2(double *A, int j1, int j2, double *B, int nA, int nB)
{
    double tmp;
    for (int i = 0; i < nA; i++)
    {
        tmp = A[i * nA + j1];
        A[i * nA + j1] = A[i * nA + j2];
        A[i * nA + j2] = tmp;
    }
    for (int i = 0; i < nB; i++)
    {
        tmp = B[i * nB + j1];
        B[i * nB + j1] = B[i * nB + j2];
        B[i * nB + j2] = tmp;
    }
}
double *transpose(double *A, int n, int m)
{
    if (n == 0)
        return NULL;
    double *B = (double *)calloc(n * n, sizeof(double));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            B[i * n + j] = A[i + j * n];
        }
    }
    return B;
}
void mult2(double *A, int j, double t, double *B, int nA, int nB)
{
    for (int i = 0; i < nA; i++)
    {
        A[i * nA + j] = A[i * nA + j] * t;
    }
    for (int i = 0; i < nB; i++)
    {
        B[i * nB + j] = B[i * nB + j] * t;
    }
}
void div2(double *A, int j, double t, double *B, int nA, int nB)
{
    for (int i = 0; i < nA; i++)
    {
        A[i * nA + j] = euclidean_division(A[i * nA + j], t);
    }
    for (int i = 0; i < nB; i++)
    {
        B[i * nB + j] = euclidean_division(B[i * nB + j], t);
    }
}
void add2(double *A, int j1, double t, int j2, double *B, int nA, int nB)
{
    for (int i = 0; i < nA; i++)
    {
        A[i * nA + j1] = A[i * nA + j1] + t * A[i * nA + j2];
    }
    for (int i = 0; i < nB; i++)
    {
        B[i * nB + j1] = B[i * nB + j1] + t * B[i * nB + j2];
    }
}
double gcd(double a, double b)
{
    // printf("%f mod %f = ", a, b);
    double tmp;
    while (b != 0)
    {
        tmp = a;
        a = b;
        // printf("%f mod %f = ", tmp, b);
        b = modulo(tmp, b);
        // printf("%f\n", b);
    }
    return a;
}
void divide(double *A, int j, double t, int n)
{
    for (int i = 0; i < n; i++)
    {
        A[i * n + j] = euclidean_division(A[i * n + j], t);
    }
}
double gcdcol(double *A, int j, int n)
{
    double d = 0;
    for (int i = 0; i < n; i++)
    {
        // printf("%f %f \n",d,A[i * n + j]);
        // printf("gcd(%f,%f)=", d, A[i * n + j]);
        d = gcd(d, A[i * n + j]);
        // printf("%f\n", d);
    }
    return d;
}
void normaliser_colonne(double *A, int j, int n)
{
    double d = gcdcol(A, j, n);
    if (d != 0)
    {
        divide(A, j, d, n);
    }
}
void pivoter2(double *A, int l, int k, double *B, int nA, int nB)
{
    double t = A[l * nA + k];
    for (int j = k + 1; j < nA; j++)
    {
        double u = A[l * nA + j];
        // printf("%f %f\n", t, u);
        double d = gcd(t, u);
        // printf("%f rat\n", d);
        for (int i = 0; i < nA; i++)
        {
            // printf("%f %f connasse\n", u, d);

            A[i * nA + j] = euclidean_division(t, d) * A[i * nA + j] - euclidean_division(u, d) * A[i * nA + k];
        }
        for (int i = 0; i < nB; i++)
        {
            B[i * nB + j] = euclidean_division(t, d) * B[i * nB + j] - euclidean_division(u, d) * B[i * nB + k];
        }
    }
}
pivot_t chercher_pivot(double *A, int k, int n)
{
    pivot_t pivot;
    pivot.i = -1;
    pivot.j = -1;
    pivot.ok = false;
    for (int i = k; i < n; i++)
    {
        for (int j = k; j < n; j++)
        {
            if (A[i * n + j] != 0)
            {
                pivot.i = i;
                pivot.j = j;
                pivot.ok = true;
                return pivot;
            }
        }
    }
    return pivot;
}
double *eye(int n)
{
    double *tmp = (double *)calloc(n * n, sizeof(double));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (j == i)
            {
                tmp[i * n + j] = 1;
            }
        }
    }
    return tmp;
}
imker_t imker(double *A, int n)
{
    double *copieA = (double *)malloc(n * n * sizeof(double));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            copieA[i * n + j] = A[i * n + j];
        }
    }
    double *B = eye(n);
    int nB = n;
    int k = 0;
    while (true)
    {
        pivot_t res = chercher_pivot(copieA, k, n);
        if (!res.ok)
            break;
        if (res.j != k)
            echanger2(copieA, res.j, k, B, n, nB);

        // printf("%d %d\n", res.i, k);
        // __repr__(B,4,4);
        pivoter2(copieA, res.i, k, B, n, nB);
        k++;
    }
    //   __repr__(B,4,4);
    for (int i = 0; i < n; i++)
    {
        normaliser_colonne(copieA, i, n);
    }
    for (int i = 0; i < n; i++)
    {
        normaliser_colonne(B, i, nB);
    }
    imker_t ret;
    ret.A = copieA;
    ret.B = B;
    return ret;
}
double *col(double *A, int j, int n)
{
    double *res = (double *)calloc(n * n, sizeof(double));
    for (int i = 0; i < n; i++)
    {
        res[i] = A[i * n + j];
    }
    return res;
}
bool est_col_nulle(double *A, int j, int n)
{
    for (int i = 0; i < n; i++)
    {
        if (A[i * n + j] != 0)
            return false;
    }
    return true;
}
col_t col_nulles(double *A, int n)
{
    int *cols = (int *)calloc(n, sizeof(int));
    int size = 0;
    for (int i = 0; i < n; i++)
    {
        if (est_col_nulle(A, i, n))
        {
            cols[i] = i + 1;
            size++;
        }
    }
    int *col = (int *)calloc(size, sizeof(int));
    int iter = 0;
    for (int i = 0; i < n; i++)
    {
        if (cols[i] != 0)
        {
            col[iter] = cols[i] - 1;
            iter++;
        }
    }
    free(cols);
    col_t ret;
    ret.col = col;
    ret.size = size;
    return ret;
}
noyau_t noyau(double *A, int n)
{
    imker_t res = imker(A, n);
    col_t nulles = col_nulles(res.A, n);
    double *tmp;
    double *ret = (double *)malloc(n * n * sizeof(double));
    for (int i = 0; i < nulles.size; i++)
    {
        tmp = col(res.B, nulles.col[i], n);
        for (int j = 0; j < n; j++)
        {
            ret[i * n + j] = tmp[j];
        }
    }
    double *transp = transpose(ret, n, nulles.size);
    noyau_t ker;
    ker.vectors = transp;
    ker.size = nulles.size;
    free(ret);
    return ker;
}