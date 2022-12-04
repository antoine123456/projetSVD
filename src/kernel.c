#include <kernel.h>
int euclidean_division(double a, double b)
{
    int q;
    double r = a;
    while (r >= b)
    {
        r = r - b;
        q += 1;
    }
    return q;
}
double modulo(double x, double y)
{
    /*x modulo y*/
    x -= y * abs(euclidean_division(x, y));
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
    double tmp;
    while (b != 0)
    {
        tmp = a;
        a = b;
        b = modulo(tmp, b);
    }
}
void div(double *A, int j, double t, int n)
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
        d = gcd(d, A[i * n + j]);
    }
    return d;
}
void normaliser_colonne(double *A, int j, int n)
{
    double d = gcdcol(A, j, n);
    if (d != 0)
    {
        div(A, j, d, n);
    }
}
void pivoter2(double *A, int l, int k, double *B, int nA, int nB)
{
    double t = A[l * nA + k];
    for (int j = k + 1; nA; j++)
    {
        double u = A[l * nA + j];
        double d = gcd(t, u);
        for (int i = 0; i < nA; i++)
        {
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
    double *tmp = (double *)calloc(n, sizeof(double));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (j == i)
            {
                tmp[i * n + j] = 0;
            }
        }
    }
    return tmp;
}
double *imker(double *A, int n)
{
    double *copieA = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            copieA[i * n + j] = A[i * n + j];
        }
    }
    double *B = eye(n);
    int k = 0;
    while (true)
    {
        pivot_t res = chercher_pivot(A, k, n);
        if (!res.ok)
        {
            break;
        }
    }
}
double *transpose(double *A, int n);
double *col(double *A, int j, int n);
bool est_col_nulle(double *A, int j, int n);
int *col_nulles(double *A, int n);
double *noyau(double *A, int n);