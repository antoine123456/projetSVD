#include <math.h>
#include <string.h>

#include <gram_schmidt.h>
#include <io_matrix.h>
#include "test_mat.h"
#include "test_utils.h"

int nerr = 0;

void Test_GramSchmidt(void (*f)(QR_t*, double*, int), int n, double eps, const char *name) {
    QR_t QR;
    QR.Q = (double *) malloc(sizeof(double )* n*n);
    QR.R = (double *) malloc(sizeof(double )* n*n);
    double *C = (double*) malloc(sizeof(double) * n*n);
    char test_name[30] = "QR Decomposition ";
    strcat(test_name, name);

    double *givens_mat = givens(n,n);
    f(&QR, givens_mat, n);
    MatMul(C, QR.Q, QR.R, n);
    for (int i=0 ; i<n*n ; i++)
        if (fabs(C[i] - givens_mat[i]) > eps) {
            Print_Error_Mat(test_name, "Givens");
            nerr++;
            break;
        }

    double *aegerter_mat = aegerter(n);
    f(&QR, aegerter_mat, n);
    MatMul(C, QR.Q, QR.R, n);
    for (int i=0 ; i<n*n ; i++)
        if (fabs(C[i] - aegerter_mat[i]) > eps) {
            Print_Error_Mat(test_name, "Aegerter");
            nerr++;
            break;
        }

    double *orth_symm_mat = orth_symm(n);
    f(&QR, orth_symm_mat, n);
    MatMul(C, QR.Q, QR.R, n);
    for (int i=0 ; i<n*n ; i++)
        if (fabs(C[i] - orth_symm_mat[i]) > eps) {
            Print_Error_Mat(test_name, "Orth_Symm");
            nerr++;
            break;
        }
}

void Test_GramSchmidt_Tridiag() {
    // Test parameters
    int n = 10;        // matrix order
    double eps = 0.01; // error threshold

    QR_t QR;
    QR.Q = (double *) malloc(sizeof(double )* n*n);
    QR.R = (double *) malloc(sizeof(double )* n*n);
    double *C = (double*) malloc(sizeof(double) * n*n);

    double *bab_mat = bab(n, 5, 2);
    GramSchmidtMod_Tridiag(&QR, bab_mat, n);
    MatMul(C, QR.Q, QR.R, n);
    for (int i=0 ; i<n*n ; i++)
        if (fabs(C[i] - bab_mat[i]) > eps) {
            Print_Error_Mat("QR Decomposition (Tridiag)", "Bab");
            nerr++;
            break;
        }   

    double *lesp_mat = lesp(n, n);
    GramSchmidtMod_Tridiag(&QR, lesp_mat, n);
    MatMul(C, QR.Q, QR.R, n);
    for (int i=0 ; i<n*n ; i++)
        if (fabs(C[i] - lesp_mat[i]) > eps) {
            Print_Error_Mat("QR Decomposition (Tridiag)", "Lesp");
            nerr++;
            break;
        }   
}

int main() {
    int n = 50;
    double eps = 0.01;

    Test_GramSchmidt(GramSchmidt, n, eps, "CGS");
    Test_GramSchmidt(GramSchmidtMod, n, eps, "MGS");

    Test_GramSchmidt_Tridiag();

    if (nerr == 0)
        Print_Success("QR Decomposition");

    return 0;
}