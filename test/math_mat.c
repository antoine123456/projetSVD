#include <assert.h>
#include <stdio.h>

#include <math_mat.h>

void Test_MatMul() {
    double A[4] = {1, 2, 3, 4};
    double B[4] = {4, 5, 6, 7};
    double C[4];

    MatMul(C, A, B, 2);

    assert(C[0] == 16);
    assert(C[1] == 19);
    assert(C[2] == 36);
    assert(C[3] == 43);
}

void Test_MatMulTrans() {
    double A[4] = {1, 2, 3, 4};
    double B[4] = {4, 5, 6, 7};
    double C[4];

    MatMulTrans(C, A, B, 2);

    assert(C[0] == 22);
    assert(C[1] == 26);
    assert(C[2] == 32);
    assert(C[3] == 38);
}

void Test_DotProdLC() {
    double A[4] = {1, 2, 3, 4};
    double B[4] = {4, 5, 6, 7};

    assert(DotProdLC(A, 0, B, 0, 2, 2) == 16.);
}

void Test_DotProdCC() {
    double A[6] = {1, 2, 3, 4, 5, 6};
    double B[4] = {1, 2, 3, 4};

    assert(DotProdCC(A, 1, 3, B, 0, 2, 2) == 17);
}


int main() {
    Test_MatMul();
    Test_MatMulTrans();
    Test_DotProdLC();
    Test_DotProdCC();

    return 0;
}