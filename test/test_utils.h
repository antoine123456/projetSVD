#pragma once

#include <stdio.h>
#include <limits.h>
#include <math.h>

int FindEl(double *v, double val, double eps, int n) {
    double errmin = INT_MAX;
    for (int i=0 ; i<n ; i++) {
        double err = fabs(v[i] - val);
        if (err < errmin)
            errmin = err;
    }

    // err_mean += errmin;

    if (errmin < eps)
        return 1;
    else
        return 0;
}

int Test_Eigenvalues(
    double *known_ev,
    double *ev,
    double err_threshold,
    int n) {
    
    int passed = 1;
    for (int i=0 ; i<n ; i++)
        if (!FindEl(ev, known_ev[i], err_threshold, n))
            passed = 0;

    return passed;
}

// Transform eigenvalues vector to singular values vector
void Eigen_to_Singular(double *ev, int n) {
    for (int i=0 ; i<n ; i++)
        ev[i] = ev[i]*ev[i];
}

void Print_Error(const char *category) {
    printf("\033[0;31m");
    printf("- ERROR (%s)\n", category);
    printf("\033[0m");
}

void Print_Error_Mat(const char *category, const char *matrix) {
    printf("\033[0;31m");
    printf("- ERROR (%s) (Matrix : %s)\n", category, matrix);
    printf("\033[0m");
}

void Print_Success(const char *category) {
    printf("\033[0;32m");
    printf("- PASSED (%s)\n", category);
    printf("\033[0m");
}