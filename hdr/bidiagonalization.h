#pragma once

#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <math_mat.h>

// B = U*AV Golub-Kahan-Lanczos Bidiagonalization Procedure
double *Bidiagonalization(double *A, int m, int n);