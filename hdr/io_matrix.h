/* Contains IO utility functions */

#pragma once

#include <stdio.h>
#include <stdlib.h>

/* Returns a matrix red from the file "filename" */
double *ReadMat(const char* filename, int *n, int *m);

/* Returns a vector red from the file "filename" */
double *ReadVec(const char* filename, int *n);

/* Display in terminal the matrix "mat" with "n" ligns and "m" cols */
void PrintMat(double* mat, int n, int m);

/* Display in terminal the vector "x" of size "n" */
void PrintVec(double *x, int n);