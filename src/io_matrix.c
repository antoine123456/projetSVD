#include <io_matrix.h>

/* Returns a matrix red from the file "filename" */
double *ReadMat(const char* filename, int *n, int *m) {
    // Open file
    FILE *fp = fopen(filename, "r");

    // Read header : matrix dimensions (n,m)
    fscanf(fp, "%d", n);
    fscanf(fp, "%d", m);

    // Declaration
    double *mat = (double*) malloc(sizeof(double) * (*n)*(*m));

    // Read mat elements
    double el;
    int k=0;
    while (fscanf(fp, "%lf", &el) != EOF)
        mat[k++] = el;

    return mat;
}

/* Returns a vector red from the file "filename" */
double *ReadVec(const char* filename, int *n) {
    // Open file
    FILE *fp = fopen(filename, "r");

    // Read header : vector size
    fscanf(fp, "%d", n);

    // Declaration
    double *vec = (double*) malloc(sizeof(double) * (*n));

    // Read vec elements
    double el;
    int k=0;
    while (fscanf(fp, "%lf", &el) != EOF)
        vec[k++] = el;

    return vec;
}

/* Display in terminal the matrix "mat" with "n" ligns and "m" cols */
void PrintMat(double* mat, int n, int m) {
    for (int i=0 ; i<n ; i++) {
        for (int j=0 ; j<m ; j++)
            printf("%lf ", mat[i*m + j]);
        printf("\n");
    }
}

/* Display in terminal the vector "x" of size "n" */
void PrintVec(double *x, int n) {
    for (int i=0 ; i<n ; i++)
        printf("%lf\n", x[i]);
}