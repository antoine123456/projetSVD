#include <SVD.h>
#include "trigonalisation.h"
#include <io_matrix.h>
#include <string.h>
#include<omp.h>
#include<stdlib.h>
#include<math_mat.h>
#include<gen_mat.h>


double* SVD_1(double *A, int m, int n) {
    double *B;  // will store A^TA or AA^T
    int nB;    // B dimension

    // Compute A^TA or AA^T depending on the comparaison between n and m

    if (n < m) {
        B = get_AtA(A, m, n);
        nB = n;
    } else {
        B = get_AAt(A, m, n);
        nB = m;
    }

    eigen_t Beigens = QR_method(B, nB);

    free(B);
    free(Beigens.vectors);

    // Beigen values are A singular values
    return Beigens.values;

}

double* SVD_1_Tridiag(double *A, int m, int n) {
    double *B;  // will store A^TA or AA^T
    int nB;    // B dimension
    eigen_t Beigens;
    // Compute A^TA or AA^T depending on the comparaison between n and m
/*	#pragma omp parallel private(B,nB,Beigens)
    {
	#pragma omp single
	    { */
    if (n < m) {
        B = get_AtA(A, m, n);
        nB = n;
    } else {
        B = get_AAt(A, m, n);
        nB = m;
    }

     Beigens = QR_method_Tridiag(B, nB);
//	}

    free(B);
    free(Beigens.vectors);
	//}
    // Beigen values are A singular values
    return Beigens.values;
}

// Call SVD_1 with A transformed into a bidiagonale matrix
double* SVD_3(double *A, int m, int n) {
    // Bidiagonalization (GKL)
    GKL_Bidiag_t b = Bidiagonalization(A, m, n);

    double* singular_values = SVD_1_Tridiag(b.B, m, n);

    freeGKL(b);

    return singular_values;
}

double* SVD_Hessenberg(double *A, int m, int n) {
    assert(A);
    assert(m>1);
    assert(n>1);
    //double *B ;  // will store A^TA or AA^T
  //  double *H ; // hessenberg form of B as B = P'HP
    int dimB;    // B dimension

    if(n<m){
      dimB = n;
    }else {
      dimB = m;
    }
    double *B = (double*) malloc(dimB*dimB*sizeof(double));
    test_malloc(B);

    get_symB(A,m,n,B,&dimB);

    Hess_Reduction2(B, dimB);

    eigen_t Beigens = QR_method_Tridiag(B, dimB);
    free(B);

    free(Beigens.vectors);

    // Beigen values are A singular values
    return Beigens.values;
}

void copy_A_to_B(double *A, int ma,int na, double *B, int mb,int nb,int start_i,int start_j){
  int offset = 0;
  int limit_row =start_i+mb;
  int limit_colum = start_j+nb;
  if(B==NULL) B =malloc(sizeof(double)*mb*nb);

  for( int i = start_i; i<limit_row ; i++ ){
    for( int j = start_j; j<limit_colum; j++){
      offset = i*na+j;
  //    printf("A[%d] = %lf\n",offset,A[offset]);
      B[i*nb+j]= A[offset];
  //    printf("B[%d] = %lf\n",i*nb+j,B[i*nb+j]);
    }
  }
}

void concatenate( double *array, double *even,double *odd, int ne, int no,int *index){
  int loop,  e_len, o_len;

  e_len =ne;
  o_len = no;



  for(loop = 0; loop < e_len; loop++) {
    array[*index] = even[loop];
    (*index)=(*index)+1;
  }

  for(loop = 0; loop < o_len; loop++) {
    array[*index] = odd[loop];
    (*index)=(*index)+1;
  }

}
