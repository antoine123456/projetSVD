#include <SVD.h>
#include "trigonalisation.h"
#include <io_matrix.h>
#include <errno.h>
#include <string.h>
//#include<omp.h>
#include<stdlib.h>

#ifndef _TEST_MALLOC_
#define _TEST_MALLOC_
void test_malloc(double *buff){
  int errornum;
  if (buff == NULL) {
      errornum = errno;
      fprintf(stderr, "The Value of errno: %d\n", errno);
      perror("Error message that is printed by perror");
      fprintf(stderr, "Error message for allocation of B : %s\n", strerror( errornum ));
  }
}
#endif

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

    // Compute A^TA or AA^T depending on the comparaison between n and m
    if (n < m) {
        B = get_AtA(A, m, n);
        nB = n;
    } else {
        B = get_AAt(A, m, n);
        nB = m;
    }

    eigen_t Beigens = QR_method_Tridiag(B, nB);

    free(B);
    free(Beigens.vectors);

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
    double *B ;  // will store A^TA or AA^T
    double *H ; // hessenberg form of B as B = P'HP
    double *P; // householder matrix associate : B  =P'HP
    int nB;    // B dimension
    int errornum;


    //PrintMat(A,m,n);
  // Compute A^TA or AA^T depending on the comparaison between n and m
    if (n < m) {
        //printf("A:\n");
      //  B = get_AtA(A, m, n);
        B = (double*) malloc(n*n*sizeof(double));
               test_malloc(B);
      //  double *C =NULL;
        for (int i=0; i<n ; i++)
            for (int j=0 ; j<n ; j++) {
                double r = 0;
                for (int k=0 ; k<m ; k++){
                    r += A[k*n + i] * A[k*n + j];
                  }
                B[i*n + j] = r;
            }
        nB = n;

    } else {
       printf("A:\n");
       B = (double*) malloc(m*m*sizeof(double));
       test_malloc(B);

       for (int i=0; i<m ; i++)
           for (int j=0 ; j<m ; j++) {
               double r = 0;
               for (int k=0 ; k<n ; k++){
                   r += B[i*n + k] * B[j*n + k];
                 }
               B[i*m + j] = r;
           }
      //  B = get_AAt(A, m, n);
        nB = m;
    }
    H = calloc(nB*nB,sizeof(double*));
    test_malloc(H);

    P = calloc(nB*nB,sizeof(double*));
    test_malloc(P);

    Hess_Reduction(B, nB, H,P);
  //  printf("H:\n");
  //  PrintMat(H,nB,nB);
  //  printf("P:\n");
//    PrintMat(P,nB,nB);

    eigen_t Beigens = QR_method_Tridiag(H, nB);
  //  printf("H apres QR:\n");
  //  PrintMat(H,nB,nB);


    free(B);
    free(H);
    free(P);
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

void concatenate( double *array, double *even,double *odd, int ne, int no,int index){
  int loop,  e_len, o_len;

  e_len =ne;
  o_len = no;

  index = 0;

  for(loop = 0; loop < e_len; loop++) {
    array[index] = even[loop];
    index++;
  }

  for(loop = 0; loop < o_len; loop++) {
    array[index] = odd[loop];
    index++;
  }

}
double* SVD_blocks(double *M, int m, int n) {
   int na,ma,nb,mb,nc,mc,nd,md;
   double* A;
   double* B;
   double* C;
   double* D;
   int si,sj;
   if (m%2==0){
     ma = m/2;
     mb=ma;
   }else{
     ma = (m+1)/2;
     mb = (m-1)/2;
   }
   if (n%2==0){
     na = n/2;
     nb=na;
   }else{
     na = (n+1)/2;
     nb = (n-1)/2;
   }
   mc=ma;
   nc=na;
   md=mb;
   nd=nb;
   double *singvalsA;
   double *singvalsB;
   double *singvalsC;
   double *singvalsD;
   double *singvals;
   /*
   singvals = malloc(sizeof(double)*m*4);
  singvalsA = malloc(sizeof(double)*ma);
  singvalsB = malloc(sizeof(double)*mb);
  singvalsC = malloc(sizeof(double)*mc);
  singvalsD = malloc(sizeof(double)*md); */

   int index = 0;

//   omp_set_nested(1);
//  #pragma omp parallel sections num_threads(4) private(A,B,C,D)
  //{
    //  #pragma omp section
    //  {
        A =malloc(sizeof(double)*ma*na);
        copy_A_to_B(M,m,n,A,ma,na,0,0);
        singvalsA = SVD_1(A,  ma,  na);
        PrintVec(singvalsA,index);
        free(A);
    //  }

    //  #pragma omp section
    //  {
        B =malloc(sizeof(double)*mb*nb);
        copy_A_to_B(M,m,n,B,mb,nb,0,na);
        singvalsB = SVD_1(B,  mb, nb);
        PrintVec(singvalsB,index);
        free(B);
    //  }


    //  #pragma omp section
    //  {
        C =malloc(sizeof(double)*mc*nc);
        copy_A_to_B(M,m,n,C,mc,nc,ma,0);
        singvalsC = SVD_1(C,  mc, nc);
        PrintVec(singvalsC,index);
        free(C);
  //    }

    //  #pragma omp section
    //  {
         D =malloc(sizeof(double)*md*nd);
        copy_A_to_B(M,m,n,D,md,nd,ma,na);
        singvalsD = SVD_1(D,  md,  nd);
        PrintVec(singvalsD,index);
        free(D);
      //}
//  }
  concatenate(singvals,singvalsA,singvalsB,ma,mb,index);
  concatenate(singvals,singvalsC,singvalsD,mc,md,index);
  PrintVec(singvals,index);
  return singvals;
}

double* SVD_1_Blocks(double *A, int m, int n) {
    double *M;  // will store A^TA or AA^T
    int nB;    // B dimension
    int na,ma,nb,mb,nc,mc,nd,md;
    double* A1;
    double* B;
    double* C;
    double* D;
    PrintVec(A,m*n);

    // Compute A^TA or AA^T depending on the comparaison between n and m
    if (n < m) {
        M = get_AtA(A, m, n);
        PrintVec(M,n);
        nB = n;
    } else {
        M = get_AAt(A, m, n);
        PrintVec(M,m);
        nB = m;
    }
    if (nB%2==0){
      ma = nB/2;
      mb=ma;
    }else{
      ma = (nB+1)/2;
      mb = (nB-1)/2;
    }
    if (nB%2==0){
      na = nB/2;
      nb=na;
    }else{
      na = (nB+1)/2;
      nb = (nB-1)/2;
    }
    mc=ma;
    nc=na;
    md=mb;
    nd=nb;
    double *singvalsA;
    double *singvalsB;
    double *singvalsC;
    double *singvalsD;
    double *singvals;

    singvals = malloc(sizeof(double)*m*4);
  /* singvalsA = malloc(sizeof(double)*ma);
   singvalsB = malloc(sizeof(double)*mb);
   singvalsC = malloc(sizeof(double)*mc);
   singvalsD = malloc(sizeof(double)*md); */

    int index = 0;

 //   omp_set_nested(1);
 //  #pragma omp parallel sections num_threads(4) private(A,B,C,D)
   //{
     //  #pragma omp section
     //  {
         A1 =malloc(sizeof(double)*ma*na);
         copy_A_to_B(M,nB,nB,A1,ma,na,0,0);
         eigen_t BeigensA = QR_method_Tridiag(A1, na);
         PrintVec(BeigensA.values,na);
         free(A1);
     //  }

     //  #pragma omp section
     //  {
         B =malloc(sizeof(double)*mb*nb);
         copy_A_to_B(M,m,n,B,mb,nb,0,na);
         eigen_t BeigensB = QR_method_Tridiag(B, nb);
         PrintVec(BeigensB.values,nb);
      //   free(B);
     //  }


     //  #pragma omp section
     //  {
  /*       C =malloc(sizeof(double)*mc*nc);
         copy_A_to_B(M,m,n,C,mc,nc,ma,0);
          eigen_t BeigensC = QR_method_Tridiag(C, nc);
         PrintVec(BeigensC.values,nc);
       free(C);*/
   //    }

  //    #pragma omp section
     //  {
  /*        D =malloc(sizeof(double)*md*nd);
         copy_A_to_B(M,m,n,D,md,nd,ma,na);
         eigen_t BeigensD = QR_method_Tridiag(D, nd);
         PrintVec(BeigensD.values,nd);
         free(D); */
       //}
 //  }
   concatenate(singvals,BeigensA.values,BeigensB.values,na,nb,index);
//   concatenate(singvals,BeigensC.values,BeigensD.values,nc,nd,index);
   PrintVec(singvals,index);





//    free(M);
//    free(Beigens.vectors);

    // Beigen values are A singular values
    return singvals;
}
