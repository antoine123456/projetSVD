#include <SVD.h>
#include "hessenberg.h"
#include <io_matrix.h>
#include <string.h>
#include<omp.h>
#include<stdlib.h>
#include<math_mat.h>
#include<gen_mat.h>
#include<io_matrix.h>
#include<QR.h>
#include <gram_schmidt.h>


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

  //  printf("Matrice rectangulaire A\n");
  //  PrintMat(A,m,n);

    double *B = (double*) calloc(dimB*dimB,sizeof(double*));
    my_test_malloc(B);
  //  printf("Matrice symmetrique B\n");

    if (n < m) {
        B = get_AtA(A, m, n);
        dimB = n;
    } else {
        B = get_AAt(A, m, n);
        dimB = m;
    }

    //get_symB(A,m,n,B,&dimB);

    //PrintMat(B,dimB,dimB);

  //  printf("forme d'Hessenberg de B\n");
    Hess_Reduction(B, dimB);
  //  PrintMat(B,dimB,dimB);

    eigen_t Beigens = QR_method_Tridiag(B, dimB);
    for(int i=0; i<dimB; i++){
      if(Beigens.values[i]>0) Beigens.values[i] = sqrt(Beigens.values[i]);
    }
    free(B);

    free(Beigens.vectors);

    // Beigen values are A singular values
    return Beigens.values;
}

//Calcul la SVD en block poour une sous matrice  Ai qui compose une matrice A
double* SVD_Hess_blocks(double *A, int m, int n, int mg,int ng) {
    assert(A);
    assert(m>1);
    assert(n>1);

  //  double *H ; // hessenberg form of B as B = P'HP
    int dimB;    // B dimension
    //if(choice ==1){
    if(ng<mg){
      dimB = n;
    }else {
      dimB = m;
    }
  //}else{
  //  dimB=n;
  //}
  //  printf("Matrice rectangulaire A\n");
  //  PrintMat(A,m,n);
  //double *B ;  // will store Ai^TAi or AiAi^T
    double *B = (double*) calloc(dimB*dimB,sizeof(double*));
  // Calcule de  Bi = Ai^TAi or AiAi^T
    my_test_malloc(B);
   //printf("Matrice symmetrique B\n");
  //if(choice ==1){
    if (ng < mg) {
        B = get_AtA(A, m, n);
        dimB = n;
    } else {
        B = get_AAt(A, m, n);
        dimB = m;
    }
  //}else{
  //  B = get_AtA(A, m, n);
  //  dimB = n;
  //}
    //get_symB(A,m,n,B,&dimB);

    //PrintMat(B,dimB,dimB);

   //printf("forme d'Hessenberg de B\n");
    Hess_Reduction(B, dimB); //Reduction sous la forme d'Hessenberg de Bi
  // PrintMat(B,dimB,dimB);

    eigen_t Beigens = QR_method_Tridiag(B, dimB); //Decomostion en factorisation QR
  //  printf("Eigens values la sous matrice Bo\n ");
  //  PrintMat(Beigens.values,1,dimB);
  //  printf(" Est modifié Bo\n ");
  //  PrintMat(B,dimB,dimB);
    // les éléments diagonaux de Ri sont les valeurs propres de Bi
    for(int i=0; i<dimB; i++){
     if(Beigens.values[i]>0) Beigens.values[i] = sqrt(Beigens.values[i]);
    }
  //  printf("SIngular values la sous matrice Ai\n ");
  //  PrintMat(Beigens.values,dimB,dimB);
  //  printf("Eigens vectors of B de la sous matrice Ai\n ");
  //  PrintMat(Beigens.vectors,dimB,dimB);
    //Les valeurs propres de Bi sont le carrées de valeurs singulieres de Ai
    // Le Q de Ai est le V ou le U selon comment on a calculer Bi
    // Ui si Bi=Ai*Ai^t  sont les vecteur propres de B donc Beigens.vectors
    // Vi si Bi= Ai^t*Ai sont les vecteur propres de B donc Beigens.vectors
//    if (ng < mg) {
//        B = get_AtA(A, m, n);
//        dimB = n;
//    } else {
//        B = get_AAt(A, m, n);
//        dimB = m;
//    }
    B = get_AAt(A, m, n);
    dimB = m;
  //  printf(" Rettour de B original\n ");
    //PrintMat(B,dimB,dimB);
    // Ui fait partie de la concatenation des autres Ui des matrices Ai pour X
    // Vi fait partie de la concatenation des autres Vi des matrices Ai pour Y
    QR_t QR_B;
    QR_t QR_X ;
    QR_t QR_Y;
    QR_X.Q = (double *) calloc(dimB*dimB, sizeof(double ));
    QR_X.R = (double *) calloc(dimB*dimB, sizeof(double ));

    QR_B.Q = (double *) calloc(dimB*dimB, sizeof(double ));
    QR_B.R = (double *) calloc(dimB*dimB, sizeof(double ));
    GramSchmidt(&QR_B,B,dimB); //On calcule soit Rx ou Ry
  //  GramSchmidt(&QR_XY,Beigens.vectors,dimB); //On calcule soit Rx ou Ry
    //printf(" R pour  B\n ");
    //PrintMat(QR_B.R,dimB,dimB);
    //printf("  Q PourB\n ");
    //PrintMat(QR_B.Q,dimB,dimB);
    qr_decomposition(QR_B.Q, dimB, QR_X.Q, QR_X.R);
    //GramSchmidt(&QR_XY,QR_B.Q,dimB); //On calcule soit Rx ou Ry
    //printf("Matrice Rx du U de la sous matrice Ai \n ");
    //PrintMat(QR_X.R,dimB,dimB);
    //printf("Matrice Qx U ou de la sous matrice Ai \n ");
    //PrintMat(QR_X.Q,dimB,dimB);

    QR_Y.Q = (double *) calloc(dimB*dimB, sizeof(double ));
    QR_Y.R = (double *) calloc(dimB*dimB, sizeof(double ));
    QR_B.Q = (double *) calloc(dimB*dimB, sizeof(double ));
    QR_B.R = (double *) calloc(dimB*dimB, sizeof(double ));
    B = get_AtA(A, m, n);
    dimB = n;
    GramSchmidt(&QR_B,B,dimB); //On calcule soit Rx ou Ry
  //  GramSchmidt(&QR_XY,Beigens.vectors,dimB); //On calcule soit Rx ou Ry
  //  printf(" R pour  B\n ");
  //  PrintMat(QR_B.R,dimB,dimB);
  //  printf("  Q PourB\n ");
  //  PrintMat(QR_B.Q,dimB,dimB);
    qr_decomposition(QR_B.Q, dimB, QR_Y.Q, QR_Y.R);
    //GramSchmidt(&QR_XY,QR_B.Q,dimB); //On calcule soit Rx ou Ry
  //  printf("Matrice Ry du V de la sous matrice Ai \n ");
  //  PrintMat(QR_Y.R,dimB,dimB);
  //  printf("Matrice Qy du V  de la sous matrice Ai \n ");
  //  PrintMat(QR_Y.Q,dimB,dimB);


    double *W = (double*) calloc(dimB*dimB,sizeof(double*));
    for(int i = 0 ; i< dimB; i ++){
      cblas_dger(CblasRowMajor,dimB,dimB,Beigens.values[i],&QR_X.R[i*dimB],1,&QR_Y.R[i*dimB],1,W,dimB);
    }

  //  printf("W = RxS²Rx^t ou RyS²Ry^t\n ");
  //  PrintMat(W,dimB,dimB);
    /*W = get_AAt(W, dimB, dimB);
    Hess_Reduction(W, dimB);
    printf("W reduit\n ");
    PrintMat(W,dimB,dimB);*/

    //Beigens.values = SVD_Hessenberg(W, dimB,dimB); //Decomostion en factorisation QR
  //  for(int i=0; i<dimB; i++){
  //     if(Beigens.values[i]>0) Beigens.values[i] = sqrt(Beigens.values[i]);
      // Beigens.values[i] = sqrt(Beigens.values[i]);
  //    }
  //  printf("Eigens values la sous matrice W\n ");
  //  PrintMat(Beigens.values,1,dimB);

    free(B);
    free(QR_B.Q);
    free(QR_B.R);
    free(QR_X.Q);
    free(QR_X.R);
    free(QR_Y.Q);
    free(QR_Y.R);
    free(W);

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
