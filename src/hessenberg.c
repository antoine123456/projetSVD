#ifndef TRIGO_H
#define TRIGO_H
#include "hessenberg.h"
#include "math_mat.h"
#include "gen_mat.h"
#include<string.h>
#include <errno.h>
#include <assert.h>

/*
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
*/
double signe( double scalar){
  if(scalar<0){
    return -1.0;
  }else{
    return 1.0;
  }
}
double *sub_Copy_vec(double *dest, double *src, int a, int b){
  for(int i = a ; i<b ; i++){
    dest[i] = src[i];
  }
  return dest;
}

void get_symB(double * A, int m, int n, double* B, int *dimb){
    assert(A);
    assert(B);
    assert(m>1);
    assert(n>1);
    assert(dimb);

    if(m>n) {
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, m, 1.0, A, m, A, m, 0.0, B,n);
        if(*dimb != n)*dimb=n;

	}else{
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, m, n, 1.0, A, n, A, n, 0.0, B, m);
        if(*dimb !=m )*dimb=m;
     }
}


void Hess_Reduction(double *A, int n ){
  assert(n>=2);
  assert(A);

  // On itere sur chaque colonne de la matrice A
  for (int k = 0; k < n - 1; k++) {
    // Allocation de mémoire pour les vecteurs x, v et w
    double* x = (double*) calloc(n,sizeof(double*));
    double* v = (double*) calloc(n,sizeof(double*));
    double* w = (double*) calloc(n,sizeof(double*));
    double norm_x = 0;
    double norm_v = 0;

     // Calcul de x en extraisant les éléments de la colonne k du sous-bloc courant
    for (int i = 0; i < n - k - 1; i++) {
      x[i] = A[(k + i + 1) * n + k];
    }

    // Calcul de la norme de x
    for (int i = 0; i < n - k - 1; i++) {
      norm_x += x[i] * x[i];
    }
    norm_x = sqrt(norm_x);

    // Calcul de v en ajoutant la valeur absolue de la norme de x à l'élément x[0]
    v[0] = x[0] + copysign(norm_x, x[0]);


    // Copie des éléments de x dans
    for (int i = 1; i < n - k - 1; i++) {
      v[i] = x[i];
    }

    // Calcul de la norme de v
    for (int i = 0; i < n - k - 1; i++) {
      norm_v += v[i] * v[i];
    }
    norm_v = sqrt(norm_v);

    // Normalisation de v
    for (int i = 0; i < n - k - 1; i++) {
      v[i] = v[i] / norm_v;
    }


   // Calcul de w en multipliant la matrice A(k+1:n,:) par v
   for (int j = k + 1; j < n; j++) {
      for (int i = 0; i < n; i++) {
       w[i] = w[i] + A[j * n + i] * v[j - k - 1];
      }
   }

    //printMat(w,n,1,"w");
     // Mise à jour de la matrice A(k+1:n,:) en soustrayant 2 * w * v
    for (int i = 0; i < n; i++) {
      for (int j = k + 1; j < n; j++) {
        A[j * n + i] = A[j * n + i] - 2 * w[i] * v[j - k - 1];
      }
    }
    //reinitialisation de w
    for(int i = 0; i<n; i++){
      w[i]= 0;
    }
    // w = A(:,k+1:n) * v^t
    for (int i = 0; i < n; i++) {
      for (int j = k + 1; j < n; j++) {
        w[i] = w[i] + A[i * n + j] * v[j - k - 1];
       }
    }
    // Mise à jour de la matrice A(:,k+1:n) en soustrayant 2 * v * w
    for (int i = 0; i < n; i++) {
      for (int j = k + 1; j < n; j++) {
       A[i * n + j] = A[i * n + j] -2 * v[j - k - 1] * w[i];
      }
    }
    free(w);
    free(v);
    free(x);
  }
}

void Hess_Reduction2(double *A,int n,double* U,double* V) {
  assert(n>=2);
  assert(A);

  // On itere sur chaque colonne de la matrice A
  for (int k = 0; k < n - 2; ++k) {
    // Allocation mémoire des vecteurs de calcul intemédiaires
    double* x = (double*) calloc(n,sizeof(double*));
    double* v = (double*) calloc(n,sizeof(double*));
    double *w = (double*) calloc(n,sizeof(double*));

    //Norme des vecteurs intermédiaires
    double norm_v =0.0;
    double norm_x=0.0;

    //Copie du sous-vecteur colonne de x <- A(k+1:n,k) de taille n-k-1
    cblas_dcopy(n - k - 1, &A[(k ) * n + k+1], 1, x, 1);

    norm_x = cblas_dnrm2(n - k - 1, x, 1); //Calcul de la norme 2 euclidienne de x

    x[0] = x[0] + copysign(norm_x, x[0]); //Modification du premier element

    cblas_dcopy(n - k - 1, x, 1, v, 1); //Copy du x modifié dans le v ,vecteur de taille n-k-1 x 1

    norm_v = cblas_dnrm2(n - k - 1, v, 1); //Calcul de la norme 2 euclidienne de v

    cblas_dscal(n - k - 1, 1.0 / norm_v, v, 1); //Normalisation de v

    //Multiplication matriciel : A(k+1:n,:) <-  A(k+1:n,:) - 2*v*v^t* A(k+1:n,:)
    // w <- v^t*A(k+1:n,:) avec w 1 x n
    cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, 1,n,n-k-1,1.0,v,n-k-1,&A[(k+1)*n],n,0,w,n);
    // A(k+1:n,: ) <- v*w
    cblas_dger(CblasRowMajor,n-k-1,n,-2.0,v,1,w,1,&A[(k+1)*n],n);

    // Multiplication matriciel A(:,K+1:n) <- A(:,K+1:n) - 2 *A(:,K+1:n)*v*v^t
    //reinitialisation de w
    for(int i = 0; i<n; i++){
      w[i]= 0;
    }
    // w = A(:,k+1:n) * v^t
    for (int i = 0; i < n; i++) {
       w[i] = w[i] + cblas_ddot(n-k-1,&A[i*n+k+1],1,v,1);
    }
    // Mise a jour par  lignge
    for (int i = 0; i < n; i++) {
      cblas_daxpy(n-k-1,-2*w[i],v,1,&A[i*n+k+1],1);
    }

    free(w);
    free(v);
    free(x);
  }
}
#endif
