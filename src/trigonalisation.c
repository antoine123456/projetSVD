#include "trigonalisation.h"
#include "math_mat.h"
#include "gen_mat.h"
#include<string.h>
#include <errno.h>

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
}

// Hessenberg Reduction
// A matrix which will apply the reduction
//Q  Unitary matrix of Householder
//H Hessenberg Matrix
void Hess_Reduction(double *A, int n, double* H, double* Q){
  double *u;  //vector as
  double *v;
  double *save;
  double *save2;
  double *save3;
  int errornum;

  double norm_u;

  // Q = eye(n,n)
  GenIdentityMatrix(n,Q); // Unitary matrix of Householder

  // H =A
  memcpy(H,A,n*n*sizeof(double*));  // Hessenberg Matrix
  //PrintMat(H,n,n);
  u = calloc(n,sizeof(double));
  test_malloc(u);

  v = calloc(n,sizeof(double));
  test_malloc(v);

  save = calloc(n,sizeof(double));
  test_malloc(save);

  save2 = calloc(n,sizeof(double));
  test_malloc(save2);

  save3 = calloc(n,sizeof(double));
  test_malloc(save3);

  for(int j = 0; j<n-2; j++){
    // Find W  = I - 2vv'
    //vecteurs intermediaire calculée

    //vecteur u = H(j+1:n,j)
    for (int i = j+1; i<n;i++){
      u[i-(j+1)] = H[i*n+j];

    }
    //Calcul norme de u
    norm_u = seq_d_norm(n,u);

    // u(1) = u(1) + signe(u(1))*norm(u)
    u[0] = u[0] + signe(u[0])*norm_u;

    //Calcul v = u/norm(u)
    norm_u = seq_d_norm(n,u);
    seq_d_axpy(n,1/norm_u,u,v);

    // --Find W = I - 2vv' to put zeros below (H(j+1,j))
    //H(j+1:n,:) = H(j+1:n,:)-2*v*(v'*H(j+1:n,:))
    // Compute save = v'*H(j+1:n,:)
    init_vect(n,save);
    for(int k = j+1 ; k<n ; k++){
      for(int l=0; l<n; l++){
        save[l] =  save[l]+v[k-(j+1)]*H[k*n+l];
      }
    }

    //Compute H(j+1:n,:) = H(j+1:n,:)-2*v(1:n-j,1)*save(1,1:n-j)
    for(int p = j+1; p<n ; p++){
      for(int q = 0; q<n ; q++){
        H[p*n+q] = H[p*n+q] - 2.0*v[p-(j+1)]*save[q];
      }
    }

    //  H(:,j+1:n) = H(:,j+1:n)-(H(:,j+1:n)*(2*v))* v'
    init_vect(n,save2);
    for(int k = j+1 ; k<n ; k++){
      for(int l=0; l<n; l++){
        save2[l] = save2[l]+v[k-(j+1)]*H[l*n+k];
      }
    }

    //  H(:,j+1:n) = H(:,j+1:n)-save(1:n-j,1)* v'(1,1:n-j)
    for(int p = j+1; p<n ; p++){
      for(int q = 0; q<n ; q++){
        H[q*n+p] = H[q*n+p] - 2.0*v[p-(j+1)]*save2[q];
      //  printf("%lf=%lf-%lf,%lf\n",H[q*n+p] ,H[q*n+p] ,v[q],save[q]);
      }
    }

    // Q(:,j+1:n) = Q(:,j+1:n)-(Q(:,j+1:n)*(2*v))* v
    init_vect(n,save3);
    for(int k = j+1 ; k<n ; k++){
      for(int l=0; l<n; l++){
        save3[l] =save3[l] +v[k-(j+1)]*Q[l*n+k];
      }
    }

    for(int p = j+1; p<n ; p++){
      for(int q = 0; q<n ; q++){
        Q[q*n+p] = Q[q*n+p] - 2.0*v[p-(j+1)]*save3[q];
      }
    }


  }


  free(u);
  free(v);
  free(save);
  free(save2);
  free(save3);

}
