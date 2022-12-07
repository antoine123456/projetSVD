#include "trigonalisation.h"
#include "math_mat.h"
#include "gen_mat.h"
#include<string.h>

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

  double norm_u;
  //printf("n=%d \n",n);
  // Q = eye(n,n)
  GenIdentityMatrix(n,Q); // Unitary matrix of Householder
//  PrintMat(Q,n,n);
  // H =A
  memcpy(H,A,n*n*sizeof(double*));  // Hessenberg Matrix
  //PrintMat(H,n,n);
  u = calloc(n,sizeof(double));
  v = calloc(n,sizeof(double));
  save = calloc(n,sizeof(double));
  save2 = calloc(n,sizeof(double));
  save3 = calloc(n,sizeof(double));
  for(int j = 0; j<n-2; j++){
    // Find W  = I - 2vv'
    //vecteurs intermediaire calculée

    //vecteur u = H(j+1:n,j)
    for (int i = j+1; i<n;i++){
      u[i-(j+1)] = H[i*n+j];
      //printf("u[%d] = %lf \n",i-(j+1),u[i-(j+1)]);
      //printf("%lf \n",H[i*n+j]);
    }
    //Calcul norme de u
    norm_u = seq_d_norm(n,u);
    //printf("norm = %lf \n ",norm_u);

    // u(1) = u(1) + signe(u(1))*norm(u)
    u[0] = u[0] + signe(u[0])*norm_u;
    //printf("vecteur u \n");
    //PrintVec(u,n);
    //Calcul v = u/norm(u)
    norm_u = seq_d_norm(n,u);
    seq_d_axpy(n,1/norm_u,u,v);
  //  printf("vecteur v ! \n");
    //PrintVec(v,n);

    // --Find W = I - 2vv' to put zeros below (H(j+1,j))
    //H(j+1:n,:) = H(j+1:n,:)-2*v*(v'*H(j+1:n,:))
    // Compute save = v'*H(j+1:n,:)
    init_vect(n,save);
    for(int k = j+1 ; k<n ; k++){
      for(int l=0; l<n; l++){
        //  printf("%lf=+%lf*lf\n",save[k],v[k],H[k*n+l]);
        save[l] =  save[l]+v[k-(j+1)]*H[k*n+l];
        //printf("%lf=+%lf*lf\n",save[k],v[k],H[k*n+l]);
      }
    }
    //printf("save :\n");
    //PrintVec(save,n);
    //Compute H(j+1:n,:) = H(j+1:n,:)-2*v(1:n-j,1)*save(1,1:n-j)

    for(int p = j+1; p<n ; p++){
      for(int q = 0; q<n ; q++){
        H[p*n+q] = H[p*n+q] - 2.0*v[p-(j+1)]*save[q];
      //  printf("%lf=%lf-%lf,%lf\n",H[p*n+q] ,H[p*n+q] ,v[q],save[q]);
      }
    }
    //printf("H modified: \n");
    //PrintMat(H,n,n);
    //  H(:,j+1:n) = H(:,j+1:n)-(H(:,j+1:n)*(2*v))* v'
    init_vect(n,save2);
    for(int k = j+1 ; k<n ; k++){
      for(int l=0; l<n; l++){
        save2[l] = save2[l]+v[k-(j+1)]*H[l*n+k];
      //  printf("%lf=+%lf*lf\n",save[k],v[k],H[l*n+k]);
      }
    }
    //printf("save2 :\n");
    //PrintVec(save2,n);
    //  H(:,j+1:n) = H(:,j+1:n)-save(1:n-j,1)* v'(1,1:n-j)
    for(int p = j+1; p<n ; p++){
      for(int q = 0; q<n ; q++){
        H[q*n+p] = H[q*n+p] - 2.0*v[p-(j+1)]*save2[q];
      //  printf("%lf=%lf-%lf,%lf\n",H[q*n+p] ,H[q*n+p] ,v[q],save[q]);
      }
    }
    //printf("H modified 2 :\n");
    //PrintMat(H,n,n);
    // Q(:,j+1:n) = Q(:,j+1:n)-(Q(:,j+1:n)*(2*v))* v
    init_vect(n,save3);
    for(int k = j+1 ; k<n ; k++){
      for(int l=0; l<n; l++){
        //printf("k= %d (%d), l =%d\n",k,k-(j+1),l);
        //printf("save[%d] = %lf\n",l,save3[l]);
        save3[l] =save3[l] +v[k-(j+1)]*Q[l*n+k];
        //printf("save[%d]=+v[%d(%d)]*Q[%d,%d]=%lf = %lf*%lf\n",l,k,k-(j+1),l,k,save3[l],v[k-(j+1)],Q[l*n+k]);
      }
    }
    //printf("save3 :\n");
    //PrintVec(save3,n);
    for(int p = j+1; p<n ; p++){
      for(int q = 0; q<n ; q++){
        //printf("p= %d (%d), q =%d\n",p,p-(j+1),q);
        Q[q*n+p] = Q[q*n+p] - 2.0*v[p-(j+1)]*save3[q];
      }
    }
  //  printf("Q : \n");
  //  PrintMat(Q,n,n);
  //  printf("H : \n");
//    PrintMat(H,n,n);

  }
  ///  printf("Q : \n");
  //  PrintMat(Q,n,n);
  //  printf("H : \n");
  //  PrintMat(H,n,n);

  free(u);
  free(v);
  free(save);
  free(save2);
  free(save3);

}
