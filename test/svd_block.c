#include<stdio.h>
#include<stdlib.h>
#include"SVD.h"
#include"hessenberg.h"
#include "test_mat.h"
#include "test_utils.h"
#include <io_matrix.h>
#include <gen_mat.h>
#include <lapacke_utils.h>


//#include<cblas.h>
/*
int search_mult(int n, int* tab_n){
    int ref_mult[] = {97,89,83,79,73,71,67,61,59,53,47,43,41,37,31,29,23,19,17,13,11,7,5,3,2};
    int mult;
    for(int i = 0 ;i<25; i++){
        if(n%ref_mult[i] == 0){
           tab_n = (int*) malloc(ref_mult[i]*sizeof(int*));
           int c2 = 0;
           do{
            tab_n[c2] = n/ref_mult[i];
            c2++;
            }while(c2<ref_mult[i]);
            for(int k = 0; k<c2; k++){
                printf("n[%d]=%d\n",k,tab_n[k]);
            }
          mult = ref_mult[i];
          return mult;
        }else{
        int tab_n[] = {(n+1)/2, (n-1)/2};
        int c1=2;
        for(int k = 0; k<c1; k++){
        printf("n[%d]=%d\n",k,tab_n[k]);
        }
        mult = 0;
    }
    if(n%2==1){
        int tab_n[] = {(n+1)/2, (n-1)/2};
        int c1=2;
        for(int k = 0; k<c1; k++){
        printf("n[%d]=%d\n",k,tab_n[k]);
        }
        mult = 0;
    }
    return mult;
}
*/




typedef struct matrix_Block_S {
    double **matrix_blocks;
    int *t_m;
    int *t_n;
    int subd_r;
    int subd_c;
} mat_blk_t;

/* Display in terminal the matrix "mat" with "n" ligns and "m" cols */

void decomp_dim(int n,int* t_n){
    //int p = n/2;
    if(n%2==0){
        t_n[0] = n/2;
        t_n[1] = n/2;
        printf("n1=%d , n2= %d\n",t_n[0],t_n[1]);
    }else{
        t_n[0] = (n+1)/2;
        t_n[1] = (n-1)/2;
        printf("n1=%d , n2= %d\n",t_n[0],t_n[1]);
    }
}

void decomp_dim3(int n,int* t_n){
    //assert(n>8);
    if(n%3==0){
        t_n[0] = n/3;
        t_n[1] = n/3;
        t_n[2] = n/3;
    //    printf("n1=%d , n2= %d\n",t_n[0],t_n[1],t_n[2]);
    }else if(n%3==1){
        t_n[0] = (n+2)/3;
        t_n[1] = (n-1)/3;
        t_n[2] = (n-1)/3;
    //    printf("n1=%d , n2= %d n3=%d\n",t_n[0],t_n[1],t_n[2]);
    }else if(n%3==2){
        t_n[0] = (n+1)/3;
        t_n[1] = (n+1)/3;
        t_n[2] = (n-2)/3;
    //    printf("n1=%d , n2= %d n3=%d\n",t_n[0],t_n[1],t_n[2]);
    }
}
void init_partition_matrix(int m,int n,mat_blk_t *partition){
    int tab_n[3];
    int tab_m[3];
    printf("m=%d,n=%d\n",m,n);
    decomp_dim3(m,tab_m);
    decomp_dim3(n,tab_n);
    partition->t_m = (int*) malloc(4*sizeof(int*));
    partition->t_n = (int*) malloc(4*sizeof(int*));
    int count = 0;
    for(int i =0; i<3;i++){
        for(int j = 0; j<3; j++){
            partition->t_m[count] = tab_m[i];
            printf("partition->t_m[%d] = tab_m[%d]=%d\n", count, i, tab_m[i]);
            partition->t_n[count] = tab_n[j];
            printf("partition->t_n[%d] = tab_n[%d]=%d\n", count, j, tab_n[j]);
            count++;
        }
    }
    partition->matrix_blocks = (double**) malloc(9*sizeof(double**));
    for (int i = 0; i<9; i++){
        partition->matrix_blocks[i] = (double*) malloc(partition->t_n[i]* partition->t_m[i] * sizeof(double*));
    }
}

/* Display in terminal the matrix "mat" with "n" ligns and "m" cols */
/*void PrintMat(double* mat, int n, int m) {
    for (int i=0 ; i<n ; i++) {
        for (int j=0 ; j<m ; j++)
            printf("%lf ", mat[i*m + j]);
        printf("\n");
    }
}*/
void partition_matrix(int m, int n, double* mat,mat_blk_t *partition){
    int tab_n[3];
    int tab_m[3];
    decomp_dim3(m,tab_m);
    decomp_dim3(n,tab_n);
   // init_partition_matrix( m, n,partition);
    partition->t_m = (int*) malloc(9*sizeof(int*));
    partition->t_n = (int*) malloc(9*sizeof(int*));

    int count = 0;
    for(int i =0; i<3;i++){
        for(int j = 0; j<3; j++){
            partition->t_m[count] = tab_m[i];
            partition->t_n[count] = tab_n[j];
            count++;
        }
    }
    //m[0],m[0],m[0],m[1]
    //n[0],n[1],n[2],n[0]
    //partition->t_m[] = {tab_m[0],tab_m[0],tab_m[1],tab_m[1]} ;
    //partition->t_n[] = {tab_n[0],tab_n[1],tab_n[0],tab_n[1]};
    int start_row =0;
    int end_row=0;
    int start_col=0;
    int end_col =0;
    int r[]={0,0,0,1,1,1,2,2,2};
    int c[]={0,1,2,0,1,2,0,1,2};
     count = 0;
   // printf("matrice entiere\n");
   // printf("hello");
   // PrintMat(mat,m,n);
    partition->matrix_blocks = (double**) calloc(9,sizeof(double**));
    for(int i =0; i<3;i++){
      end_row = start_row+ partition->t_m[count];
        for(int j = 0; j<3; j++){
          if(j==0) {
            start_col = 0;
            end_col = 0;
          }
          assert(partition->matrix_blocks );
          assert(partition->matrix_blocks[i*3+j]=!NULL);
          //printf("sous matrice m = %d, n = %d\n",partition->t_m[i],partition->t_n[j]);
          partition->matrix_blocks[i*3+j] = (double*) calloc(2*partition->t_n[i]* partition->t_m[j] , sizeof(double*));




          end_col =  start_col + partition->t_n[i*3+j];
        //  printf("sous matrice m = %d, n = %d\n",partition->t_m[i*3+j],partition->t_n[i*3+j]);
      //    printf("start_row =%d,end_row=%d\n",start_row,end_row);
        //  printf("start_col =%d,end_col=%d\n",start_col,end_col);
          for( int k = start_row ;k<end_row ;k++){
            for(int l = start_col; l<end_col; l++){
          //      printf("i*3+j= %d ; (k = %d ,l=%d),offset local=%d ,k*n+l= %d\n",i*3+j,k,l,( k - start_row ) * partition->t_n[i*3+j] + l-start_col, k*n+l);
                partition->matrix_blocks[i*3+j][( k - start_row ) * partition->t_n[j] + l-start_col] = mat[k*n+l];
            }
          }
          start_col = start_col+ partition->t_n[i*3+j];
        //  PrintMat(partition->matrix_blocks[i*3+j],partition->t_m[i*3+j],partition->t_n[i*3+j]);
          count = 3*i+j;
        }
        start_row =start_row + partition->t_m[count];

      }
       // printf("sous matrice m = %d, n = %d\n",partition->t_m[i],partition->t_n[i]);
       // PrintMat(partition->matrix_blocks[i],partition->t_m[i],partition->t_n[i]);
}


void tri_selection_decroissant(double tableau[], int n) {
    int i, j, indice_max;
    double temp;
    for (i = 0; i < n - 1; i++) {
        indice_max = i;
        for (j = i + 1; j < n; j++) {
            if (tableau[j] > tableau[indice_max]) {
                indice_max = j;
            }
        }
        if (indice_max != i) {
            temp = tableau[i];
            tableau[i] = tableau[indice_max];
            tableau[indice_max] = temp;
        }
    }
}




int main()
{
    int n = 20;
    int m =40;
    double* mat = (double*) calloc(2*m*n,sizeof(double*));
    int n1=8;
    int m1=10;
    double* mat1 = (double*) malloc(m*n*sizeof(double*));
    double* Bm = (double*) calloc(m*m,sizeof(double*));
    double* Bn = (double*) calloc(n*n,sizeof(double*));
    mat_blk_t submat ={NULL,NULL,NULL, 0,0};
    double *singval;// = (double*) malloc(sizeof(double*)*n);

    GenRandMat2(m,n,mat);
  //  GenRandMat2(m1,n1,mat1);
    //get_symB(mat,m,n,Bn,&n);
  //  printf("Hello World Bn\n");
  //  PrintMat(Bn,n,n);
    //get_symB(mat1,m1,n1,Bm,&m1);
  //  printf("Hello World Bm\n");
    //PrintMat(Bm,m1,n1);

    printf("Hello World mat entier\n");
    PrintMat(mat,m,n);
    singval= SVD_Hessenberg(mat,m,n);
    //ingval= SVD_Hessenberg(mat,m,n);
    printf("singval entier\n");
    PrintMat(singval,n,1);
  //   double norm_mat = cblas_dnrm2(m*n,mat,1);
  //  cblas_dscal(m*n,(double) 1.0/norm_mat,mat,1);
  //  PrintMat(mat,m,n);
  //  printf("M entier * 1/norm\n");
  /*  double norm_mat = cblas_dnrm2(m*n,mat,1);




    singval= SVD_Hessenberg(mat,m,n);
    printf("singval entier /norm de mat\n");
    PrintMat(singval,n,1); */
      //PrintMat(mat,m,n);

    double **singvals = (double**) malloc(sizeof(double**)*9);
    double **Wi= (double**) malloc(sizeof(double**)*9);
    partition_matrix( m, n, mat,&submat);
    int place=0;
    int sumlength=0;
    for(int i =0; i<9; i++){
      sumlength = sumlength + MIN(submat.t_m[i],submat.t_n[i]);
    }
    double *result = (double*) calloc(sumlength,sizeof(double*));
    for(int i =0 ; i< 9;i++){
        singvals[i] = (double*) malloc(sizeof(double*)*submat.t_n[i]);
        Wi[i] = (double*) malloc(sizeof(double*)*submat.t_n[i]);
    //   printf("sous matrice m = %d, n = %d\n",submat.t_m[i],submat.t_n[i]);
       //PrintMat(submat.matrix_blocks[i],submat.t_m[i],submat.t_n[i]);
      //  Wi[i]= SVD_Hess_blocks(submat.matrix_blocks[i],submat.t_m[i],submat.t_n[i],m,n);
        singvals[i]= SVD_Hess_blocks(submat.matrix_blocks[i],submat.t_m[i],submat.t_n[i],m,n);
      //  printf("singvals\n");
      //  PrintMat(singvals[i],1,submat.t_n[i]);
      //  cblas_dscal(submat.t_n[i],(double) norm_mat,singvals[i],1);
        cblas_daxpy(submat.t_n[i],1.0,singvals[i],1,&result[place],1);
      //  printf("resultats W \n");
      //  PrintMat(Wi[i],submat.t_n[i],submat.t_n[i]);

        //free(singvals[i]);
        //free(submat.matrix_blocks[i]);
        place = place + MIN(submat.t_m[i],submat.t_n[i]);

    }
  //  result[0] = Wi[0][0] + Wi[1][0] + Wi[2][0];
  //  result[1] = Wi[3][0] + Wi[4][0] + Wi[5][0];
  //  result[2] = Wi[6][0] + Wi[7][0] + Wi[8][0];
  //  printf("resultat = %lf, %lf , %lf\n",result[0],result[1],result[2]);
    tri_selection_decroissant(result, sumlength);
    PrintMat(result,1,MIN(m,n));

    free(singvals);
    free(mat);
    free(Bm);
    free(Bn);
    free(submat.matrix_blocks);
    free(submat.t_n);
    free(submat.t_m);
    free(mat1);
    free(result);
    //decomp_dim( m, n);
    return 0;
}
