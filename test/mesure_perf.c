#include<stdio.h>
#include<stdlib.h>
#include<lapacke.h>
#include<stdarg.h>
#include<stddef.h>
#include<setjmp.h>
#include<cmocka.h>
#include"hessenberg.h"
#include"SVD.h"
#include<math.h>
#include "QR.h"
#include <assert.h>
#include <limits.h>
#include"io_matrix.h"
#include "test_mat.h"
#include "test_utils.h"
#include "io_matrix.h"
#include "gen_mat.h"
#include <clapack.h>
#include<time.h>

int main(int argc, char* argv[]){
/*
  double *singvals;
  double *s = (double*) calloc(3,sizeof(double*));
  double *work = (double*) calloc(3*5,sizeof(double*));
  singvals = SVD_Hessenberg(M2 ,4,3) ;
  int info;
  LAPACKE_dgesvd(LAPACK_ROW_MAJOR,'N','N',4,3,M1c,4,s,NULL,1,NULL,1,work);
  printf("eigenvalues  \n");
  printf("%lf , %lf, %lf\n", sqrt(fabs(singvals[0])),sqrt(fabs(singvals[1])),sqrt(fabs(singvals[2])));
  printf("%lf, %lf, %lf,\n",s[0],s[1],s[2]);

  singvals = SVD_Hessenberg(M1 ,3,4) ;
  printf("eigenvalues  \n");
  printf("%lf , %lf, %lf\n", sqrt(fabs(singvals[0])),sqrt(fabs(singvals[1])),sqrt(fabs(singvals[2])));*/
  char filename[] = "mesure_perf.dat";
  FILE *file_pointer = fopen(filename,"w+");
  clock_t t_init;
  clock_t t_final;
  float t_cpu;

  fprintf(file_pointer,"m\t;n\t;temps\n");
  for (int i =1 ;i<3;i++){
    double* s_v = (double*) calloc(i*3,sizeof(double*));
    //double* s_v2 = (double*) calloc(i*3,sizeof(double*));
    double* rand_mat = (double*) calloc((i*3)*(i*4),sizeof(double*));
    GenRandMat2(i*3,i*4,rand_mat);
    t_init = clock();
    for(int k = 0; k<5; k++){
      s_v = SVD_Hessenberg(rand_mat,i*3,i*4);
      //printf("%lf,%lf,%lf,%lf",s_v[0],s_v[1],s_v[2],s_v[3]);
    }
    t_final = clock();
    t_cpu = 0.001*(t_final - t_init)*1e-6;
    fprintf(file_pointer,"%d;\t%d;\t;%f\n",i*3,i*4,t_cpu);
    //SingularValueDecomposition(i*3,i*5,i*3,rand_mat,s_v2);
    printf("singular values pour n = %d\n",i*3);
    //for(int j = 0; j<i*3;j++){
      //s_v[j] = sqrt(fabs(s_v[j]));
  //    printf("s_v[%d]=%lf\n",j,s_v[j]);
  //    printf("s_v2[%d]=%lf\n",j,s_v2[j]);
      //free(rand_mat);
    //}

  }
fclose(file_pointer);

  return EXIT_SUCCESS;
  //return cmocka_run_group_tests(tests, NULL, NULL);
}
