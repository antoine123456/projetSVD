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
    // Matrix A,P,H as : A = P*H*P'
double A1[] = {5.4 , -4.0, 7.7,
    3.5, -0.7, 2.8,
    -3.2, 5.1, 0.8} ;

double H_sol[] = {5.4,     8.1478381 ,   2.9837451,
      -4.7423623,  -3.9512228,   1.5439751,
      4.441e-16,  -0.7560249,   4.0512228};

double P_sol[] =  {   1.,   0.,          0.,
          0.,  -0.7380288 ,  0.6747692,
          0.,   0.6747692,   0.7380288 };

double M[] = {3.0, 2.0, 2.0,
            2.0, 3.0, -2.0};
double U[] = {1.0/sqrt(2), -1.0/sqrt(2),
              1.0/sqrt(2), 1.0/sqrt(2)};
double sigma[] = {5.0, 0.0, 0.0,
                0.0, 3.0, 0.0};
double V[] = { 1.0/sqrt(2.0), -1/(3.0*sqrt(2.0)), -2.0/3.0,
              1.0/sqrt(2.0), 1/(3.0*sqrt(2.0)), 2.0/3.0,
              0.0 ,         2.0*sqrt(2)/3.0,    1.0/3.0};

double M1[] = {11,-5,4,6,
              -2,10,5,7,
              3,-6 ,7, 9};
double M1c[] = {11,-5,4,6,
              -2,10,5,7,
              3,-6 ,7, 9};
double B1sol[] = {198 ,-10, 145
                  -10, 178, 32,
                    145,32,175};

double M2[] = {11,-5, 4,  6,  6,
              -2, 10, 5,  7,  9,
              3,  -6, 7, 9, 10,
              3,  -10, 6,-6,-5};

double B2sol[] = {234,44,205,41,
                  44,259,122,-163,
                  205,122,275,7,
                  41,-163,7,206};

static void get_symB_Test(void **state){

    double *Btest = (double*) calloc(9,sizeof(double*));
    double residual;
    int ldb = 0;

    get_symB(M1,3,4,Btest,&ldb);
   /* for(int i=0; i<9; i++){
        residual = residual + (Btest[i] -B1sol[i])*(Btest[i] -B1sol[i]);
    }
    residual = sqrt(residual);*/
    cblas_daxpy(ldb*ldb,-1.0,Btest,1,B1sol,1);
    residual =cblas_dnrm2(9,B1sol,1);
    assert_true((residual < 1e-6)==0);
}
static void Hessenberg_Reduction_Test(void **state){

    double *Atest = (double*) calloc(9,sizeof(double*));

    double residual;

    cblas_dcopy(9,A1,1,Atest,1);
    Hess_Reduction(A1,3);
    /*for(int i =0; i<9 ; i++){
        residual = residual + (Atest[i] - H_sol[i])*(Atest[i] - H_sol[i]);
    }
    residual = sqrt(residual);*/
    cblas_daxpy(9,-1.0,Atest,1,H_sol,1);
    residual =cblas_dnrm2(9,H_sol,1);
    assert_true((residual < 1e-6)==0);
    free(Atest);
}
static void Hessenberg_Reduction2_Test(void **state){

    double *Atest = (double*) calloc(9,sizeof(double*));
    double residual;

    cblas_dcopy(9,A1,1,Atest,1);
    Hess_Reduction2(A1,3);
    /*for(int i =0; i<9 ; i++){
        residual = residual + (Atest[i] - H_sol[i])*(Atest[i] - H_sol[i]);
    }
    residual = sqrt(residual);*/
    cblas_daxpy(9,-1.0,Atest,1,H_sol,1);
    residual =cblas_dnrm2(9,H_sol,1);
    assert_true((residual < 1e-6)==0);
    free(Atest);
}

static void SVD_Hessenberg_test(void **state){
  double *singvals;
  singvals = SVD_Hessenberg(M ,2,3) ;
//  printf("eigenvalues = 5,3 \n");
//  printf("%lf , %lf\n", sqrt(singvals[0]),sqrt(singvals[1]));
  double res1= sqrt(singvals[0]) -5;
  double res2= sqrt(singvals[1]) -3;
  assert_true((res1< 0.00001));
  assert_true((res2< 0.00001));
}
static void SVD_Hessenberg_test2(void **state){
  double *singvals;
  singvals = SVD_Hessenberg(M1 ,3,4) ;
//  printf("eigenvalues  \n");
  double svsol[] = {18.256839 , 13.530458, 5.883411};
  double res1;

  //printf("%lf , %lf, %lf\n", sqrt(fabs(singvals[0])),sqrt(fabs(singvals[1])),sqrt(fabs(singvals[2])));
  //singvals[0] = sqrt(fabs(singvals[0]));
  //singvals[1]= sqrt(fabs(singvals[1]));
  //singvals[2]= sqrt(fabs(singvals[2]));
  //printf("%lf , %lf, %lf\n", fabs(singvals[0]),fabs(singvals[1]),fabs(singvals[2]));
  cblas_daxpy(3,-1.0,svsol,1,singvals,1);
  //printf("%lf , %lf, %lf\n", fabs(singvals[0]),fabs(singvals[1]),fabs(singvals[2]));
  res1 = cblas_dnrm2(3,singvals,1);
  assert_true((res1< 0.00001));
  //assert_true((res2< 0.00001));
}
/*
static void SVD_Hessenberg_test3(void **state){
  double *singvals;
  double *s = (double*) calloc(3,sizeof(double*));
  double *work = (double*) calloc(3*5,sizeof(double*));
  singvals = SVD_Hessenberg(M2 ,4,3) ;
  int info;
  LAPACK_dgesvd('N','N',4,3,M1c,4,s,NULL,1,NULL,1,work,15,&info);
  printf("eigenvalues  \n");
  printf("%lf , %lf, %lf\n", sqrt(fabs(singvals[0])),sqrt(fabs(singvals[1])),sqrt(fabs(singvals[2])));
  printf("%lf, %lf, %lf,\n",s[0],s[1],s[2]);
  double res1= sqrt(singvals[0]) -s[0];
  double res2= sqrt(singvals[1]) -s[1];
  assert_true((res1< 0.00001));
  assert_true((res2< 0.00001));
}
*/
double *singvals;
double* M3 ;



void SingularValueDecomposition(int m,     // number of rows in matrix
                                int n,     // number of columns in matrix
                                int lda,   // leading dimension of matrix
                                double *a, double *buff) // pointer to top-left corner
{
    // Setup a buffer to hold the singular values:
    int numberOfSingularValues = m < n ? m : n;
    double *s = malloc(numberOfSingularValues * sizeof s[0]);

    // Setup buffers to hold the matrices U and Vt:
    double *u = malloc(m*m * sizeof u[0]);
    double *vt = malloc(n*n * sizeof vt[0]);

    // Workspace and status variables:
    double workSize;
    double *work = &workSize;
    int lwork = -1;
    int *iwork = malloc(8 * numberOfSingularValues * sizeof iwork[0]);
    int info = 0;

    // Call dgesdd_ with lwork = -1 to query optimal workspace size:
    dgesdd_("A", &m, &n, a, &lda, s, u, &m, vt, &n, work, &lwork, iwork, &info);
    if (info) // handle error conditions here

    // Optimal workspace size is returned in work[0].
    lwork = workSize;
    work = malloc(lwork * sizeof work[0]);

    // Call dgesdd_ to do the actual computation:
    dgesdd_("A", &m, &n, a, &lda, s, u, &m, vt, &n, work, &lwork, iwork, &info);
    if (info) // handle error conditions here

    // Cleanup workspace:
    free(work);
    free(iwork);

    // do something useful with U, S, Vt ...
    cblas_dcopy(numberOfSingularValues,s,1,buff,1);
    // and then clean them up too:
    free(s);
    free(u);
    free(vt);
}
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
  for (int i =1 ;i<5;i++){
    double* s_v = (double*) calloc(i*3,sizeof(double*));
    //double* s_v2 = (double*) calloc(i*3,sizeof(double*));
    double* rand_mat = (double*) calloc((i*3)*(i*5),sizeof(double*));
    GenRandMat2(i*3,i*5,rand_mat);
    s_v = SVD_Hessenberg(rand_mat,i*3,i*5);
    //SingularValueDecomposition(i*3,i*5,i*3,rand_mat,s_v2);
    printf("singular values pour n = %d\n",i*3);
    for(int j = 0; j<i*3;j++){
      //s_v[j] = sqrt(fabs(s_v[j]));
      printf("s_v[%d]=%lf\n",j,s_v[j]);
      //printf("s_v2[%d]=%lf\n",j,s_v2[j]);

    }
    free(rand_mat);
  }

  const struct CMUnitTest tests[] = {cmocka_unit_test(Hessenberg_Reduction_Test),cmocka_unit_test(Hessenberg_Reduction2_Test),
                        cmocka_unit_test(get_symB_Test),cmocka_unit_test(SVD_Hessenberg_test),cmocka_unit_test(SVD_Hessenberg_test2)};
  //return EXIT_SUCCESS;
  return cmocka_run_group_tests(tests, NULL, NULL);
}
