#include<stdio.h>
#include<stdlib.h>
#include<lapacke.h>
#include<stdarg.h>
#include<stddef.h>
#include<setjmp.h>
#include<cmocka.h>
#include"trigonalisation.h"
#include"SVD.h"
#include<math.h>
#include "QR.h"
#include <assert.h>
#include <limits.h>
#include"io_matrix.h"
#include "test_mat.h"
#include "test_utils.h"

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

double M2[] = {11,-5, 4,  6,  6,
              -2, 10, 5,  7,  9,
              3,  -6, 7, 9, 10,
              3,  -10, 6,-6,-5};



static void Hess_Reduction_test_Q_mat(void **state){
    (void) state; /*unused*/

    double *verifyM1;
    double *H_test;
    double *P_test;
    double res1;
    verifyM1 = calloc(9,sizeof(double*));

    H_test = calloc(9,sizeof(double*));

    P_test = calloc(9,sizeof(double*));


      Hess_Reduction(A1, 3, H_test,P_test);


      double myzeros = 0.0;

      C_aAplusbB( 3,verifyM1,1.0,P_sol, -1.0,P_test);
      res1= norme_Mat(3, verifyM1);

       assert_true((res1 < 1e-8)==0);

}

static void Hess_Reduction_test_H_mat(void **state){
    (void) state; /*unused*/

    double *verifyM2;
    double *H_test;
    double *P_test;
    double res2;
    verifyM2 = calloc(9,sizeof(double*));

    H_test = calloc(9,sizeof(double*));

    P_test = calloc(9,sizeof(double*));


      Hess_Reduction(A1, 3, H_test,P_test);

      double myzeros = 0.0;

      C_aAplusbB( 3,verifyM2,1.0,H_sol, -1.0,H_test);
      res2= norme_Mat(3, verifyM2);

       assert_true((res2 < 1e-8)==0);
}
static void SVD_Hessenberg_test(void **state){
  double *singvals;
  singvals = SVD_Hessenberg(M,2 ,3 ) ;
  printf("eigenvalues = 5,3 \n");
  printf("%lf , %lf\n", sqrt(singvals[0]),sqrt(singvals[1]));
  double res1= sqrt(singvals[0]) -5;
  double res2= sqrt(singvals[1]) -3;
  assert_true((res1< 0.00001));
  assert_true((res2< 0.00001));
}

double *singvals;
double* M3 ;
int main(int argc, char* argv[]){



  M3 =  malloc(sizeof(double)*5*6);
  GenRandMat2(5,6, M3);
//  PrintMat(M,2,3);
  singvals = SVD_Hessenberg(M,2 ,3 ) ;
  printf("eigenvalues M\n");
  PrintVec(singvals,2);
//  Eigen_to_Singular(singvals, 4);
  for(int i =0; i<4;i++){
    singvals[i]= sqrt(singvals[i]);
  }
  printf("singular M\n");
  PrintVec(singvals,2);
  //PrintMat(M,2,3);

  singvals = SVD_Hessenberg(M1, 3,4 ) ;
  printf("eigenvalues M1\n");
  PrintVec(singvals,3);
//  Eigen_to_Singular(singvals, 4);
  for(int i =0; i<4;i++){
    singvals[i]= sqrt(singvals[i]);
  }
  printf("singular M1 \n");
  PrintVec(singvals,3);
  PrintMat(M1,3,4);


  singvals = SVD_Hessenberg(M2, 4,5) ;
  printf("eigenvalues M2\n");
  PrintVec(singvals,4);
//  Eigen_to_Singular(singvals, 4);
  for(int i =0; i<4;i++){
    singvals[i]= sqrt(singvals[i]);
  }
  printf("singular M2 \n");
  PrintVec(singvals,4);
    PrintMat(M2,4,5);
  singvals = SVD_Hessenberg(M3, 5,6 ) ;
  printf("eigenvalues M3\n");
  PrintVec(singvals,4);
//  Eigen_to_Singular(singvals, 4);
  for(int i =0; i<4;i++){
    singvals[i]= sqrt(singvals[i]);
  }
  printf("singular M3\n");
  PrintVec(singvals,4);
  PrintMat(M3,3,4);




/*  // Test parameters
  int n = 10;        // matrix order
  double eps = 0.01; // error threshold

  int nerr = 0;
  double *singvals;

  double *givens_mat = givens(n,n);
  double *givens_ev = givens_eigenvalues(n);
  //printf("eigenvalues \n");
  //PrintVec(givens_ev,n);
  Eigen_to_Singular(givens_ev, n);
  //printf("Singular values\n");
  //PrintVec(givens_ev,n);
  singvals = SVD_Hessenberg(givens_mat, n, n);
  //printf("eigen value \n");
  //PrintVec(singvals,n);
  Eigen_to_Singular(givens_ev, n);
  //printf("singular \n");
  //PrintVec(singvals,n);
  if(!Test_Eigenvalues(givens_ev, singvals, eps, n)) {
      Print_Error_Mat("SVD_1", "Givens");
      nerr++;
  }

  double *aegerter_mat = aegerter(n);
  double *aegerter_ev = aegerter_eigenvalues(n);
  //printf("eigenvalues \n");
  //PrintVec(aegerter_ev,n);
  Eigen_to_Singular(givens_ev, n);
  //printf("Singular values\n");
  //PrintVec(aegerter_ev,n);
  singvals = SVD_Hessenberg(aegerter_mat, n, n);
  Eigen_to_Singular(aegerter_ev, n);
  //printf("eigen value \n");
  //PrintVec(singvals,n);
  Eigen_to_Singular(givens_ev, n);
  //printf("singular \n");
  //PrintVec(singvals,n);

  if(!Test_Eigenvalues(aegerter_ev, singvals, eps, n)) {
      Print_Error_Mat("SVD_1", "Aegerter");
      nerr++;
  }

  // Clean memory
  free(givens_mat);
  free(givens_ev);
  free(aegerter_mat);
  free(aegerter_ev);
  free(singvals);

  if (nerr == 0)
      Print_Success("SVD 1");
*/
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(Hess_Reduction_test_Q_mat),cmocka_unit_test(Hess_Reduction_test_H_mat),
      cmocka_unit_test(SVD_Hessenberg_test)
  };

  return 0 ;//cmocka_run_group_tests(tests, NULL, NULL);
}
