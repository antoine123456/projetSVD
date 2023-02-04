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
#include"io_matrix.h"

    // Matrix A,P,H as : A = P*H*P'
double A[] = {5.4 , -4.0, 7.7,
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



static void Hess_Reduction_test_Q_mat(void **state){
    (void) state; /*unused*/

    double *verifyM1;
    double *H_test;
    double *P_test;
    double res1;
    verifyM1 = calloc(9,sizeof(double*));

    H_test = calloc(9,sizeof(double*));

    P_test = calloc(9,sizeof(double*));


      Hess_Reduction(A, 3, H_test,P_test);


      // double myzeros = 0.0;

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


      Hess_Reduction(A, 3, H_test,P_test);

      // double myzeros = 0.0;

      C_aAplusbB( 3,verifyM2,1.0,H_sol, -1.0,H_test);
      res2= norme_Mat(3, verifyM2);

       assert_true((res2 < 1e-8)==0);
}
static void SVD_Hessenberg_test(void **state){
  double *singvals;
  singvals = SVD_Hessenberg(M,2 ,3 ) ;
  // printf("%lf , %lf\n", sqrt(singvals[0]),sqrt(singvals[1]));
  double res1= sqrt(singvals[0]) -5;
  double res2= sqrt(singvals[1]) -3;
  assert_true((res1< 0.00001));
  assert_true((res2< 0.00001));

}
int main(int argc, char* argv[]){
  double *singvals;
  singvals = SVD_Hessenberg(M,2 ,3 ) ;
  // PrintVec(singvals,6);

  const struct CMUnitTest tests[] = {
      cmocka_unit_test(Hess_Reduction_test_Q_mat),cmocka_unit_test(Hess_Reduction_test_H_mat),
      cmocka_unit_test(SVD_Hessenberg_test)
  };

  return cmocka_run_group_tests(tests, NULL, NULL);
}
