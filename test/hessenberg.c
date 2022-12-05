#include<stdio.h>
#include<stdlib.h>
#include<lapacke.h>
#include<stdarg.h>
#include<stddef.h>
#include<setjmp.h>
#include<cmocka.h>
#include"trigonalisation.h"
#include"io_matrix.h"

static void Hess_Reduction_test_positive(void **state){
    (void) state; /*unused*/

    double *H_test;
    H_test = calloc(9,sizeof(double*));
    double *P_test;
    P_test = calloc(9,sizeof(double*));


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

      Hess_Reduction(A, 3, H_test,P_test);


      assert_memory_equal(H_test,H_sol,9*sizeof(double*));

}

int main(int argc, char* argv[]){

  double *H_test;
  H_test = calloc(9,sizeof(double*));
  double *P_test;
  P_test = calloc(9,sizeof(double*));
  double *verifyM1;
  double *verifyM2;
  verifyM1 = calloc(9,sizeof(double*));
  verifyM2 = calloc(9,sizeof(double*));

  // Matrix A,P,H as : A = P*H*P'
    double A[] = {5.4 , -4.0, 7.7,
        3.5, -0.7, 2.8,
        -3.2, 5.1, 0.8} ;

    double H_sol[] = {5.4,     8.1478381 ,   2.9837451,
          -4.7423623,  -3.9512228,   1.5439751,
          4.441D-16,  -0.7560249,   4.0512228};

    double P_sol[] =  {   1.,   0.,          0.,
              0.,  -0.7380288 ,  0.6747692,
              0.,   0.6747692,   0.7380288 };
    Hess_Reduction(A, 3, H_test,P_test);
    printf("H test : \n");
    PrintMat(H_test,3,3);
    printf("P_test : \n");
    PrintMat(P_test,3,3);

    C_aAplusbB( 3,verifyM1,1.0,P_sol, -1.0,P_test);
    C_aAplusbB( 3,verifyM2,1.0,H_sol, -1.0,H_test);

    printf("Hsol-Htest : \n");
    PrintMat(verifyM2,3,3);
    printf("Psol-Ptest : \n");
    PrintMat(verifyM1,3,3);

    double res1= norme_Mat(3, verifyM1);
    double res2= norme_Mat(3, verifyM2);
    printf("\n res1 = %lf, res2 = %lf \n",res1,res2);
  //  PrintMat(H_test,3,3);
  //  PrintMat(P_test,3,3);


  const struct CMUnitTest tests[] = {
      cmocka_unit_test(Hess_Reduction_test_positive),
  };

  return cmocka_run_group_tests(tests, NULL, NULL);
}
