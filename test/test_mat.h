/* Function isolated from John Burkardt's code base whose purpose is to generate test matrices */

#pragma once

#include <stdlib.h>
#include <math.h>


int i4_max ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.
*/
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}


double *givens ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    GIVENS returns the GIVENS matrix.

  Discussion:

    Note that this is NOT the "Givens rotation matrix".  This
    seems to be more commonly known as the Moler matrix

  Formula:

    A(I,J) = 2 * min ( I, J ) - 1

  Example:

    N = 5

    1 1 1 1 1
    1 3 3 3 3
    1 3 5 5 5
    1 3 5 7 7
    1 3 5 7 9

  Properties:

    A is integral, therefore det ( A ) is integral, and 
    det ( A ) * inverse ( A ) is integral.

    A is positive definite.

    A is symmetric: A' = A.

    Because A is symmetric, it is normal.

    Because A is normal, it is diagonalizable.

    The inverse of A is tridiagonal.

    A has a simple Cholesky factorization.

    A has eigenvalues

      LAMBDA(I) = 0.5 * sec ( ( 2 * I - 1 ) * PI / ( 4 * N ) )^2

    The condition number P(A) is approximately 16 N^2 / PI^2.

    The family of matrices is nested as a function of N.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 October 2010

  Author:

    John Burkardt

  Reference:

    Morris Newman, John Todd,
    Example A9,
    The evaluation of matrix inversion programs,
    Journal of the Society for Industrial and Applied Mathematics,
    Volume 6, Number 4, pages 466-476, 1958.

    John Todd,
    Basic Numerical Mathematics,
    Volume 2: Numerical Algebra,
    Birkhauser, 1980,
    ISBN: 0817608117,
    LC: QA297.T58.

    Joan Westlake,
    A Handbook of Numerical Matrix Inversion and Solution of 
    Linear Equations,
    John Wiley, 1968,
    ISBN13: 978-0471936756,
    LC: QA263.W47.

  Parameters:

    Input, int M, N, the order of the matrix.

    Output, double GIVENS[M*N], the matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = ( double ) ( 2 * i4_min ( i, j ) + 1 );
    }
  }
  return a;
}
/******************************************************************************/

double *givens_eigenvalues ( int n )

/******************************************************************************/
/*
  Purpose:

    GIVENS_EIGENVALUES returns the eigenvalues of the GIVENS matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Output, double GIVENS_EIGENVALUES[N], the eigenvalues.
*/
{
  double angle;
  int i;
  double *lambda;
  const double r8_pi = 3.141592653589793;

  lambda = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    angle = ( double ) ( 2 * i + 1 ) * r8_pi / ( double ) ( 4 * n );
    lambda[i] = 0.5 / pow ( cos ( angle ), 2 );
  }
  return lambda;
}
/******************************************************************************/


double *aegerter ( int n )

/******************************************************************************/
/*
  Purpose:

    AEGERTER returns the AEGERTER matrix.

  Formula:

    if ( I == N )
      A(I,J) = J
    else if ( J == N )
      A(I,J) = I
    else if ( I == J )
      A(I,J) = 1
    else
      A(I,J) = 0

  Example:

    N = 5

    1  0  0  0  1
    0  1  0  0  2
    0  0  1  0  3
    0  0  0  1  4
    1  2  3  4  5

  Properties:

    A is integral, therefore det ( A ) is integral, and 
    det ( A ) * inverse ( A ) is integral.

    A is symmetric: A' = A.

    Because A is symmetric, it is normal.

    Because A is normal, it is diagonalizable.

    A is border-banded.

    det ( A ) = N * ( - 2 * N * N + 3 * N + 5 ) / 6

    A has N-2 eigenvalues equal to 1.

    The other two eigenvalues are

      ( N + 1 + sqrt ( ( N + 1 )^2 - 4 * det ( A ) ) ) / 2
      ( N + 1 - sqrt ( ( N + 1 )^2 - 4 * det ( A ) ) ) / 2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 July 2008

  Author:

    John Burkardt

  Reference:

    MJ Aegerter,
    Construction of a Set of Test Matrices,
    Communications of the ACM,
    Volume 2, Number 8, August 1959, pages 10-12.

  Parameters:

    Input, int N, the order of the matrix.

    Output, double AEGERTER[N*N], the matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( i = 1; i <= n; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      if ( i == n )
      {
        a[i-1+(j-1)*n] = ( double ) ( j );
      }
      else if ( j == n )
      {
        a[i-1+(j-1)*n] = ( double ) ( i );
      }
      else if ( i == j )
      {
        a[i-1+(j-1)*n] = 1.0;
      }
      else
      {
        a[i-1+(j-1)*n] = 0.0;
      }
    }
  }
  return a;
}
/******************************************************************************/


double *aegerter_eigenvalues ( int n )

/******************************************************************************/
/*
  Purpose:

    AEGERTER_EIGENVALUES returns the eigenvalues of the AEGERTER matrix.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Output, double AEGERTER_EIGENVALUES[N], the eigenvalues.
*/
{
  double determ;
  int i;
  double *lambda;
  double np1;

  lambda = ( double * ) malloc ( n * sizeof ( double ) );

  determ = ( double ) ( n - ( ( n - 1 ) * n * ( 2 * n - 1 ) ) / 6 );
  np1 = ( double ) ( n + 1 );

  lambda[0]     = 0.5 * ( np1 - sqrt ( np1 * np1 - 4.0 * determ ) );
  for ( i = 1; i < n - 1; i++ )
  {
    lambda[i] = 1.0;
  }
  lambda[n-1] = 0.5 * ( np1 + sqrt ( np1 * np1 - 4.0 * determ ) );

  return lambda;
}


double *orth_symm ( int n )

/******************************************************************************/
/*
  Purpose:

    ORTH_SYMM returns the ORTH_SYMM matrix.

  Formula:

    A(I,J) = sqrt ( 2 ) * sin ( I * J * pi / ( N + 1 ) ) / sqrt ( N + 1 )

  Example:

    N = 5

    0.326019   0.548529   0.596885   0.455734   0.169891
    0.548529   0.455734  -0.169891  -0.596885  -0.326019
    0.596885  -0.169891  -0.548529   0.326019   0.455734
    0.455734  -0.596885   0.326019   0.169891  -0.548528
    0.169891  -0.326019   0.455734  -0.548528   0.596885

  Properties:

    A is orthogonal: A' * A = A * A' = I.

    A is symmetric: A' = A.

    A is not positive definite (unless N = 1 ).

    Because A is symmetric, it is normal.

    Because A is symmetric, its eigenvalues are real.

    Because A is orthogonal, its eigenvalues have unit norm.

    Only +1 and -1 can be eigenvalues of A.

    Because A is normal, it is diagonalizable.

    A is involutional: A * A = I.

    If N is even, trace ( A ) = 0; if N is odd, trace ( A ) = 1.

    LAMBDA(1:(N+1)/2) = 1; LAMBDA((N+1)/2+1:N) = -1.

    A is the left and right eigenvector matrix for the
    second difference matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2011

  Author:

    John Burkardt

  Reference:

    Morris Newman, John Todd,
    The evaluation of matrix inversion programs,
    Journal of the Society for Industrial and Applied Mathematics,
    Volume 6, Number 4, pages 466-476, 1958.

  Parameters:

    Input, int N, the order of the matrix.

    Output, double ORTH_SYMM[N*N], the matrix.
*/
{
  double *a;
  double angle;
  int i;
  int j;
  const double r8_pi = 3.141592653589793;

  a = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      angle = 2.0 * ( double ) ( ( i + 1 ) * ( j + 1 ) ) * r8_pi 
                  / ( double ) ( 2 * n + 1 );
      a[i+j*n] = 2.0 * sin ( angle ) / sqrt ( ( double ) ( 2 * n + 1 ) );
    }
  }
  return a;
}
/******************************************************************************/

double *orth_symm_eigenvalues ( int n )

/******************************************************************************/
/*
  Purpose:

    ORTH_SYMM_EIGENVALUES returns eigenvalues of the ORTH_SYMM matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Output, double ORTH_SYMM_EIGENVALUES[N], the eigenvalues.
*/
{
  int i;
  double *lambda;

  lambda = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < ( n + 1 ) / 2; i++ )
  {
    lambda[i] = +1.0;
  }
  for ( i = ( n + 1 ) / 2; i < n; i++ )
  {
    lambda[i] = -1.0;
  }
  return lambda;
}
/******************************************************************************/


double *bab ( int n, double alpha, double beta )

/******************************************************************************/
/*
  Purpose:

    BAB returns the BAB matrix.

  Discussion:

    The name is meant to suggest the pattern "B  A  B" formed by
    the nonzero entries in a general row of the matrix.

  Example:

    N = 5
    ALPHA = 5, BETA = 2

    5  2  .  .  .
    2  5  2  .  .
    .  2  5  2  .
    .  .  2  5  2
    .  .  .  2  5

  Properties:

    A is banded, with bandwidth 3.

    A is tridiagonal.

    Because A is tridiagonal, it has property A (bipartite).

    A is Toeplitz: constant along diagonals.

    A is symmetric: A' = A.

    Because A is symmetric, it is normal.

    Because A is normal, it is diagonalizable.

    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).

    The family of matrices is nested as a function of N.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 September 2010

  Author:

    John Burkardt

  Reference:

    CM da Fonseca, J Petronilho,
    Explicit Inverses of Some Tridiagonal Matrices,
    Linear Algebra and Its Applications,
    Volume 325, 2001, pages 7-21.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double ALPHA, BETA, the parameters.

    Output, double BAB[N*N], the matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[i+j*n] = alpha;
      }
      else if ( i == j + 1 || i == j - 1 )
      {
        a[i+j*n] = beta;
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }
  return a;
}
/******************************************************************************/

double *bab_eigenvalues ( int n, double alpha, double beta )

/******************************************************************************/
/*
  Purpose:

    BAB_EIGENVALUES returns the eigenvalues of the BAB matrix.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 September 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Input, double ALPHA, BETA, the parameters.

    Output, double BAB_EIGENVALUES[N], the eigenvalues.
*/
{
  double angle;
  int i;
  double *lambda;
  const double r8_pi = 3.141592653589793;

  lambda = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    angle = ( double ) ( i + 1 ) * r8_pi / ( double ) ( n + 1 );
    lambda[i] = alpha + 2.0 * beta * cos ( angle );
  }

  return lambda;
}
/******************************************************************************/

double *lesp ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    LESP returns the LESP matrix.

  Formula:

    if ( I - J == 1 ) then
      A(I,J) = 1 / I
    else if ( I - J == 0 ) then
      A(I,J) = - ( 2*I+3 )
    else if ( I - J == 1 ) then
      A(I,J) = J
    else
      A(I,J) = 0.0

  Example:

    M = 5, N = 5

     -5    2    .    .     .
     1/2  -7    3    .     .
      .   1/3  -9    4     .
      .    .   1/4 -11     5
      .    .    .   1/5  -13


  Properties:

    The matrix is tridiagonal.

    Because A is tridiagonal, it has property A (bipartite).

    A is generally not symmetric: A' /= A.

    The eigenvalues are real, and smoothly distributed in [-2*N-3.5, -4.5].

    The eigenvalues are sensitive.

    The matrix is similar to the symmetric tridiagonal matrix with
    the same diagonal entries and with off-diagonal entries 1,
    via a similarity transformation using the diagonal matrix
    D = diagonal ( 1, 2, ..., N ).

    The family of matrices is nested as a function of N.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Reference:

    Wim Lenferink, MN Spijker,
    On the use of stability regions in the numerical analysis of initial
    value problems,
    Mathematics of Computation,
    Volume 57, 1991, pages 221-237.

    Lloyd Trefethen,
    Pseudospectra of matrices,
    in Numerical Analysis 1991,
    Proceedings of the 14th Dundee Conference,
    D F Griffiths and G A Watson, editors,
    Pitman Research Notes in Mathematics, volume 260,
    Longman Scientific and Technical, Essex, UK, 1992, pages 234-266.

  Parameters:

    Input, int M, N, the order of the matrix.

    Output, double LESP[M*N], the matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( i - j == 1 )
      {
        a[i+j*m] = 1.0 / ( double ) ( i + 1 );
      }
      else if ( i - j == 0 )
      {
        a[i+j*m] = - ( double ) ( 2 * i + 5 );
      }
      else if ( i - j == -1 )
      {
        a[i+j*m] = ( double ) ( j + 1 );
      }
      else
      {
        a[i+j*m] = 0.0;
      }
    }
  }
  return a;
}
/******************************************************************************/
