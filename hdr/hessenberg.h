#ifndef TRIGO_H
#define TRIGO_H
#include<stdlib.h>
#include<stdio.h>
#include<cblas.h>
#include<string.h>
#include<lapack.h>
#include<math.h>
#include<assert.h>


/*************
function [H,Q] = hessred(A)
    [n,m] = size(A);
    Q = eye(n,n)
    H = A

    for j = 1:n-2
        // -- Find W = I-2vv’ to put zeros below H(j+1,j)
        u =  H(j+1:n,j)
        u(1) = u(1) + sign(u(1))* norm(u)
        v =  u/norm(u)
        // --Find W = I - 2vv' to put zeros below H(j+1,j)
        H(j+1:n,:) = H(j+1:n,:)-2*v*(v'*H(j+1:n,:))
        H(:,j+1:n) = H(:,j+1:n)-(H(:,j+1:n)*(2*v))* v'
        Q(:,j+1:n) = Q(:,j+1:n)-(Q(:,j+1:n)*(2*v))* v'
    end
endfunction
**************/
void get_symB(double * A, int m, int n, double* B, int *dimb);
void Hess_Reduction( double *A, int n);
void Hess_Reduction2( double *A, int n);
void parallel_Hess_Reduction( double *A, int n);
double signe(double scalar);
double *sub_Copy(double *dest, double *src,int n, int a,int b );

#endif
