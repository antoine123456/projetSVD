#ifndef TRIGO_H
#define TRIGO_H
#include <math.h>
#include <stdlib.h>
#include<stdio.h>
#include <math_mat.h>
#include <lapacke.h>


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

void Hess_Reduction( double *A, int n, double *H, double *Q);
double signe(double scalar);
void sub_Copy(double *dest, double *src,int n, int a,int b );

#endif
