//Auteur: Jean-francois Brandon
// contact :brandon.jean-francois@ens-paris-saclay.fr
// Projet SVD
// Householder transformation


// From a matrix A we tridiagonalize

function [Ak] = my_householder(A)
    alpha = -sign(A(2,1))*norm(A)
    n=1
    m= size(A)(1)
    r = sqrt(0.5*(alpha**2-A(2,1)*alpha))
    v = zeros(size(A)(1),1)
    v(2)= 0.5*(A(2,1)-alpha)/r
    for i = 3:size(v)(1)
        v(i) = 0.5*A(i,1)/r
    end
    P = eye(m,m) - 2*v*(v')
    disp(P)
    Ak= P*A
    disp(P*A)
    for k = 2:m-1
        alpha = -sign(Ak(k+1,k))*norm(Ak)
        r = sqrt(0.5*(alpha**2-A(k+1,k)*alpha))
        v = zeros(size(A)(1),1)
        v(k+1)= 0.5*(A(k,k)-alpha)/r
        for i = k+2:m
            v(i) = 0.5*A(i,k)/r
        end
        P = eye(m,m) - 2*v*(v')
        disp(P)
        Ak= P*Ak
        disp(Ak)
    end    
endfunction
