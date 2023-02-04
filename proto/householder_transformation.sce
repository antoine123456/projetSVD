//Auteur: Jean-francois Brandon
// contact :brandon.jean-francois@ens-paris-saclay.fr
// Projet SVD
// Householder transformation


// From a matrix A we tridiagonalize

A = [5.4 , -4.0, 7.7;
    3.5, -0.7, 2.8;
    -3.2, 5.1, 0.8 ]
function [Ak] = my_householder(A)
    alpha = -sign(A(2,1))*norm(A)
    n=1
    m= size(A)(1)
    disp("taille de A ",m)
    r = sqrt(0.5*(alpha**2-A(2,1)*alpha))
    disp("r=",r)
    v = zeros(size(A)(1),1)
    v(2)= 0.5*(A(2,1)-alpha)/r
    for i = 3:size(v)(1)
        v(i) = 0.5*A(i,1)/r
    end
    P = eye(m,m) - 2*v*(v')
    disp(P)
    Ak= P*A
    disp(Ak)
    Ak= Ak*P
    disp(Ak)
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
        Ak= Ak*P
        disp(Ak)
    end
endfunction

function [Ak] = other_householder(A)
    n = size(A)(1)
    x =rand(n,1)
    alpha = -sign(A(1,1))*norm(A(:,1),2)
    v = A(:,1) - alpha*[1;zeros(n-1,1)]
    beta = 2/ (norm(v,2)**2)
    P = eye(n,n) -beta*v*(v')
    Ak = P*A
    Ak= Ak*P
    for k = 2:n-1
        alpha = -sign(Ak(k,k))*norm(Ak(k:n,k),2)
        disp(Ak(k:n,k))
        v(1:k-1)=zeros(k-1,1)
        disp(v)
        v(k:n) = Ak(k:n,k) - alpha*[1;zeros(n-k,1)]
        disp(v)
        beta = 2/ (norm(v,2)**2)
        P = eye(n,n) - beta*v*(v')
        disp(P)
        Ak = P*A
        disp(Ak)
        Ak= Ak*P
        disp(Ak)
    end
endfunction

function [Q,H] = tri_householder(A)
    [n,m] = size(A)
    R = A
    Q = eye(m,m)
    for k = 1:n-1
        v(k:n) = R(k:n,k)
        normvk = norm(v(k:n))
        v(1) = v(1) + sign(v(1))*normvk
        v(k:n) = (1/normvk)*v(k:n)
        disp(R(k:n,k:n))
        R(k:n,k:n) = R(k:n,k:n) - 2*v(k:n)*(v(k:n)'*R(k:n,k:n))
        disp(R)
        Q(k:n) = Q(k:n) - 2*v(k:n)*(v(k:n)'*Q(k:n))
        disp(Q)
    end
endfunction

function [H,Q] = Hessenberg(A)
    [n,m] = size(A)
    H = A
    Q = eye(m,m)
    for k = 1:n-2
        v = H(k+1:n,k)
        v(1) = v(1) + sign(v(1))*norm(v,2)
        v =(1/norm(v,2))*v
        H(k+1:n,k:n) = H(k+1:n,k) - 2*v(k:n)*(v(k:n)'*H(k+1:n,k))
        H(:,k+1:n) = H(:,k+1:n) - 2*(H(:,k+1:n)*v(k+1:n))*(v(k+1:n)')
        Q(k+1:n) = Q(k+1:n) - -2*v(k+1)*((v(k+1:n))'*Q(k+1:n))
    end
endfunction

//For any sqare matrix A, we can find

function [H,Q] = hessred(A)
    [n,m] = size(A);
    Q = eye(n,n)
    H = A

    for j = 1:n-2
        // -- Find W = I-2vv’ to put zeros below H(j+1,j)
        u =  H(j+1:n,j)
        disp("u=H(j+1:n,j)")
        disp(u)
        disp("norme de u:")
        disp(norm(u))
        u(1) = u(1) + sign(u(1))* norm(u)
        disp("u modified :")
        disp(u)
        disp("norme de u")
        disp(norm(u))
        v =  u/norm(u)
        disp("v = u/norm(u) : ")
        disp(v)
        // --Find W = I - 2vv' to put zeros below H(j+1,j)
        H(j+1:n,:) = H(j+1:n,:)-2*v*(v'*H(j+1:n,:))
        disp("vt*H(j+1:n,:)")
        disp(v'*H(j+1:n,:))
        disp("v*vt*H(j+1:n,:)")
        disp(v*v'*H(j+1:n,:))
        disp("H(j+1:n,:)")
        disp(H(j+1:n,:))
        H(:,j+1:n) = H(:,j+1:n)-(H(:,j+1:n)*(2*v))* v'
        disp("H(:,j+1:n)")
        disp(H(:,j+1:n))
        Q(:,j+1:n) = Q(:,j+1:n)-(Q(:,j+1:n)*(2*v))* v'
        disp("Q(:,j+1:n)")
        disp(Q(:,j+1:n))
    end
endfunction

function [H,Q] = trans_householder(A)
    [n,m] = size(A)
    for i= 1:n-2
    end
endfunction
