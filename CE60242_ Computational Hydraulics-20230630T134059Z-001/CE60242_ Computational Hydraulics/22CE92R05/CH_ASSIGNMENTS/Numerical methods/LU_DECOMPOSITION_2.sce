clc
clear
function phi=ludcomp(A,r)
//Linear Equation:A*phi=r

[nr_A,nc_A]=size(A)//size of A
[nr_r,nc_r]=size (r)//size of r

if nc_A<>nr_A then
    error('A is not a square matrix')
    abort;
end
if nc_A<>nr_r then
    error('Not compatible matrices')
    abort;
end
n=nc_A
//Decomposition
for k=1:1:n-1
    for i=k+1:1:n
        gam=A(i,k)/A(k,k)
        A(i,k)=gam
        for j=k+1:n
            A(i,j)=A(i,j)-gam*A(k,j)
        end
    end
end

//Forward substitution
psi(1)=r(1)
for i=2:1:n
    sumj=r(i)
    for j=1:i-1
        sumj=sumj-A(i,j)*psi(j)
    end
    psi(i)=sumj
end

//Backward substitution
phi(n)=psi(n)/A(n,n)
for i=n-1:-1:1
    sumj=psi(i)
    for j=i+1:n
        sumj=sumj-A(i,j)*phi(j)
    end
    phi(i)=sumj/A(i,i)
end
endfunction
A=[1 0 0 0 0
 1 2 1 0 0
 0 1 3 -1 0
 0 0 1 2 1
 0 0 0 0 1];
 r=[1
 12
 11
 28
 9];
phi=ludcomp(A,r)
