clc
clear
function phi=gausselim(A,r)
//Linear Equation: A*phi=r
[nr_A,nc_A]=size(A)//size of A
[nr_r,nc_r]=size(r)//size of r
//" [nr_A,nc_A]=size(A) " retrieves the number of rows and columns of a matrix A and assigns them to the variables "nr_A" and "nc_A", respectively.

if nc_A<>nr_A then
    error('A is not a square matrix')
    abort;
 //if nc_A is not equal to nr_A , it means that matrix A is not square
end
if nc_A<>nr_r then
    error('Not compatible matrices')
    abort;
end
n=nc_A
//Forward Elimination
for k=1:1:n-1
   for i=k+1:1:n
    gam=A(i,k)/A(k,k)
    for j=k+1:n
        A(i,j)=A(i,j)-gam*A(k,j)
    end
    r(i)=r(i)-gam*r(k)
end
end
//Backward Substitution
phi(n)=r(n)/A(n,n)
for i=n-1:-1:1
    sumj=r(i)
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
 phi=gausselim(A,r)
disp(phi)
