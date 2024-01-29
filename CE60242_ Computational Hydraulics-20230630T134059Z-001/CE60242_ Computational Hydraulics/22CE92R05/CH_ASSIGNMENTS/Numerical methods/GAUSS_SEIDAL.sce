clc
clear
function[count,rmse,phi]=gseidal(A,r,phio,eps_max,omega)
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

count=0
rmse=1
phi=phio
while rmse>eps_max
    rmse=0
    for i=1:1:n
        resi=r(i)
        for j=1:n
            resi=resi-A(i,j)*phi(j)
        end
        phi(i)=phi(i)+omega*resi/A(i,i)
        rmse=rmse+(omega*resi/A(i,i)).^2
    end
    rmse=sqrt(rmse/n)
    count=count+1;
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
// A2=[1 2 -3 4 5
// 0 3 -5 -7 9
// 5 -4 3 -2 1
// 1 4 -7 -10 13
// -15 13 11 -9 2];
// r2=[37
// 8
// 3
// 13
// 18];
 phio=[0
0
0
0
0];
eps_max=1e-6;
omega=1
[count,rmse,phi]=gseidal(A,r,phio,eps_max,omega)
disp(phi)
disp(rmse)
disp(count)
