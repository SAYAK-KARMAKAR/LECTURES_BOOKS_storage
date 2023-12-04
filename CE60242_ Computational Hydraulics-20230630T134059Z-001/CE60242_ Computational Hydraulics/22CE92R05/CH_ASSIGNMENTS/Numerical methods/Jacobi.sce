clc
clear
function[count,rmse,phi]= jacobi(A,r,phio,eps_max)
//Linear Equation: A*phi=r

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

while rmse>eps_max
    rmse=0
    for i=1:n
        res(i)=r(i)
        for j=1:n
            res(i)=res(i)-A(i,j)*phio(j)
        end
        phi(i)=phio(i)+res(i)/A(i,i)
        rmse=rmse+(res(i)/A(i,i)).^2
    end
    phio=phi
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
phio=[0
0
0
0
0];
eps_max=1e-6
[count,rmse,phi]=jacobi(A,r,phio,eps_max)
disp(phi)





//Explanation of the loop
//phio = [0; 0; 0]
//res1 = r(1) - A(1,1)*phio(1) - A(1,2)*phio(2) - A(1,3)*phio(3)
//    = 1 - 2*0 - 1*0 - 1*0
//    = 1
//phi(1) = phio(1) + res1/A(1,1)
//       = 0 + 1/2
//       = 0.5
//
//res2 = r(2) - A(2,1)*phio(1) - A(2,2)*phio(2) - A(2,3)*phio(3)
//    = 1 - 1*0 - 2*0 - 1*0
//    = 1
//phi(2) = phio(2) + res2/A(2,2)
//       = 0 + 1/2
//       = 0.5
//
//res3 = r(3) - A(3,1)*phio(1) - A(3,2)*phio(2) - A(3,3)*phio(3)
//    = 1 - 1*0 - 1*0 - 3*0
//    = 1
//phi(3) = phio(3) + res3/A(3,3)
//       = 0 + 1/3
//       = 0.333
//
//RMSE = (res1/A(1,1))^2 + (res2/A(2,2))^2 + (res3/A(3,3))^2
//     = (1/2)^2 + (1/2)^2 + (1/3)^2
//     = 7/36
//
