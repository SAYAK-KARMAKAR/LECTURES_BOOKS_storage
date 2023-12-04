//Groundwater 1d using Gauss elimination method
clc
clear
function phi=gausselim(A,r)
//Linear Equation: A*phi=r
[nr_A,nc_A]=size(A)//size of A
[nr_r,nc_r]=size(r)//size of r

if nc_A<>nr_A then
    error('A is not a square matrix')
    abort;
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

//....Problem dependent parameters
nnode=21; //number of nodes
xl=0;//coordinate of left end boundary xL
xr=1000;
cconf=1e-11;
T=2e-5;
hs=90;
c0=90;
c1=0.06;
c2=-0.00003;
//....
x=linspace(xl,xr,nnode);//In between left and right boundary, nnode numbers of x coordinates has been generated.
delx=x(2)-x(1);//gird size
//.....Initialization of matrices...
h=zeros(nnode,1);
A=zeros(nnode,nnode);
r=zeros(nnode,1);
//Left boundary
A(1,1)=1.0;
r(1)=hs;
//Interior nodes
for i=2:nnode-1
    A(i,i-1)=1.0/(delx^2);
    A(i,i)=-((cconf/T)+2.0/(delx^2));
    A(i,i+1)=1.0/(delx^2);
    r(i)=-(cconf/T)*(c0+c1*x(i)+c2*x(i)^2);
end
//Right boundary
//Two point
A(nnode,nnode)=1/delx;
A(nnode,nnode-1)=-1/delx;
r(nnode)=0;

disp(A)
disp(r)
h=gausselim(A,r)
plot(x,h','-r')
disp(max(h))







