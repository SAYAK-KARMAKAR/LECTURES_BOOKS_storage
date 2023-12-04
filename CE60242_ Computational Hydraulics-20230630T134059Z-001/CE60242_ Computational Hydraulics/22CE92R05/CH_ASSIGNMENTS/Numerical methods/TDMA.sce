clc
clear
function phi=tdmasolv(b,d,a,r)
    
//n:number of rows
n=length(d);

a(1)=a(1)/d(1);
r(1)=r(1)/d(1);
//forward elimination
for i=2:n-1
    fact=d(i)-b(i)*a(i-1);
    a(i)=a(i)/fact;
    r(i)=(r(i)-b(i)*r(i-1))/fact;
end

r(n)=(r(n)-b(n)*r(n-1))/(d(n)-b(n)*a(n-1));

//Backward substitution
phi(n)=r(n);
for i=n-1:-1:1
    phi(i)=r(i)-a(i)*phi(i+1);
end
endfunction
d=[1 2 3 2 1];
b=[0 1 1 1 0];
a=[0 1 -1 1 0];
r=[1 12 11 28 9];
phi=tdmasolv(b,d,a,r)
disp(phi)