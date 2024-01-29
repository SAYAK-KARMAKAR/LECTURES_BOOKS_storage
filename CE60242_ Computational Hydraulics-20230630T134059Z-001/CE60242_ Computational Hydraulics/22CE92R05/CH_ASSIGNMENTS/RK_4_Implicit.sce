clear
//Runge kutta 4th order implicit
x0=0;y0=0.8
xf=200,h=10
n=(xf-x0)/h
function dydx=f(x,y)
    dydx=(0.0008-0.09*(15+2*y)^(4/3)/(15*y)^(10/3))/(1-611.62/(15*y)^3)
endfunction 
x(1)=x0
y(1)=y0
for i=1:n
    x(i+1)=x(i)+h
    k1=h*f(x(i)+0.211325*h,y(i)+0.25*k1-0.038675k2)
    k2=h*f(x(i)+0.788675*h,y(i)+0.538675*k1+0.25*k2)
    y(i+1)=y(i)+(1/2)*(k1+k2)
end
z=[x,y]
disp(z)
plot(x,y)
