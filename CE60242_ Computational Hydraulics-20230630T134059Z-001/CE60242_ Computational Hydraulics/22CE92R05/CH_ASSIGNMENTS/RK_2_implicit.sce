clear
//Runge kutta 2nd order implicit
x0=0;y0=0.8
xf=200,h=10
n=(xf-x0)/h
function dydx=f(x,y)
    dydx=(0.0008-0.09*(15+2*y)^(4/3)/(15*y)^(10/3))/(1-611.62/(15*y)^3)
endfunction
function coff1=k1(x,y,k1)
    coff1=
endfunction
x(1)=x0
y(1)=y0
for i=1:n
    x(i+1)=x(i)+h
    k1=h*f(x(i)+0.5*h,y(i)+0.5*k1)
    y(i+1)=y(i)+k1
end
z=[x,y]
disp(z)
plot(x,y)
