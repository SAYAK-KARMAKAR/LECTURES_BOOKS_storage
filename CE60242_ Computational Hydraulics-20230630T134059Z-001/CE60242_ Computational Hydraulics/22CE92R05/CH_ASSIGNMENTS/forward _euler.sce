clear
//forward Eulers method
function dydx=f(x,y)
    dydx=(0.0008-0.09*(15+2*y)^(4/3)/(15*y)^(10/3))/(1-611.62/(15*y)^3)
endfunction
x0=0;y0=0.8
xf=200,h=10
n=(xf-x0)/h
x=zeros(n,1)
y=zeros(n,1)
x(1)=x0
y(1)=y0
for i=1:n
    x(i+1)=x(i)+h
    y(i+1)=y(i)+h*f(x(i),y(i))
end
z=[x,y]
disp(z)
plot(x,y)
