clear
function dydx = f(x,yn)
    dydx = (0.0008 - (0.09/(15*yn)/(15+2*yn)^4/3 * (15 * yn*yn))) / (1 - (0.181/yn^3))
endfunction
x0 = 0;
y0 = 0.8;
xf = 200;
delx = 10;
x(1)=x0;
y(1)=y0;
n=(xf-x0)/delx;
for i= 1:n
    x(i+1)=x(i)+ delx
     k1 = delx * f(x(i),y(i))
    k2= delx * f(x(i)+ delx, y(i) + k1)
   
    y(i+1)=y(i) + (k1+k2)/2
    
end
x0 = 0;
yc0 = ((20/15)^2 /9.81)^1/3;
xf = 200;
delx = 10;
x(1)=x0;
yc(1)=yc0;
n=(xf-x0)/delx;
for i= 1:n
    x(i+1)=x(i)+ delx
    yc(i+1)=yc(i)

    end
plot2d(x,y,21)
disp(x,y)
plot2d(x,yc,16)
disp(x,yc)
