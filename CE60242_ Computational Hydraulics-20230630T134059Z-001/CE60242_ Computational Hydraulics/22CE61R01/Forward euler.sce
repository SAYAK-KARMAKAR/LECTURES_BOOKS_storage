clear
global yn0
yn0=0.8
function dydx = f(x,yn)
    dydx = (0.0008 - (0.09/(15*yn)/(15+2*yn)^4/3 * (15 * yn*yn))) / (1 - (0.181/yn^3))
endfunction
/*function Gyn = f(yn0)
    Gyn = ((0.0008^1/2 * 15^ 5/3 /0.015)* (yn0/15+2*yn0)^2/3 * yn0) - 20
endfunction
function Gpyn= f(yn0)
    Gpyn= (0.0008^1/2 * 15^5/3 * yn0 ^2/3 *(5*15 + 6*yn0))/(3*0.015*(15+2*yn0)^5/3)
endfunction
y=yn0;
    err=1;
    count=0;
    while abs(err)>1
        count=count+1;
        Gyn = f(yn0);
        Gpyn= f(yn0);
        disp([Gyn Gpyn])
        err=-Gyn/Gpyn;
        y=y+err;
        disp([count y])    
    end
disp(y)*/
x0 = 0;
y0 = 0.8;
xf = 200;
delx = 10;
x(1)=x0;
y(1)=y0;
n=(xf-x0)/delx;
for i= 1:n
    x(i+1)=x(i)+ delx
    y(i+1)=y(i) + delx * f(x(i),y(i))

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
plot2d(x,y,12)
disp(x,y)
plot2d(x,yc,16)
disp(x,yc)
plot2d(x,yn0)
