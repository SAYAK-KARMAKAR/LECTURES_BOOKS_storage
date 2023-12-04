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
    k2= delx * f(x(i)+ 0.5* delx, y(i) +0.5 *k1)
   
    y(i+1)=y(i) + k2
    
end

plot2d(x,y,5)
disp(x,y)
for i= 1:n
    x(i+1)=x(i)+ delx;
    k1= delx * f(x(i),y(i))
    k2= delx * (f(x(i)+ 2/3 *delx ,y(i)+ 2/3 *k1))
    y(i+1)=y(i) + (k1+3*k2)/4
end
plot2d(x,y,6)
disp(x,y)
for i= 1:n
    x(i+1)=x(i)+ delx;
    k1= delx * f(x(i),y(i))
    k2= delx * f((x(i)+ delx/2) ,(y(i)+ k1*0.5))
     k3= delx * f((x(i)+ delx/2) ,(y(i)+ k2*0.5))
     k4= delx * f(x(i)+ delx, y(i) +k3)
    y(i+1)=y(i) + (k1+ 2*k2 + 2*k3+ k4)/6 
end
plot2d(x,y,11)
disp(x,y)
for i= 1:n
    x(i+1)=x(i)+ delx
    k1 =delx * f(x(i+1),y(i))
    y(i+1)= y(i) + k1
    
    end

plot2d(x,y,14)
disp(x,y)
for i= 1:n
    x(i+1)=x(i)+ delx;
    k1= delx * f(x(i),y(i))
    k2= delx * (f(x(i)+ 0.5 *delx ,y(i)+ 0.5 *k1))
    k3= delx * f(x(i)+ delx , y(i) - k1+ 2*k2)
    y(i+1)=y(i) + (k1+ 4*k2 + k3)/6
end
plot2d(x,y,3)
disp(x,y)
for i= 1:n
    x(i+1)=x(i)+ delx
     k1 = delx * f(x(i),y(i))
    k2= delx * f(x(i)+ delx, y(i) + k1)
   
    y(i+1)=y(i) + (k1+k2)/2
    
end

plot2d(x,y,21)
disp(x,y)
for i= 1:n
    x(i+1)=x(i)+ delx
     k1 = delx * f(x(i),y(i))
    k2= delx * f(x(i)+ delx, y(i) + k1)
   
    y(i+1)=y(i) + (k1+k2)/2
    
end

plot2d(x,y,21)
disp(x,y)
for i= 1:n
    x(i+1)=x(i)+ delx;
    k1= delx * f(x(i),y(i))
    k2 = delx * f(x(i) + 0.5*delx, y(i) + 0.5*( delx * f(x(i),y(i))))
    k3 = delx * f(x(i) +0.5* delx, y(i) + 0.5* (delx * f(x(i) + 0.5*delx, y(i) +0.5)))
    k4 = delx * f(x(i) + delx, y(i) + (delx * f(x(i) +0.5* delx, y(i) + 0.5* (delx * f(x(i) + 0.5*delx, y(i) +0.5)))))
     
    y(i+1)=y(i) + k1/6 + k2/3 +k3/3 +k4/6
end
plot2d(x,y,11)
disp(x,y)
for i= 1:n
    x(i+1)=x(i)+ delx;
    k1= delx * f(x(i)+ 0.5 * delx  ,y(i)+ 0.5*(f(x(i),y(i))))

    y(i+1)=y(i) + k1
end
plot2d(x,y,20)
disp(x,y)
