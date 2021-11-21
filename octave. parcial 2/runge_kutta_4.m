function ejemplo
  
  clc;clear;warning off;
  a=2;
  b=10;
  y0=4;
  m=50;
  f='(x+y)/x';
  [x_v,y_v]=runge_kutta_4(f,a,b,y0,m)
  
  # Prueba de ejemplo (descomentar para comparar)
  #fs = "x*(log(x/2)+2)";
  #fs = matlabFunction(sym(fs));
  #xg = a:0.001:b;
  #yg = fs(xg);
  
  #hold on;
  #plot(xg,yg,'g');
  
end

function [x_v,y_v]=runge_kutta_4(f,a,b,y0,m)
  
  pkg load symbolic;
  
  syms x;
  syms y;
  
  fs=matlabFunction(sym(f));
  
  h=(b-a)/(m-1);
  x_v=a:h:b
  y_v=[y0];
  
  for n=1:m-1
    k1 = fs(x_v(n),y_v(n));
    k2 = fs(x_v(n)+h/2,y_v(n)+h*(k1/2));
    k3 = fs(x_v(n)+h/2,y_v(n)+h*(k2/2));
    k4 = fs(x_v(n)+h,y_v(n)+h*k3);
    y_v_1 = y_v(n)+(h/6)*(k1+2*k2+2*k3+k4);
    y_v = [y_v y_v_1];
  endfor
  
  hold on
  plot(x_v,y_v,'r')
  title('Metodo de runge kutta 4')
  xlabel('x_v')
  ylabel('y_v')
  grid on
  
  
  
end