function ejemplo
  
  clc;clear;
  a=2;
  b=10;
  y0=4;
  m=50;
  [x_v,y_v]=runge_kutta_3(a,b,y0,m)
  
  # Prueba de ejemplo (descomentar para comparar)
  #fn = "x*(log(x/2)+2)";
  #fn = matlabFunction(sym(fn));
  #xg = a:0.001:b;
  #yg = fn(xg);
  
  #hold on;
  #plot(xg,yg,'g');
  
end

function [x_v,y_v]=runge_kutta_3(a,b,y0,m)
  
  pkg load symbolic;
  
  syms x;
  syms y;
  f='(x+y)/x';
  fn=matlabFunction(sym(f));
  
  h=(b-a)/(m-1);
  x_v=a:h:b
  y_v=[y0];
  
  for n=1:m-1
    k1 = fn(x_v(n),y_v(n));
    k2 = fn(x_v(n)+h/2,y_v(n)+h*(k1/2));
    k3 = fn(x_v(n)+h,y_v(n)+h*(2*k2-k1));
    y_v_1 = y_v(n)+(h/6)*(k1+4*k2+k3);
    y_v = [y_v y_v_1];
  endfor
  
  hold on
  plot(x_v,y_v,'r')
  title('Metodo de runge kutta 3')
  xlabel('x_v')
  ylabel('y_v')
  grid on
  stem(x_v,y_v,'LineStyle','None', 'r')  # Para ver cada punto.

  
end