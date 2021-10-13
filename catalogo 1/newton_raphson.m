function ejemplo_newton_raphson
  
  clc; clear;
  
  f='exp(x)+x-2';
  x0=5;
  tol=10**-9;
  iterMax=1000;
  
  [xk err]=newton_raphson(f,x0,tol,iterMax)
  
end

function [xk err]=newton_raphson(f,x0,tol,iterMax)
%Esta funcion aproxima la solucion de la ecuacion f(x)=0, 
%utilizando el metodo de Newton Raphson
    %
    %Sintaxis:  [xk err]=newton_raphson(f,x0,tol,iterMax)
    % 
    %Parametros Iniciales: 
    %            f = una  cadena de caracteres (string) que representa a la funcion f
    %            x0 = valor inicial
    %            tol = un numero positivo que representa a la tolerancia para el criterio |f(xk)|<tol
    %            iterMax = cantidad de iteraciones maximas
    %            
    %Parametros de Salida:                           
    %            xk = aproximacion del cero de la funcion f
    %            err =  |f(xk)|
  pkg load symbolic
  f1=matlabFunction(sym(f));
  df = matlabFunction(diff (sym(f)));
  
  er = [];
  err = tol+1;
  
  k=0;
  xk=x0;
  
  while (err>tol) && (k<iterMax)
    
    k=k+1;
    
    %Se verifica si el valor de la derivada es menor que una tolerancia dada, 
    %esto para evitar que se caiga el programar debido a una division entre 0
    if abs(df(xk)) < 10**-15
      disp("Se detuvo el calculo porque se llega a un valor menor a 10**-15 en el denominador")
      break
    endif
    
    xk=xk-f1(xk)/df(xk);
    err=abs(f1(xk));
    er = [er err];
    
  endwhile
  
  
     % Graficar
     plot(0:length(er)-1,er,'g','LineWidth',2)
     set(gca, "fontsize", 20)
     title('Metodo de Newton-Rapshon (Iteraciones vs Error)')
     xlabel('Iteracion (k)')
     ylabel('Error |f(x_k)|')
     grid on
  
end