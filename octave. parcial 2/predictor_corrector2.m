# Ejemplo
function ejemplo
     # Previos
     clc; clear;  close all;  warning ("off");
     pkg load symbolic

     # Parametros
     a=0; b=5;
     N=11;
     y0 = 0.5;     
     
     # Funcion a derivar
     f=@(x,y) y-x.^2+1;
     
     # Aproximacion.
     [pares_ordenados, pol]=predictor_corrector(f, a, b, y0, N)
          
     # 6. Conociendo la solucion analitica se realiza la 
     #  comparacion para ver el error graficamente.
     
     # Solucion analitica (wolfram alpha)
     df=@(x) (x+1).^2-0.5*exp(x);
     %Graficar la solucion
     xs=a:0.0001:b;
     ys=df(xs);
     plot(xs,ys,'g')

endfunction

function [pares_ordenados, pol] = predictor_corrector(f, a, b, y0, N)
     #{
    Esta funcion utiliza el metodo Predictor Corrector para aproximar
    la solucion de ecuaciones diferenciales ordinarias.
    
    Sintaxis:  predictor_corrector(f, a, b, yo, N)
    
    Parametros Iniciales: 
        f : una cadena de caracteres (string) que representa
            a la ecuacion diferencial ordinaria.
        a : limite inferior del intervalo.
        b : limite superior del intervalo.
        y0 : valor inicial de y correspondiente a x0.
        N: cantidad de puntos.
        
        
    Parametros de Salida: 
        pares_ordenados: vector de pares ordenados (xk,yk)
        pol : es el polinomio de forma simbolica obtenido a partir
           del metodo de diferencias divididas de Newton.
    #}
    
     # 1. Calculo de h.
     h = (b-a)/(N-1);

     # 2. Calculo de los xk y yk.
     xk = a:h:b;
     yk = [y0];
     z = [];
     for n=1:N-1
          # 2.1 Calculo del Predictor (euler).
          z(n+1) = yk(n) + h * f(xk(n), yk(n));
          # 2.1 Calculo del Corrector.
          yk(n+1) = yk(n) + (h/2) * (f(xk(n), yk(n))+ f(xk(n+1), z(n+1)));
     end
     
     # 3. Calculo de los pares ordenados.
     n= length(xk);
     pares_ordenados = zeros(n, 2);
     
     for i =1: n
          pares_ordenados(i,1) = [xk(i)];
          pares_ordenados(i,2) = [yk(i)];
     end
     pares_ordenados;
     
     # 4. Calculo del polinomio de interpolacion
     # haciendo uso del metodo de lagrange.
     format short;
     DD = dd_newton(xk, yk);
     pol = calc_polynomial(DD, xk, yk);

     
     # 5. Graficacion de la aproximacion.
     hold on
     plot(xk, yk, 'r')
     title('Metodo Predictor-Corrector (Grafica de polinomio de interpolacion)')
     xlabel('xk')
     ylabel('yk')
     grid on
     
     

endfunction


function DD =dd_newton(x,y)
     n = length(x)-1;   %  Grado maximo del polinomio de interpolacion.
     DD=zeros(n+1);    % Generar matriz de ceros para almacenar diferencias divididas.
     DD(:,1)=y;               % La columna 1 corresponde a los valores de y. Por la regla f(xi) =yi
     % Se inicia la recursion en la columna 2, puesto que  ya fue asignada la 1.
     for j=2:n+1 
          for i=j:n+1
               % Por la regla:
               % f( xi, x[i+1], .... , x[i+j] ) =   (   f( x[i+1], .... , x[i+j] )  -  f( xi, x[i+1], .... , x[i+j-1])   ) / (x[i]  - x[i-j])
               DD(i,j) = [ DD(i,j-1) - DD(i-1,j-1) ] / [ x(i) - x(i-j+1) ];

          end
     end
endfunction

function p = calc_polynomial(DD, X)
     syms x;
     n = length(X)-1;   %  Grado maximo del polinomio de interpolacion.
     p = DD(1,1);
     poly=1;                         % Valor de 1 para que NO afecte resultados.
     for i=1:n
          poly=poly*(x-X(i));
          p = p+poly*DD(i+1,i+1);
     end
     p = expand(p);
endfunction



