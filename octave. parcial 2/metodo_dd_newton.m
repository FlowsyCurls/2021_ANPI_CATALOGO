function ejemplo()
     clc;
     pkg load symbolic;
     x=[-2 0 1];
     y=[0 1 -1];

     disp('La matriz de diferencias divididas de Newton es:');
     DD = dd_newton(x,y)
     disp('El polinomio de Newton es');
     p = calc_polynomial(DD, x,y);
     pretty(p)
end

function DD =dd_newton(x,y)
     n = length(x)-1;   %  Grado maximo del polinomio de interpolacion.
     DD=zeros(n+1);    % Generar matriz de ceros para almacenar diferencias divididas.
     DD(:,1)=y;               % La columna 1 corresponde a los valores de y. Por la regla f(xi) =yi
     % Se inicia la recursion en la columna 2, puesto que  ya fue asignada la 1.
     for j=2:n+1 
          display("Columna");   j-1
          for i=j:n+1
              % i-1
              % j-1
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



