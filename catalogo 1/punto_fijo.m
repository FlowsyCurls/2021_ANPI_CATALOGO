function funcion
     clc; clear;
     f ='log(2*x+1)';
     x0=7;
     tol=0.0000001;
     iterMax=100;
     [xk error] = punto_fijo(f, x0, tol, iterMax)
end

function [xk error] = punto_fijo(f, x0, tol, iterMax)
% Metodo numerico de punto fijo iterativo para encontrar al menos un cero  de una funcion f(x)=0
% que cumple previamente con la existencia y unicidad. 
     % Sintaxis:   [xk error] = punto_fijo(f, x0, tol, iterMax)
     % Parametros Iniciales: 
     %         f : funcion a evaluar con el metodo.
     %         x0 : punto inicial que cumple los requisitos de un punto fijo.
     %         tol : valor de la tolerancia a la equivocacion minima aceptable.
     %         iterMax : cantidad maxima de iteraciones a realizar.
    % Parametros de Salida:                           
    %          xk : aproximacion del cero de la funcion f.
    %          err :  |f(xk)|
    
     pkg load symbolic
     g=matlabFunction(sym(f));
     it=[];    %Definir vector de iteraciones. 
     er=[];    %Definir vector de error.
     for k=0:iterMax;
         xk = g(x0);     % Se evalua la funcion con el valor actual.
         error = abs(xk-x0);  % Se calcula el error.
         x0=xk;     % Actualizar el valor.
         it = [it k];   % Introducir resultado de la iteracion en vector respectivo.
         er = [er error];     % Introducir resultado del error en vector respectivo (iteracion vieja).  
         if error<tol    % Evaluar la aceptabilidad del error.
              break
          end
     end

     % Graficar
     plot(0:length(er)-1,er,'g','LineWidth',2)
     set(gca, "fontsize", 20)
     title('Metodo de Punto Fijo (Iteraciones vs Error)')
     xlabel('Iteracion (k)')
     ylabel('Error |f(x_k+1) - f(x_k)|')
     grid on
     
endfunction
