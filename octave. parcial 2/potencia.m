# Ejemplo
function ejemplo
     # Previos
     clc; clear;  close all;  warning ("off");
     pkg load symbolic

     # Parametros
     x0 = [1; 1; 1];
     A=[-3 1 0;
            1 -2 1;
            0 1 -3];

     iterMax=100;
     tol=10**-2;
     
     # Aproximacion.
     [valor, vector, k, error]=potencia(A, x0, iterMax, tol);
     
     # Print.
     fprintf('valor propio\t:  %d\n',valor);
     fprintf('vector propio\t:  [');
     fprintf('%g, ', vector(1:end-1));
     fprintf('%g]\n', vector(end));
     fprintf('iteraciones\t:  %d\n', k);
     fprintf('error aprox.\t:  %d\n', error(end));

endfunction

function [valor, vector, k, err] = potencia(A, x0, iterMax, tol)
     #{
    Esta funcion utiliza el metodo de Potencia para calcular el
    valor propio de mayor magnitud una matrix de entrada 
    simetrica definida positiva asi como el vector propio 
    respectivo a este. 
    
    Sintaxis:   potencia(A, x0, iterMax, tol)
    
    Parametros Iniciales: 
        A :  matriz simetrica definida positiva a la cual se le deben 
        calcular el vector propio.
        x0 : vector inicial aleatorio no nulo.
        tol = un numero positivo que representa a la tolerancia 
        para el criterio |xkn - xk|<tol
        iterMax = cantidad de iteraciones maximas
        
    Parametros de Salida: 
        valor : valor propio de mayor magnitud de la matriz.
        vector : vector propio respectivo al valor propio calculado.
        error : el error obtenido bajo el criterio de |xkn - xk|.
    #}
    
    xk = x0;
    k = iterMax;
    error = [];

    for n=1 : iterMax
         yk = A*xk;
         ck = norm(yk, Inf);  % Norma infito.
         xk_n = (1/ck)*yk;
         err = norm(xk_n-xk);
         error = [ error,  err ];  % Error.
         xk = xk_n;
         

         % Error.
         if err < tol
              k=n;
              break
         end
    end
    
    % Salida.
    valor  = ck;
    vector = xk';
    err = error(end);
    

    % Graficar
    plot(0:length(error)-1,error,'g','LineWidth',1)
    set(gca, "fontsize", 17)
    title('Metodo de Potencia (Iteraciones vs Error)')
    xlabel('Iteracion (k)')
    ylabel('Error |(xkn - xk)|')
    grid on

endfunction




