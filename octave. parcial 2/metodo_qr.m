# Ejemplo.
function ejemplo
     # Previos.
     clc; clear;  close all;  warning ("off");
     pkg load symbolic

     # Parametros.
     A=[0 11 -5;
            -2 17 -7;
            -4 26 -10];

     iterMax=100;
     tol=10**-10;
     
     # Aproximacion.
     [valores, vectores, k, error]=metodo_qr(A, iterMax, tol);
     
     # Print.
     n=  length(valores);
     display("METODO QR\n")
     for j=1: n
          fprintf('£f%d : %d  y  ', j, valores(j));
          fprintf('£h%d : [', j);
          fprintf('%g, ', vectores(j,1:end-1));
          fprintf("%g]\'\n\n", vectores(j,end));
    end
    fprintf('Iteraciones : %d\n', k);
    fprintf('Error : %d\n\n', error);

endfunction

function [valores, vectores, n, err] = metodo_qr(A, iterMax, tol)
     #{
    Esta funcion utiliza el metodo de Potencia para calcular el
    valor propio de mayor magnitud una matrix de entrada 
    simetrica definida positiva asi como el vector propio 
    respectivo a este.
    
    Sintaxis:   metodo_qr(A, iterMax, tol)
    
    Parametros Iniciales: 
        A :  matriz simetrica a la cual se le deben 
        calcular los valores y vectores propios.
        tol = un numero positivo que representa a la tolerancia 
        para el criterio |xkn - xk|<tol
        iterMax = cantidad de iteraciones maximas
        
    Parametros de Salida: 
        valores : valores propios de la matriz.
        vectores : matriz de vectores propio respectivo a cada valor propio calculado.
        error : el error obtenido bajo el criterio de |Ukn - Uxk|.
    #}
    n = size(A)(1);
    Ak = A;
    U = eye(n);
    error = [];

    for n=1 : iterMax
         
         [Q, R] = qr(A);
         U = U*Q;
         A = R*Q;
         
         Vp = diag(A);
         Vpk = diag(Ak);
         Ak = A;
         
         err = norm(Vpk-Vp)/norm(Vpk);  % Error.
         error = [ error, err ]; 

         % Error.
         if err < tol
              break
          endif
     endfor

    % Salida.
    valores  = diag(A)';
    vectores = U';
    err = error(end);

    % Graficar
    plot(0:length(error)-1,error,'g','LineWidth',1)
    set(gca, "fontsize", 17)
    title('Metodo de Potencia (Iteraciones vs Error)')
    xlabel('Iteracion (k)')
    ylabel('Error |(xkn - xk)|')
    grid on

endfunction




