% Ejemplo y previos.
function ejemplo
     clc; clear;
     pkg load symbolic
     
     f = 'sin((pi*x)/2)';
     ptos = [-1 0 1 2];    % Valores de soporte en Intervalo [-1, 2]
     printf("La cota de error del polinomio de interpolacion de esta funcion\
     \nen el intervalo de %d a %d es de:\n\n", ptos(1), ptos(end));
     
     cota = cota_poly_inter(f, ptos)
end

% Cota de Error Polinomio de Interpolacion
function cota = cota_poly_inter(f, ptos)
    #{
    Esta funcion calcula el error que se genera al realizar una aproximacion 
    con un polinomio de interpolacion en el intervalo [a,b] y utilizando un 
    vector soporte que contiene puntos dentro del intervalo mencionado.
    
    Sintaxis:  cota_poly_inter(f,ptos)
    
    Parametros Iniciales: 
        f : una cadena de caracteres (string) que representa
            a la funcion f
        ptos : es un vector que contiene los puntos de soporte
            dentro de un intervalo determinado.

    Parametros Salida: 
        cota :  corresponde a la cota del error
    #}
     n=length(ptos)-1;             # Grado del polinomio.
    
     fs=sym(f);     % Funcion simbolica.
     a = ptos(1); b = ptos(end);   % Extremos del intervalo.
     alpha = alpha(fs, n+1, a, b);  % Alpha
  
     x_max = polynomial_max(ptos, a, b);
     suma = 1;
     for k=ptos
          % k+1: indices en Octave empiezan en 1.
          suma *= abs(x_max - k);
     end
     cota = alpha / factorial(n+1) * suma;   % Cota
endfunction

function alpha = alpha(fs, m, a, b)
     fn = matlabFunction(fs);   % Funcion numerica.
     dfm_s = diff(fs, m);     % Derivada n+1 de f simbolica.
     dfm_n = matlabFunction(dfm_s);     % Derivada n+1 de f numerica.
     
     % Calculo del maximo utilizando funcion del minimo sobre -f(x).
     % min{ -f } en [a,b] -> max{ f } en [a,b].
     fs_aux = -1*abs(dfm_s);
     fn_aux = matlabFunction(fs_aux); 
     x_max = fminbnd(fn_aux, a, b);
     alpha = abs(dfm_n(x_max));
endfunction

function x_max = polynomial_max(ptos, a, b)
     syms x;
     
     ## Creacion  polinomio simbolico con los puntos.
     f = 1;
     for i = ptos
          f *= (x - i);
     end
     
     f = expand(f);
     fs=sym(f);     % Funcion simbolica.
     fn = matlabFunction(fs);
     f1 = diff(f, x); % Primera derivada del polinomio.
     crit_pts = flip(solve(f1)'); %Puntos criticos.
     X = [ fn(a), fn(b), crit_pts ]; % extremos relativos.
     
     Y = [];
     for k=1 : length(X)
          punto = double(X(k));
          if (a < punto) & (punto < b)
               % k+1: indices en Octave empiezan en 1.
               Y = [ Y, abs(fn(punto)) ];
          else
               X(k) = []; # Se elimina.
          endif
     end
     [y_max index] = max(Y);
     x_max = double(X(index));

endfunction
