% Ejemplo y previos.
function ejemplo
     clc; clear; 
     pkg load symbolic
     
     f='exp(x/2)';
     ptos=[1 1.5 1.75 2.15 2.4 3];  % Valores de soporte en Intervalo [1, 3]
     printf("La cota de error del trazador cubito de esta funcion \nen el intervalo de %d a %d es de:\n\n", ptos(1), ptos(end));
     cota = cota_traz_cubico(f, ptos)
end

% Cota de Error Trazador Cubico Natural
function cota = cota_traz_cubico(f, ptos)
     
     %Calculo de h
     n=length(ptos);
     dist=zeros(n-1,1);
     for i=1:n-1
       dist(i)=ptos(i+1)-ptos(i);
     end
     h=max(dist);
     
     fs=sym(f);                         % Funcion simbolica.
     df4_s = diff(fs, 4);               % Derivada cuarta de f simbolica.
     df4_n = matlabFunction(df4_s);     % Derivada cuarta de f numerica.
     a = ptos(1); b = ptos(end);        % Extremos del intervalo.
     
     % Calculo del maximo utilizando funcion del minimo sobre -f(x).
     % min{ -f } en [a,b] -> max{ f } en [a,b].
     fs_aux = -1*abs(df4_s);
     fn_aux = matlabFunction(fs_aux); 
     x_max = fminbnd(fn_aux, a, b);     % Alpha.
     alpha = df4_n(x_max);
     cota = (5*h^4/384)*alpha;          % Cota.
end