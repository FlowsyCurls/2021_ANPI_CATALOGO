% Ejemplo y previos.
function ejemplo
     clc; clear; 
     pkg load symbolic
     xv=[-2 0 1];
     yv=[0 1 -1];
     p=lagrange(xv,yv)
end

% Metodo de Lagrange
function p=lagrange(xv,yv)
     syms x  % Variable simbolica.
     n=length(xv)-1; %  Grado maximo del polinomio de interpolacion.
     p=0;
     for k=0:n % Indice empieza en 0 para simular parecido a la formula vista.
          p=p+yv(k+1)*do_Lk(xv,k); % k+1 : indices en Octave empiezan en 1.
     end
     p=expand(p);
end

# Funcion auxiliar que calcula el Lk para un k especifico.
function Lk=do_Lk(xv,k)
     syms x
     %k=0,1,....,n
     n=length(xv)-1;
     Lk=1; % Valor de 1 para que NO afecte resultados.
     for j=0:n
          if j~=k
               % k+1 y j+1 : indices en Octave empiezan en 1.
               Lk=Lk*(x-xv(j+1))/(xv(k+1)-xv(j+1));
          end
     end
     Lk=expand(Lk);
end

