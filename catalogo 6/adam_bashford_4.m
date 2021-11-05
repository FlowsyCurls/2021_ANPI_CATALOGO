function ejemplo_adam_bashford_4

  clc; clear; close all
  warning ("off")
  a=2; b=4;
  num_pt=11;
  f='1+(x-y)^2';
  y0 = 1;
  
  [xk yk p] = adam_bashford_4(f,a,b,y0,num_pt)
  
  # Prueba de ejemplo (descomentar para comparar)
  #fs = @(x) (x.^2-x-1)./(x-1);
  #xg = a:0.001:b;
  #yg = fs(xg);
  
  #hold on;
  #plot(xg,yg,'g');
  
  
end

function [xk yk p] = adam_bashford_4(f,a,b,y0,num_pt)
  
%  Metodo adam bashford de cuatro pasos
%
%  Esta funcion utiliza el metodo del adam bashford 4 para la solucion de 
%  ecuaciones diferenciales ordinarias.
%
%  Sintaxis: [xk yk p] = adam_bashford_4(f,a,b,num_pt)
%
%  Parametros Iniciales: 
%      f : ecuacion diferencial ordinaria
%      a : inicio del intervalo
%      b : fin del intervalo
%      num_pt: cantidad de puntos
%
%  Parametros Salida: 
%      xk: los valores de xk de los pares ordenados
%      yk: los valores de yk de los pares ordenados
%      p: es el polinomio de forma simbolica obtenido a partir
%         del metodo de Lagrange.
  pkg load symbolic
  f=matlabFunction(sym(f));
  h=(b-a)/(num_pt-1)
  xk=a:h:b
  
  y1_2_3=predictor_corrector(f,a,b,y0,num_pt);
  
  yk=y1_2_3(1:4);
  
  for n=4:num_pt-1  
    
    yk(n+1)=yk(n)+(h/24)*(55*f(xk(n),yk(n)) - 59*f(xk(n-1),yk(n-1)) 
               + 37*f(xk(n-2),yk(n-2)) - 9*f(xk(n-3),yk(n-3)));
    
  end
  
  p = lagrange(xk,yk);
  hold on
  plot(xk,yk,'r')
  title('Metodo de adam bashford 4')
  xlabel('xk')
  ylabel('yk')
  grid on
  
end

function yk = predictor_corrector(f,a,b,y0,num_pt)
  
%  Metodo predictor corrector
%
%  Esta funcion utiliza el metodo del predictor corrector para la solucion de 
%  ecuaciones diferenciales ordinarias.
%
%  Sintaxis: euler(f, a, b, num_pt)
%
%  Parametros Iniciales: 
%      f : ecuacion diferencial ordinaria
%      a : inicio del intervalo
%      b : fin del intervalo
%      num_pt: cantidad de puntos
%
%  Parametros Salida: 
%      xk: los valores de xk de los pares ordenados
%      yk: los valores de yk de los pares ordenados
%      p: es el polinomio de forma simbolica obtenido a partir
%         del metodo de Lagrange.
  pkg load symbolic
  f=matlabFunction(sym(f));
  
  h=(b-a)/(num_pt-1);
  xk=a:h:b;
  yk=[y0];
  zv=[y0];

  for n=1:num_pt-1  
    
    zv(n+1) = yk(n)+h*f(xk(n),yk(n));
    yk(n+1)=yk(n)+(h/2)*f(xk(n),yk(n) +f(xk(n+1),zv(n+1)));
    
  end
  
  
end

function p=lagrange(xv,yv)  
  
%  Metodo de Lagrange
%
%  Esta funcion utiliza el metodo de Lagrange para realizar una 
%  aproximacion del polinomio de interpolacion que pasa por los 
%  puntos (x0, y0),..., (xn, yn) recibidos como parametro.
%
%  Sintaxis: lagrange(xk, yk)
%
%  Parametros Iniciales: 
%      xv : es el vector de preimagenes de cada par ordenado.
%      yv : es el vector de imagenes de cada par ordenado.
%
%  Parametros Salida: 
%      p: es el polinomio de forma simbolica obtenido a partir
%         del metodo de Lagrange.
% 
  
  pkg load symbolic;
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




