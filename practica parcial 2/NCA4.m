% Pregunta 2.
function ejemplo
     % Previos
     clc; clear;
     pkg load symbolic
     
     % Parametros
     f='x^2*exp(x^2)';
     a=0;
     b=2;
     n= 50;
     
     display("\nAPROXIMACIONES :");
     
     display("\na) NCA4");
     I=NCA4(f,a,b)
   
     display("\nb) NCA4_compuesto");
     I=NCA4_compuesto(a,b, n)
     display("");
     
endfunction

% Pregunta 2.a
function I = NCA4(f, a, b)
     syms x;
     f = matlabFunction(sym(f));
          
     % 1. Calculo de h.
     h = (b-a)/5;
     
     % 2. Calculo de xj.
     xj = a+(1:4)*h;
     
     % 3. Calculo de la aproximacion de la integral definida.
     I = (5*h)/24 * [ 11*f(xj(1)) + f(xj(2)) + f(xj(3))  + 11*f(xj(4)) ];    
endfunction

% Pregunta 2.b
function I = NCA4_compuesto(a, b, n)
     f='x^2*exp(x^2)';
  
     % 1. Calculo de los puntos.
     h=(b-a)/(n-1);
     ptos = a:h:b;
     
     % 2. Calculo de aproximacion.
     I = 0;
     for i=2:n
          I = I + NCA4(f, ptos(i-1), ptos(i));
     end
endfunction