function ejemplo
  clc;clear;warning off;
  f = 'ln(x)';
  a = 2;
  b = 5;
  N = 7;
  [aprox error] = simpson_compuesto(f, a, b, N);
  
end

function [aprox error] = simpson_compuesto(f,a,b,N)
  
  pkg load symbolic;
  syms x;
  
  fs=matlabFunction(sym(f));
  h=(b-a)/(N-1)
  suma_par=0;
  suma_impar=0;
  x0=a;
  i=1;
  while i<N-1
      xi=x0+i*h
      if mod(i,2) == 0  # si es par
          suma_par += fs(xi);
      else  # si es impar
          suma_impar += fs(xi);
      end
      i+=1;
  end
  
  suma_par
  suma_impar
  
  aprox = (h/3)*(fs(x0)+2*suma_par+4*suma_impar+fs(b))
  
  # Cota del error.
  
  # 1. Calculo de la cuarta derivada.g
  
  #f4_s = abs(diff(sym(f), 4))   ## d4 simbolica. 
  
  ## malo de aqui en adelante xd
##  f4_n = lambdify(x, f4_s)   ## d4 numerica.
##
##  # 2. Calculo de las funciones auxiliares:
##  # min{ -f } en [a,b] -> max{ f } en [a,b].
##  fs_aux = -1*abs(f4_s)   ## -f simbolica.
##  fn_aux = lambdify(x, fs_aux)   ## -f numerica.
##
##  # 3. Calculo de alpha_max.`
##  x_max = optimize.fminbound(fn_aux, a, b)   ## maximo.
##  alpha = f4_n(x_max)   ## alpha.
##
##  # 4. Calculo de la cota de error.
##  numerador = (b-a)*h^4
##  error = (numerador/180)*alpha   ## cota.
  error = 0
  
  
end