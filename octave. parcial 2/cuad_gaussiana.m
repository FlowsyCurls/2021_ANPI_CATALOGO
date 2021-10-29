function ejemplo_cuad_gaussiana
  clc;
  clear;
  format long;
  
  f = 'ln(x)';
  n = 2;
  a = 2;
  b = 5;
  
  [I cota] = cuad_gaussiana(f,n,a,b)
 
endfunction

function [I cota] = cuad_gaussiana(f,n,a,b)
  
  #{
  Esta funcion aproxima el valor para la integral definida de f en el 
  intervalo [a,b] y calcula la cota de cota que se genera al realizar la aproximacion
  utilizando el metodo de la Cuadratura Gaussiana.

  Sintaxis:  [I cota] = cuad_gaussiana(f,n,a,b)

  Parametros Iniciales: 
      f : una cadena de caracteres (string) que representa
          a la funcion f, la cual es continua e integrable en el
          intervalo [a,b].
      n: orden de la cuadratura gaussiana.
      a : limite inferior del intervalo sobre el cual se aplica
          la integral definida.
      b : limite superior del intervalo sobre el cual se aplica
          la integral definida.
          
  Parametros de Salida: 
      aprox :  aproximacion de la funcion en el intervalo 
          indicado.
      cota : cota del cota de la aproximacion.
  #}
  
  pkg load symbolic
  warning('off','all');
  
  syms x;
  
  I = 0;
  cota = 0;

  fs=sym(f);
  y=((b-a)*x+(b+a))/2;
  gs=(b-a)/2*subs(fs,x,y);
  gn=matlabFunction(gs);
  
  [x, w]=ceros_pesos_cuad_gauss(n);
  
  for i=1:n  

    I += w(i)*gn(x(i)); % Cuadratura gaussiana
    
  endfor
  
  
  if n == 2

    % 1. Calculo de la cuarta derivada.
    f2_s = abs(diff(gs, 4));   # d4 simbolica.
    f2_n = matlabFunction(f2_s);   # d2 numerica.

    % 2. Calculo de las funciones auxiliares: 
    % min{ -f } en [a,b] -> max{ f } en [a,b].     
    fs_aux = -1*abs(f2_s);   # -f simbolica.
    fn_aux = matlabFunction(fs_aux);   # -f numerica.

    % 3. Calculo de alpha_max.
    x_max = fminbnd(fn_aux, -1, 1);   # maximo.
    alpha = f2_n(x_max);   # alpha max.

    % 4. Calculo de la cota de cota.
    cota=alpha/135;  # cota.
    
  end  

  
end

function [x, w]=ceros_pesos_cuad_gauss(n)
  if n>10 | n<2
    x=0; w=0;
    display('El valor de n debe ser menor o igual a 10 y mayor o igual a 2')
  else   
    x=[]; w=[];
    if n==2
      x(1)=-0.5773502692; x(2)=0.5773502692;
      w(1)=1; w(2)=1;
    elseif n==3
      x(1)=-0.7745966692; x(2)=0; x(3)=0.7745966692;
      w(1)=0.5555555555; w(2)=0.888888888; w(3)=w(1);
    elseif n==4
      x(1)=-0.86113631159405;
      x(2)=-0.339981043584856;      
      x(3)=0.339981043584856;      
      x(4)=0.86113631159405;
      
      w(1)=0.347854845137454;
      w(2)=0.652145154862546;
      w(3)=w(2);
      w(4)=w(1);
    elseif n==5
      x(1)=-0.9061798459;
      x(2)=-0.5384693101;
      x(3)=0;
      x(4)=0.5384693101;
      x(5)=0.9061798459;
      
      w(1)=0.2369268851;
      w(2)=0.4786286705;
      w(3)=0.5688888889;
      w(4)=0.4786286705;
      w(5)=0.2369268851;
    elseif n==6
      x(1)=-0.9324695142;
      x(2)=-0.6612093865;
      x(3)=-0.2386191861;
      x(4)=0.2386191861;
      x(5)=0.6612093865;
      x(6)=0.9324695142;
      
      w(1)=0.1713244924;
      w(2)=0.3607615730;
      w(3)=0.4679139346;
      w(4)=w(3);
      w(5)=w(2);
      w(6)=w(1);
    elseif n==7
      x(1)=-0.9491079123;
      x(2)=-0.7415311856;
      x(3)=-0.4058451514;
      x(4)=0;
      x(5)=0.4058451514;
      x(6)=0.7415311856;
      x(7)=0.9491079123;
      
      w(1)=0.1294849662;
      w(2)=0.2797053915;
      w(3)=0.3818300505;
      w(4)=0.4179591837;
      w(5)=w(3);
      w(6)=w(2);
      w(7)=w(1);
    elseif n==8
      x(1)=-0.9602898565;
      x(2)=-0.7966664774;
      x(3)=-0.5255324099;
      x(4)=-0.1834346425;
      x(5)=0.1834346425;
      x(6)=0.5255324099;
      x(7)=0.7966664774;
      x(8)=0.9602898565;
      
      w(1)=0.1012285363;
      w(2)=0.2223810345;
      w(3)=0.3137066459;
      w(4)=0.3626837834;
      w(5)=w(4);
      w(6)=w(3);
      w(7)=w(2);
      w(8)=w(1);
    elseif n==9
      x(1)=-0.9681602395;
      x(2)=-0.8360311073;
      x(3)=-0.6133714327;
      x(4)=-0.3242534234;
      x(5)=0;
      x(6)=0.3242534234;
      x(7)=0.6133714327;
      x(8)=0.8360311073;
      x(9)=0.9681602395;
      
      w(1)=0.0812743883;
      w(2)=0.1806481607;
      w(3)=0.2606106964;
      w(4)=0.3123470770;
      w(5)=0.3302393550;
      w(6)=w(4);
      w(7)=w(3);
      w(8)=w(2);
      w(9)=w(1);
    elseif n==10
      x(1)=-0.9739065285;
      x(2)=-0.8650633667;
      x(3)=-0.6794095683;
      x(4)=-0.4333953941;
      x(5)=-0.1488743390;
      x(6)=0.1488743390;;
      x(7)=0.4333953941;
      x(8)=0.6794095683;
      x(9)=0.8650633667;
      x(10)=0.9739065285;
      
      w(1)=0.0666713443;
      w(2)=0.1494513492;
      w(3)=0.2190863625;
      w(4)=0.2692667193;
      w(5)=0.2955242247;
      w(6)=w(5);
      w(7)=w(4);
      w(8)=w(3);
      w(9)=w(2);
      w(10)=w(1);
    end
  end
end