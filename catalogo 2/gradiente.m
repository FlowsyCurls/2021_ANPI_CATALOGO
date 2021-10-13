function ejemplo_gradiente_conjugado_no_lineal
  clc;clear;
  f="(x-2)^4+(x-2*y)^2";  %funcion
  var="x y"; %variables en la funcion
  x0="0 3"; %valor inicial
  tol=1e-5; %tolerancia
  iterMax=10; %iteraciones maximas
  [xk,error]=gradiente(f,var,x0,tol,iterMax);
  display(["aproximacion = (",num2str(xk(1)),",",num2str(xk(2)) ,")"]);
  error

end

function[xk,error]=gradiente(f,var,x0,tol,iterMax)
%Esta funcion busca resolver numericamente el problema de optimizacion
%de obtener el minimo de una funcion no lineal.
    %
    %Sintaxis:  [xk err]=gradiente(f,var,x0,tol,iterMax)
    % 
    %Parametros Iniciales: 
    %            f = una  cadena de caracteres (string) que representa a la funcion f
    %            var = cadena con las variables de f separadas por un espacio
    %            x0 = valor inicial
    %            tol = un numero positivo que representa a la tolerancia para el criterio |f(xk)|<tol
    %            iterMax = cantidad de iteraciones maximas
    %            
    %Parametros de Salida:                           
    %            xk = aproximacion del cero de la funcion f
    %            error =  |f(xk)|
  pkg load symbolic;
  warning('off', 'all');
  f=sym(f);

  eval(["syms " var ";"]);
  eval(["vars=[" var "];"]);
  eval(["xk=[" x0 "]';"]);
  
  
  grdf=gradient(f,vars); %gradiente de f
  grs=sym(grdf);  %gradiente simbolico
  gk=double(subs(grs,vars,xk)); %gradiente de f (como vector columna)
  dk=-gk; % direccion de descenso 
  
  k=0;
  error=0;
  er = [];
  
  while k<iterMax
    
    delta=rand(1);
    ak=1;
    aux=xk+ak*dk;
    izq=double(subs(f,vars,aux));
    der=delta*ak*gk'*dk;
    
    while izq>der
      
      ak=ak/2;
      aux=xk+ak*dk;
      izq=double(subs(f,vars,aux)-subs(f,vars,xk));
      der=delta*ak*gk'*dk;
      
    end
    
    xk=xk+ak*dk;
    error=norm(double(subs(grs,vars,xk)));
    er = [er error];
    
    if abs(error)<=tol
      
      error=abs(error);
      
      return
    end
    
    gkp1=double(subs(grs,vars,xk));
    beta=norm(gkp1)^2/norm(gk)^2;
    gk=gkp1;
    dk=-gk+beta*dk;
    k=k+1;
    
  end
  
  error=abs(error);
  er = [er error];
  
  plot(er)
  title('Metodo de gradiente conjugado no lineal (Iteraciones vrs Error)')
  xlabel('Iteraciones (k)')
  ylabel('|f(x_k)|')
  grid on
  
  return;
end