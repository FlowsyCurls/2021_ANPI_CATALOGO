# Ejemplo
function ejemplo
     # Previos
     clc; clear;
     pkg load symbolic
     warning('off','all');
     
     # Parametros
     f='log(x)';
     a=2;
     b=5;
     N = 500;
     [aprox, error]=trapecio_compuesto(f, N, a, b);
     
     printf(
     "El valor aproximado de la funcion \'f(x) = %s\' utilizando el metodo del Trapecio\
     \nCompuesto en el intervalo de %d a %d con un total de %d puntos es de: \n%d\
     \n\nLa cota de error correspondiente a esta aproximacion es: \n%d\n\n",
     f, a, b, N, aprox, error
     
     );
endfunction

function [I, cota] = trapecio_compuesto(f, N, a, b)
#{
    Esta funcion aproxima el valor para la integral definida de f en el 
    intervalo [a,b] y calcula el error que se genera al realizar la aproximacion
    utilizando la Regla del Trapecio Compuesto con un total de N puntos.
    
    Sintaxis:  trapecio_compuesto(f, N, a, b)
    
    Parametros Iniciales: 
        f : una cadena de caracteres (string) que representa
            a la funcion f, la cual es continua e integrable en el
            intervalo [a,b].
        N: numero total de puntos en los que se debe 
            dividir el intervalo.
        a : limite inferior del intervalo sobre el cual se aplica
            la integral definida.
        b : limite superior del intervalo sobre el cual se aplica
            la integral definida.
            
    Parametros de Salida: 
        aprox :  aproximacion de la funcion en el intervalo 
            indicado.
        cota : cota del error de la aproximacion.
    #}
    
     x=sym('x');  # variable simbolica.
     fs=sym(f);  # funcion simbolica.
     fn = matlabFunction(fs);  # funcion numerica.
     
     # Aproximacion utilizando el metodo del Trapecio Compuesto.
     
     % 1. Calculo de h.
     h=(b-a)/(N-1);
     
     % 2. Calcular el vector de puntos. 
     xv=linspace(a, b, N);
     
     % 3. Realizar la aproximacion de la integral.
     I = 0;
     for i=1:N-1
          ai = xv(i);
          bi = xv(i+1);
          fai = fn(ai);
          fbi = fn(bi);
          I+=(bi-ai)*(fai+fbi)/2;
     end
     
     # Cota del error.
     
     % 1. Calculo de la segunda derivada.
     f2_s = abs(diff(fs, 2));   # d2 simbolica.
     f2_n = matlabFunction(f2_s);   # d2 numerica.
     
     % 2. Calculo de las funciones auxiliares: 
     % min{ -f } en [a,b] -> max{ f } en [a,b].     
     fs_aux = -1*abs(f2_s);   # -f simbolica.
     fn_aux = matlabFunction(fs_aux);   # -f numerica.
     
     % 3. Calculo de alpha_max.
     x_max = fminbnd(fn_aux, a, b);   # maximo.
     alpha = f2_n(x_max)   # alpha max.
     
     % 4. Calculo de la cota de error.
     cota=(((b-a)*h^2)/12)*alpha;  # cota.
     
endfunction