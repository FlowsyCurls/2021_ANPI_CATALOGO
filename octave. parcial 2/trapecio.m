# Ejemplo
function ejemplo
     # Previos
     clc; clear;
     pkg load symbolic
     
     # Parametros
     f='ln(x)';
     a=2;
     b=5;
     [aprox, error]=trapecio(f,a,b);
     
     printf(
     "El valor aproximado de la funcion \'f(x) = %s\' utilizando el metodo del\
     \nTrapecio en el intervalo de %d a %d es de: %d\
     \n\nLa cota de error correspondiente a esta aproximacion es: %d\n\n",
     f, a, b, aprox, error
     );
endfunction

function [aprox, cota] = trapecio(f, a, b)
#{
    Esta funcion aproxima el valor para la integral definida de f en el 
    intervalo [a,b] y calcula el error que se genera al realizar la aproximacion
    utilizando la Regla del Trapecio (a traves de un polinomio de grado 1).
    
    Sintaxis:  trapecio(f, a, b)
    
    Parametros Iniciales: 
        f : una cadena de caracteres (string) que representa
            a la funcion f, la cual es continua e integrable en el
            intervalo [a,b].
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
     
     # Aproximacion.
     aprox=((b-a)/2)*(fn(a)+fn(b));  # calculo de la integral aproximada.
     
     # Cota del error de la aproximacion.
     d = 2; # Segunda derivada.
     f2_s = diff(fs, d);   # segunda derivada simbolica.
     f2_n = matlabFunction(f2_s);   # segunda derivada numerica.
     # min{ -f } en [a,b] -> max{ f } en [a,b].     
     fs_aux = -1*abs(f2_s);   # -f simbolica.
     fn_aux = matlabFunction(fs_aux);   # -f numerica.
     x_max = fminbnd(fn_aux, a, b);   # maximo.
     alpha = abs(f2_n(x_max));   # alpha.
     cota=(((b-a)^3)/12)*alpha;  # cota.
endfunction