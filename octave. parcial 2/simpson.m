# Ejemplo
function ejemplo
     # Previos
     clc; clear;
     pkg load symbolic
     
     # Parametros
     f='ln(x)';
     a=2;
     b=5;
     [aprox, error] =  simpson(f,a,b);
     
     display("\n* METODO DE SIMPSON\n");
     printf(
     "El valor aproximado de \'f(x) = %s\' utilizando el metodo de Simpson\
     en el intervalo de %d a %d es de: \n%d\n\
     \nLa cota de error correspondiente a esta aproximacion es: \n%d\n\n",
     f, a, b, aprox, error
     );
endfunction

function [aprox, cota] = simpson(f, a, b)
#{
    Esta funcion aproxima el valor para la integral definida de f en el 
    intervalo [a,b] y calcula el error que se genera al realizar la aproximacion
    utilizando la Regla de Simpson (a traves de un polinomio de grado 2).
    
    Sintaxis:  simpson(f, a, b)
    
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
     c = b-a;
     aprox = (c/6) * ( fn(a) + 4*fn((a+b)/2) + fn(b) );
     
     # Cota del error de la aproximacion.
     d = 4; # indice de 4ta derivada.
     f4_s = diff(fs, d);   # 4ta derivada simbolica.
     f4_n = matlabFunction(f4_s);   # 4ta derivada numerica.
     # min{ -f } en [a,b] -> max{ f } en [a,b].     
     fs_aux = -1*abs(f4_s);   # -f simbolica.
     fn_aux = matlabFunction(fs_aux);   # -f numerica.
     x_max = fminbnd(fn_aux, a, b);   # maximo.
     alpha = abs(f4_n(x_max));   # alpha.
     cota = ( ( c**5 )/2880 )*alpha; # cota.
endfunction