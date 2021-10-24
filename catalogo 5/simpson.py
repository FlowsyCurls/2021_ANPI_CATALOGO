from sympy import *
from scipy import optimize
import numpy as np

# Ejemplo
def ejemplo():
    # Parametros
    f='ln(x)';
    a=2;
    b=5;
    aprox, error = simpson(f,a,b)


    print("\n\n ➤   METODO DE SIMPSON\n")
    string = ("El valor aproximado de \'f(x) = {0}\' utilizando el metodo\
        \nde Simpson en el intervalo de {1} a {2} es de:\n ●  {3}\n\
        \nLa cota de error correspondiente a esta aproximacion es:\n ●  {4}\
        \n").format(f, a, b, aprox, error)
    print(string)

def simpson(f,a,b):    
    """
    Esta funcion aproxima el valor para la integral definida de f en el 
    intervalo [a,b] y calcula el error que se genera al realizar la aproximacion
    utilizando la Regla de Simpson (a traves de un polinomio de grado 2).
    
    Sintaxis:  trapecio(f, a, b)
    
    Parametros Iniciales: 
        f : una cadena de caracteres (string) que representa
            a la funcion f, la cual es continua e integrable en el
            intervalo [a,b].
        a : limite inferior del intervalo sobre el cual se aplica
            la integral definida.
        b : limite superior del intervalo sobre el cual se aplica
            la integral definida.
            
    Parametros Salida: 
        aprox :  aproximacion de la funcion en el intervalo 
            indicado.
        cota : cota del error de la aproximacion.
    """

    x = Symbol('x')         ## variable simbolica.
    fs = sympify(f)         ## funcion simbolica.
    fn = lambdify(x, fs)    ## funcion numerica.

    ## Aproximacion
    c = b-a
    aprox = (c/6) * ( fn(a) + 4*fn((a+b)/2) + fn(b) )
    
    ## Cota del error de la aproximacion.
    d = 4 # indice cuatro por 4ta derivada.
    f4_s =  fs.diff(x, d)                      ## 4ta derivada simbolica.
    f4_n = lambdify(x, f4_s)                   ## 4ta derivada numerica.

    ## min{ -f } en [a,b] -> max{ f } en [a,b].
    fs_aux = -1*abs(f4_s)                      ## -f simbolica.
    fn_aux = lambdify(x, fs_aux)               ## -f numerica.
    x_max = optimize.fminbound(fn_aux, a, b)   ## maximo.
    alpha = abs(f4_n(x_max))                        ## alpha.
    cota = ( ( c**5 )/2880 )*alpha           ## cota.
    return aprox, cota

ejemplo()