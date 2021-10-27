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
    
    # 1. Calculo de la cuarta derivada.
    f4_s = Abs(fs.diff(x, 4))   ## d4 simbolica.
    f4_n = lambdify(x, f4_s)   ## d4 numerica.
    
    # 2. Calculo de las funciones auxiliares: 
    # min{ -f } en [a,b] -> max{ f } en [a,b].
    fs_aux = -1*Abs(f4_s)   ## -f simbolica.
    fn_aux = lambdify(x, fs_aux)   ## -f numerica.
    
    # 3. Calculo de alpha_max.
    x_max = optimize.fminbound(fn_aux, a, b)   ## maximo.
    alpha = f4_n(x_max)   ## alpha.

    # 4. Calculo de la cota de error.
    cota = ( ( c**5 )/2880 )*alpha           ## cota.
    return aprox, cota

ejemplo()