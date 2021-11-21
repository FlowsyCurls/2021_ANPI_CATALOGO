import numpy as np
import sympy as sp
import scipy as sc
from scipy import optimize

# Ejemplo
def ejemplo():
    # Parametros
    f = '2^x'
    a = 0
    b = 1
    N = 9
    aprox, error =  simpson_compuesto(f, a, b, N)

    print("\n\n ➤   METODO DE SIMPSON COMPUESTO\n")
    string = ("El valor aproximado de \'f(x) = {0}\' utilizando el metodo de Simpson Compuesto\
        \nen el intervalo de {1} a {2} con un total de {3} puntos es de:\n ●  {4}\n\
        \nLa cota de error correspondiente a esta aproximacion es:\n ●  {5}\
        \n").format(f, a, b, N, aprox, error)
    print(string)
    
def simpson_compuesto(f, a, b, N):
    """
    Esta funcion aproxima el valor para la integral definida de f en el 
    intervalo [a,b] y calcula el error que se genera al realizar la aproximacion
    utilizando el metodo de Simpson Compuesto con un total de N puntos.
    
    Sintaxis:  simpson_compuesto(f, a, b, N)
    
    Parametros Iniciales: 
        f : una cadena de caracteres (string) que representa
            a la funcion f, la cual es continua e integrable en el
            intervalo [a,b].
        a : limite inferior del intervalo sobre el cual se aplica
            la integral definida.
        b : limite superior del intervalo sobre el cual se aplica
            la integral definida.
        N: numero total de puntos en los que se debe dividir el 
            intervalo.
            
    Parametros Salida: 
        I :  aproximacion de la funcion en el intervalo 
            indicado.
        cota : cota del error de la aproximacion.
    """
    x = sp.Symbol('x')
    f = sp.sympify(f)
    h = (b - a) / (N - 1)
    print(h)
    print()
    suma_par = 0
    suma_impar = 0
    x0 = a
    i = 1
    while i < N - 1:
        xi = x0 + i * h
        print(xi)
        if i % 2 == 0:  # si es par
            suma_par += sp.N(f.subs(x, xi))
        else:  # si es impar
            suma_impar += sp.N(f.subs(x, xi))
        i += 1

    print(suma_par)
    print(suma_impar)



    I = (h / 3) * (sp.N(f.subs(x, x0)) + 2 * suma_par + 4 * suma_impar + sp.N(f.subs(x, b)))


    # Cota del error.
    
    # 1. Calculo de la cuarta derivada.
    f4_s = sp.Abs(f.diff(x, 4))   ## d4 simbolica.
    print(f4_s)
    f4_n = sp.lambdify(x, f4_s)   ## d4 numerica.
    print(f4_n)

    # 2. Calculo de las funciones auxiliares:
    # min{ -f } en [a,b] -> max{ f } en [a,b].
    fs_aux = -1*sp.Abs(f4_s)   ## -f simbolica.
    fn_aux = sp.lambdify(x, fs_aux)   ## -f numerica.
    
    # 3. Calculo de alpha_max.
    x_max = optimize.fminbound(fn_aux, a, b)   ## maximo.
    alpha = f4_n(x_max)   ## alpha.

    # 4. Calculo de la cota de error.
    numerador = (b-a)*h**4
    cota = (numerador/180)*alpha   ## cota.
    
    return I, cota

ejemplo()