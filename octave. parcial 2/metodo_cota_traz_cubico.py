from sympy import *
import numpy as np
from scipy import optimize

# Ejemplo
def ejemplo():

    f='exp(x/2)';
    ptos=[1, 1.5, 1.75, 2.15, 2.4, 3];  ## Valores de soporte en Intervalo [1, 3]
    cota = cota_traz_cubico(f, ptos)
    ## Notacion cientifica.
    cota = np.format_float_scientific(cota)

    string = ("\n ➤  La cota de error del trazador cubito de esta\
        \n    funcion en el intervalo de {0} a {1} es de: \n\n\t {2} \n\n\
        ").format(ptos[0], ptos[-1], cota);
    print(string)


def cota_traz_cubico(f, ptos):
    """
        Sea f una funcion real, definida, continua y con cuarta derivada continua
        sobre un intervalo cerradora de R.
        Esta funcion calcula el error que se genera al realizar una aproximacion 
        con un trazador cubico de interpolacion en el intervalo [a,b] y
        utilizando un vector soporte que contiene puntos dentro del intervalo
        mencionado.

    Sintaxis:  cota_traz_cubico(f,ptos)
    
    Parametros Iniciales: 
        f : una cadena de caracteres (string) que representa
            a la funcion f
        ptos : es un vector que contiene los puntos de soporte
            dentro de un intervalo determinado.

    Parametros Salida: 
        cota :  corresponde a la cota del error
    """
    # Calculo de h
    n=len(ptos);                    ## Grado maximo del polinomio.
    dist=zeros(n,1);
    for i in range(len(ptos)-1):
        dist[i] = ptos[i+1] - ptos[i]
    h = max(dist)

    x = Symbol('x')                 ## Variable simbolica.
    fs = sympify(f)                 ## Funcion simbolica.
    df4_s =  fs.diff(x, 4)          ## Derivada cuarta de f simbolica.
    df4_n = lambdify(x, df4_s)      ## Derivada cuarta de f numerica.

    ## Extremos del intervalo.
    a = ptos[0]
    b = ptos[-1]

    ## Calculo del maximo utilizando funcion del minimo sobre -f(x).
    ## min{ -f } en [a,b] -> max{ f } en [a,b].
    fs_aux = -1*abs(df4_s)
    fn_aux = lambdify(x, fs_aux)
    x_max = optimize.fminbound(fn_aux, a, b)    ## X max.
    alpha = df4_n(x_max)                        ## Alpha.
    cota = (5*h**4/384)*alpha                   ## Cota.
    return cota
    
ejemplo()