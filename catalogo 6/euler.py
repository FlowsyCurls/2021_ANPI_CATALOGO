import matplotlib.pyplot as plt
import numpy as np
import sympy as sp


def euler(f, a, b, num_pt):
    """
    Metodo de euler

    Esta funcion utiliza el metodo de euler para la solucion de 
    ecuaciones diferenciales ordinarias.

    Sintaxis: euler(f, a, b, num_pt)

    Parametros Iniciales: 
        f : ecuacion diferencial ordinaria
        a : inicio del intervalo
        b : fin del intervalo
        num_pt: cantidad de puntos

    Parametros Salida: 
        p_ordenados: vector de pares ordenados (xk,yk)
        p: es el polinomio de forma simbolica obtenido a partir
           del metodo de Lagrange.
    """

    x = sp.Symbol("x")
    y = sp.Symbol("y")

    f = sp.lambdify([x, y], f)

    h = (b-a)/(num_pt-1)

    xv = np.arange(a, b+h, h)
    yv = [0.5]

    for n in range(num_pt-1):
        yv.append(yv[n]+h*f(xv[n], yv[n]))

    p = lagrange(xv, yv)
    
    print("xk:")
    print(xv)
    print("yk")
    print(yv)

    print("\n➤  El polinomio de interpolacion resultante es:\n")
    sp.pprint(p)

    fig, graf = plt.subplots()
    graf.plot(xv, yv, 'r')
    graf.set_xlabel('xk')
    graf.set_ylabel('yk')
    graf.set_title('Metodo de euler (Grafica de polinomio de interpolacion)')
    graf.grid(True)
    plt.show()

    # Se retorna una lista con los pares ordenados (xk,yk) y el polinomio de interpolacion
    return [xv, yv, p]


def lagrange(xk, yk):
    """
    Metodo de Lagrange

    Esta funcion utiliza el metodo de Lagrange para realizar una 
    aproximacion del polinomio de interpolacion que pasa por los 
    puntos (x0, y0),..., (xn, yn) recibidos como parametro.

    Sintaxis: lagrange(xk, yk)

    Parametros Iniciales: 
        xk : es el vector de preimagenes de cada par ordenado.
        yk : es el vector de imagenes de cada par ordenado.

    Parametros Salida: 
        p: es el polinomio de forma simbolica obtenido a partir
           del metodo de Lagrange.
    """
    x = sp.symbols('x')  # Variable simbolica.

    n = len(xk)-1  # Grado maximo del polinomio de interpolacion.
    p = 0  # Para el calculo del primer polinomio.
    for k in range(n+1):
        # Formula.
        p = p+yk[k]*do_Lk(xk, k)
    p = sp.expand(p)
    return p


def do_Lk(xk, k):
    """
        Funcion que calcula el Lk para un k especifico.
    """
    x = sp.symbols('x')  # Variable simbolica.
    # k= 0,1,2,3,...,n
    n = len(xk)-1  # Grado maximo del polinomio de interpolacion.
    Lk = 1  # Valor de 1 para que NO afecte resultados.
    for j in range(0, n+1):  # Es excluyente.
        if j != k:
            # Formula.
            Lk = Lk*(x-xk[j])/(xk[k]-xk[j])
    Lk = sp.expand(Lk)
    return Lk


f = "y-x^2+1"
a = 0
b = 5
num_pt = 11

euler(f, a, b, num_pt)

