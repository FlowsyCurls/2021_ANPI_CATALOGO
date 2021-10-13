import sympy as sp
import numpy as np
import matplotlib.pyplot as plt


def muller(f, valores_iniciales, tol, iterMax):
    """
    Esta funcion aproxima la solucion de la ecuacion f(x)=0,
    utilizando el metodo de Muller.

    Sintaxis:  muller(f,valores_iniciales,tol,iterMax)

    Parametros Iniciales:
        f = una  cadena de caracteres (string) que representa a la funcion f
        valores_iniciales = valores iniciales (x0,x1,x2)
        tol = un numero positivo que representa a la tolerancia
                para el criterio |f(xk)|<tol
        iterMax = cantidad de iteraciones maximas

    Parametros de Salida:
        xk = aproximacion del cero de la funcion f
        error =  |f(xk)|
    """

    x = sp.Symbol('x')
    f1 = sp.sympify(f)

    x0 = valores_iniciales[0]
    x1 = valores_iniciales[1]
    x2 = valores_iniciales[2]

    fx0 = sp.N(f1.subs(x, x0))
    fx1 = sp.N(f1.subs(x, x1))
    fx2 = sp.N(f1.subs(x, x2))

    er = []
    err = tol + 1

    k = 0

    while err > tol and k < iterMax:

        k = k + 1

        a = (fx0 - fx2) / ((x0 - x2) * (x0 - x1)) \
            - ((fx1 - fx2) / ((x1 - x2) * (x0 - x1)))
        b = ((x0 - x2) * (fx1 - fx2)) / ((x1 - x2) * (x0 - x1)) - \
            ((x1 - x2) * (fx0 - fx2) / ((x0 - x2) * (x0 - x1)))
        d = (b + sp.sqrt((b ** 2) - 4 * a * fx2))
        e = (b - sp.sqrt((b ** 2) - 4 * a * fx2))

        if abs(d) < abs(e):
            x3 = x2 - ((2 * fx2) / e)
        else:
            x3 = x2 - ((2 * fx2) / d)

        err = abs(x3 - x2) / abs(x3)
        er.append(err)

        if abs(x3 - x2) < tol:
            fig, graf = plt.subplots()
            ejex = np.arange(1, k + 1, 1)
            graf.plot(ejex, er)
            graf.set_xlabel('Iteraciones ($k$)')
            graf.set_ylabel('$|f(x_k)|$')
            graf.set_title('Metodo de Muller (Iteraciones vrs Error)');
            graf.grid(True)
            plt.show()
            return x3, err

        x0 = x1
        x1 = x2
        x2 = x3
        fx0 = sp.N(f1.subs(x, x0))
        fx1 = sp.N(f1.subs(x, x1))
        fx2 = sp.N(f1.subs(x, x2))


# Ejemplo numerico
f = 'sin(x)-x/2'
valores_iniciales = [2, 2.2, 1.8]
tol = 10 ** -9
iterMax = 100

y = muller(f, valores_iniciales, tol, iterMax)
print("xk: ", y[0])
print("error: ", y[1])
