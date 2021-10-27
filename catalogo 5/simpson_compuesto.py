import numpy as np
import sympy as sp
import scipy as sc


def simpson_compuesto(f, a, b, N):
    x = sp.Symbol('x')
    f = sp.sympify(f)
    h = (b - a) / (N - 1)
    suma_par = 0
    suma_impar = 0
    x0 = a
    i = 1
    while i < N - 1:
        xi = x0 + i * h
        if i % 2 == 0:  # si es par
            suma_par += sp.N(f.subs(x, xi))
        else:  # si es impar
            suma_impar += sp.N(f.subs(x, xi))
        i += 1

    I = (h / 3) * (sp.N(f.subs(x, x0)) + 2 * suma_par + 4 * suma_impar + sp.N(f.subs(x, b)))

    #error = np.

    return I


# Ejemplo
f = 'ln(x)'
a = 2
b = 5
N = 7

I = simpson_compuesto(f, a, b, N)
print("Aproximacion:")
print(I)
print("Cota de error:")
