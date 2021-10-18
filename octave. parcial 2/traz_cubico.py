import numpy as np
import sympy as sp


# ----------------------- Auxiliares de Thomas ----------------------- #

def is_cuadrada(A):
    # Funcion que verifica si una matriz de entrada es cuadrada.
    #
    # Sintaxis: is_cuadrada(A)
    #
    # Parametros Iniciales:
    #   A: Matriz de tamano mxn.
    #
    # Parametros de Salida:
    #   True: si la matriz corresponde a una cuadrada.
    #   False: en caso contrario.
    return A.shape[0] == A.shape[1]


def is_tridiagonal(A):
    # Funcion que verifica si una matriz de entrada es tridiagonal.
    #
    # Sintaxis: is_tridiagonal(A)
    #
    # Parametros Iniciales:
    #   A: Matriz de tamano mxn.
    #
    # Parametros de Salida:
    #   True: si la matriz corresponde a una tridiagonal.
    #   False: en caso contrario.

    if not is_cuadrada(A):
        return False

    m = len(A)

    for i in range(m):
        for j in range(m):
            if j > i + 1 or j < i - 1:
                if A[i, j] != 0:
                    return False

    return True


# ----------------------- Auxiliares de Thomas ----------------------- #

def thomas(A, b):
    # Esta funcion encuentra las soluciones de un sistema de ecuaciones
    # utilizando el metodo de Thomas.
    #
    # Sintaxis:  thomas(A,b)
    #
    # Parametros Iniciales:
    #    A: Matriz A de tamano m (cuadrada, triangular superior e invertible).
    #    b: Vector b de tamano m.
    #
    # Parametros Salida:
    #    x : Vector x de tamano m, que es solucion del sistema Ax=b.

    if not is_tridiagonal(A):
        print("Matriz no es tridiagonal.")
        return

    m = len(A)
    b = b.T
    x = np.zeros((m, 1))
    p = np.zeros((m - 1, 1))
    q = np.zeros((m, 1))

    for i in range(m):
        if i == 0:
            p[i] = A[i, i + 1] / A[i, i]
            q[i] = b[i] / A[i, i]
        else:
            aux = A[i, i] - p[i - 1] * A[i, i - 1]
            q[i] = (b[i] - q[i - 1] * A[i, i - 1]) / aux

            if i < (m - 1):
                p[i] = A[i, i + 1] / aux

    n = m - 1
    x[n] = q[n]
    for i in range(n - 1, -1, -1):
        x[i] = q[i] - p[i] * x[i + 1]

    return x


# ----------------------- Auxiliares de traz_cubico ----------------------- #

def tridiagonal(h):
    # Funcion auxiliar de traz_cubico que crea la matriz tridiagonal A utilizada
    # en el metodo de Frontera Natural o Libre.
    #
    # Sintaxis: tridiagonal(h)
    #
    # Parametros Iniciales:
    #    h = vector sujeto con elementos hi = x_(i+1) - x_i
    #
    # Parametros Salida:
    #    A = Matriz tridiagonal

    n = len(h)
    A = np.zeros((n - 1, n - 1))

    for i in range(n - 1):
        A[i, i] = 2 * (h[i] + h[i + 1])
        if 0 < i <= n - 2:
            A[i, i - 1] = h[i]
            A[i - 1, i] = h[i]

    return A


def calc_vect_u(yk, h):
    # Funcion auxiliar de traz_cubico que crea el vector u
    # utilizado en el metodo de Frontera Natural o Libre.
    #
    # Sintaxis:  calc_vect_u(yk, h)
    #
    # Parametros Iniciales:
    #   yk = vector de las imagenes del polinomio
    #   h = vector h calculado en la funcion calc_h
    #
    # Parametros Salida:
    #   u = vector creado a partir de un criterio establecido
    #       por el metodo del trazado cubico de Frontera Libre

    n = len(h)

    # Se suma uno ya que el tamaño deseado toma en cuenta el cero
    u = np.zeros((1, (n - 2) + 1))[0]  # Ej: n=3 -> tamano = n-2 = 1 = [0,1]

    # Se pone (n-2)+1 ya que el for de python va desde 0 hasta n-1
    for i in range(0, (n - 2) + 1):
        aux = 6 * ((yk[i + 2] - yk[i + 1]) / h[i + 1] - (yk[i + 1] - yk[i]) / h[i])
        u[i] = aux

    return u


def calc_h(xk, n):
    # Funcion auxiliar de traz_cubico que crea el vector h
    # utilizado en el metodo de Frontera Natural o Libre.
    #
    # Sintaxis:  calc_h(xk, n)
    #
    # Parametros Iniciales:
    #   xk = vector de las preimagenes de los puntos del polinomio
    #   n = tamano del vector xk
    #
    # Parametros Salida:
    #   h = vector creado a partir de un criterio establecido
    #       por el metodo del trazado cubico de Frontera Libre
    h = np.zeros((1, n - 1))[0]
    for i in range(n - 1):
        h[i] = xk[i + 1] - xk[i]

    return h


# ----------------------- Auxiliares de traz_cubico ----------------------- #

def traz_cubico(xk, yk):
    # Funcion utilizada para construir una funcion a trozos a partir
    # del metodo del Trazador Cubico con los puntos (x0,y0) ,..., (xn,yn)
    #
    # Sintaxis: traz_cubico(xk, yk)
    #
    # Parametros Iniciales:
    #   xk = vector de las preimagenes de los puntos del polinomio
    #   yk = vector de las imagenes del polinomio
    #
    # Parametros Salida:
    #   Sk = vector donde cada entrada es cada polinomio del trazador
    #        cubico.
    x = sp.Symbol('x')

    n = len(xk)

    # P1: Calcular h
    h = calc_h(xk, n)

    # P2: Calcular el vector u
    u = calc_vect_u(yk, h)  # Vector U del metodo

    # P3: Crear la matriz tridiagonal A
    A = tridiagonal(h)

    # P4: Resolver el sistema Am=u
    # m = M1 ... M_(n-1)
    m = thomas(A, u)

    # P5: Construir M a partir de m
    M = m.T[0]
    M2 = np.zeros((1, len(M) + 2))[0]

    for i in range(1, len(M2) - 1):
        M2[i] = M[i - 1]

    M = M2

    # P6 y P7:

    Sk = []

    for i in range(0, n - 1):
        # P6: Calcular constantes
        a = (M[i + 1] - M[i]) / (6 * h[i])
        b = M[i] / 2
        c = (yk[i + 1] - yk[i]) / h[i] - (h[i] * (M[i + 1] + 2 * M[i])) / 6
        d = yk[i]

        # P7: Calcular los trazadores
        S_i = a * (x - xk[i]) ** 3 + b * (x - xk[i]) ** 2 + c * (x - xk[i]) + d
        S_i = sp.expand(S_i)
        Sk.append(S_i)

    return Sk


# Ejemplo presentacion 9 diapositiva 27:
# Puntos:
#       xk      yk
#  k=0 (1,    2.718282)
#  k=1 (1.05, 3.286299)
#  k=2 (1.07, 3.527609)
#  k=3 (1.1,  3.905416)

xk = np.array([1, 1.05, 1.07, 1.1])
yk = np.array([2.718282, 3.286299, 3.527609, 3.905416])
print()
print("Ejemplo para los pares ordenados:")
print("\nxk: " + str(xk))
print("yk: " + str(yk))

Sk = traz_cubico(xk, yk)

print("\nTrazadores:\n")

for i in range(len(Sk)):
    print("S" + str(i) + "(x) = ", end="")
    print(Sk[i])
print("\n\n")