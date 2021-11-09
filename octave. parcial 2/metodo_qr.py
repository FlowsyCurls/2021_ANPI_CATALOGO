from sympy import *
from scipy import optimize
import numpy as np
from math import *

# Ejemplo
def ejemplo():
    # Parametros.
    A = np.matrix(([0, 11, -5], [-2, 17, -7], [-4, 26, -10]))
    iterMax=25
    tol=10**-10

    # Aproximacion.
    valores, vectores = metodo_qr(A, iterMax, tol)
    print("\n\n ➤   METODO QR\n")

    for i in range(len(vectores)):
        string = "λ{0} : {1}   y   ν{0} : [ {2}".format(i, valores[i], vectores[i][0])
        for j in range(1, len(vectores)):
            string += ", {0}".format(vectores[i][j])
        string += " ]'\n"
        print(string)
    print()

def metodo_qr(A, iterMax, tol):    
    """
    Esta funcion utiliza el metodo de Potencia para calcular el
    valor propio de mayor magnitud una matrix de entrada 
    simetrica definida positiva asi como el vector propio 
    respectivo a este.
    
    Sintaxis:   metodo_qr(A, iterMax, tol)
    
    Parametros Iniciales:
        A :  matriz simetrica a la cual se le deben 
        calcular los valores y vectores propios.
        tol = un numero positivo que representa a la tolerancia 
        para el criterio |xkn - xk|<tol
        iterMax = cantidad de iteraciones maximas
        
    Parametros de Salida:
        valores : valores propios de la matriz.
        vectores : matriz de vectores propio respectivo a cada valor 
            propio calculado.
        error : el error obtenido bajo el criterio de |Ukn - Uxk|.
    """
    n = A.shape[0]
    U = eye(n)
    error = []

    for i in range(iterMax):
        [Q, R] = np.linalg.qr(A)
        A = R*Q
        U = U*Q

##         error = [ error, norm(Uk_n-Uk) ];  % Error.
##         Uk = Uk_n;
##         % Error.
##         if error < tol
##              break

    # Salida.
    valores  = np.diag(A)
    vectores = np.transpose(U)
    err = 0
    
    return valores, vectores


ejemplo()