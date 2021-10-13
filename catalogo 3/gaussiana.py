import numpy as np
from math import *


# Input: Matriz A (cuadrada e invertible).
#        Vector b de tamano m.
# Output: Matriz At (cuadrada y triangular superior),
#         Vector bt de tamano m.
def triangular_sup(A,b):
    A_aumt = np.concatenate((A,b), axis=1)    # Matriz aumentada.
    m = len(b)                      # Tamano del vector b.
    for k in range(0,m-1):          # Iteracion hasta la penultima col.
        for i in range(k+1, m):     # Iteracion en las filas.
            m_ik = A_aumt[i,k]/A_aumt[k,k]
            for j in range(k, m+1): # Iteracion en valores de las filas.
                A_aumt[i,j] = A_aumt[i,j] - m_ik*A_aumt[k,j];
    At = A_aumt[:,0:m]
    bt = A_aumt[:,m]
    return (At, bt)


# Input: Matriz A (cuadrada, triangular superior e invertible).
#        Vector b de tamano m.
# Output: Vector x de tamano n, que es solucion del sistema Ax=b.
def sustitucion_atras(A, b):
    m = len(b)                      # Tamano del vector b.
    x=np.zeros((1,m))             # Vector solucion relleno con 0s.
    for i in reversed(range(m)):    # Iteracion de adelante hacia atras.
        aux=0                       # Variable para la sumatoria.
        for j in range(i+1, m):
            aux += A[i,j]*x[0,j]
        x[0,i] = (1/A[i,i])*(b[i] - aux)
    return x

# Input: Matriz A (cuadrada, triangular superior e invertible).
#        Vector b de tamano m.
# Output: Vector x de tamano n, que es solucion del sistema Ax=b.
def gaussiana(A,b):
    if(np.linalg.det(A)!=0):
        (At, bt) = triangular_sup(A,b)      # Obtener matriz triangular.
        return sustitucion_atras(At, bt)    # Realizar sustitucion
    return ([None])




# Ejemplo:
A = np.matrix(([2,-6,12,16],
               [1,-2,6,6],
               [-1,3,-3,-7],
               [0,4,3,-6]))
b = np.matrix('70;26;-30;-26')

x=gaussiana(A,b)
if(x is not None):
    print("\n▶ A :\n", A)
    print("\n▶ b :\n", b)
    print("\n▶ Solucion del sistema \n • x: \n", np.transpose(x))