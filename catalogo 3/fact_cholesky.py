import numpy as np
from math import *

# Input: Matriz A (cuadrada, triangular e invertible).
#        Vector b de tamano m.
# Output: Vector x de tamano m, que es solucion del sistema Ax=b.
def sustitucion_atras(A, b):
    m = len(b)
    x=np.zeros(m)   #Crea un vector de ceros de tamano n
    for i in range(m-1,-1,-1):
        aux=0
        for j in range(i+1,m):
            aux+=A[i,j]*x[j]
        aux=(1/ A[i,i])*(b[i] - aux)
        x[i]=aux
    return x

# Input: Matriz A (cuadrada, triangular e invertible).
#        Vector b de tamano m.
# Output: Vector x de tamano m, que es solucion del sistema Ax=b.
def sustitucion_adelante(A, b):
    m=len(b)
    x = np.zeros(m)
    for i in range(m):
        aux=0
        for j in range(i):
            aux+=A[i,j]*x[j]
        aux= (1/A[i,i])*(b[i] - aux)
        x[i]=aux
    return x

# Input: Matriz A (cuadrada)
# Output: True si la matriz A es simetrica y False de lo contrario
def is_symetric(A):

    a_t=A.transpose()
    c= A == a_t # Compara las matrices (a y a_t)mxm y devuelve una matriz mxm 
                # con todos los valores en True si son iguales

    if c.__contains__(False):
        print("la matriz no es simetrica")
        return False
    return True

# Input: Matriz A
# Output: True si la matriz A es cuadrada y False de lo contrario
def is_squared(A):
    if A.shape[0] == A.shape[1]:
        return True
    print("la matriz A no es cuadrada")
    return False

# Input: Matriz A
# Output: True si la matriz es positiva definida y False de lo contrario
def is_positive(A):

    m = A.shape[0]

    for i in range(m):
        aux = A[0:i+1, 0:i+1]
        if np.linalg.det(A)<=0: 
            print("la matriz no es positiva definida")
            return False
    return True

# Input: Matriz A (cuadrada, simetrica y positiva)
# Output: Matriz L
def l_cholesky(A):
    m = len(A)
    L = np.zeros(A.shape)

    # Los valores de la columna 0 de L
    j = 0
    L[j, j] = np.sqrt(A[j, j])
    for i in range(1, m):
        L[i, j] = A[i, j] / L[j, j]
    # Los valores de la columna 1,...,m-1 de L
    for j in range(1, m):
        aux1 = 0
        for k in range(0, j):
            aux1 = aux1 + np.power((L[j, k]),2)
        L[j, j] = sqrt(A[j, j] - aux1)

        for i in range(j + 1, m):
            aux2 = 0
            for k in range(0, j):
                aux2 = aux2 + L[i, k] * L[j, k]
            L[i, j] = (A[i, j] - aux2) / (L[j, j])

    return L

# Input: Matriz A , vector b.
# Output: Vector x que es solucion del sistema Ax=b.
def fact_cholesky(A, b):

    if is_squared(A) and is_symetric(A) and is_positive(A):

        L=l_cholesky(A)
        print("L:")
        print(L)
        L_t=L.transpose()
        y=sustitucion_adelante(L, b)
        x=sustitucion_atras(L_t,y)

        print("La solucion del sistema es:")
        print(x)
        return x

    else:
        return None

a=np.matrix(([25,15,-5,-10],
             [15,10,1,-7],
             [-5,1,21,4],
             [-10,-7,4,18]))
b=np.matrix(([-25],
             [-1],
             [-21],
             [-5]))

fact_cholesky(a,b)
