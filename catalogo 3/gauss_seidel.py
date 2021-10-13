import matplotlib.pyplot as plt
import  sympy as sp
import numpy as np


def sustitucion_atras(A, b):
    # Esta funcion aplica la sustitucion hacia atras de una funcion Ax=b
    # Parametros de entrada: Matriz A (cuadrada, triangular e invertible).
    #        Vector b de tamano m.
    # Parametros de salida: Vector x de tamano m, que es solucion del sistema Ax=b.
    m = len(b)
    x=np.zeros((1,m))
    for i in range(m-1,-1,-1):
        aux=0
        for j in range(i+1,m):
            aux+=A[i,j]*x[0,j]
        aux=(1/ A[i,i])*(b[i] - aux)
        x[0,i]=aux
    return x[0]


def sustitucion_adelante(A, b):
    # Esta funcion aplica la sustitucion hacia adelante de una funcion Ax=b
    # Parametros de entrada: Matriz A (cuadrada, triangular e invertible).
    #        Vector b de tamano m.
    # Parametros de salida: Vector x de tamano m, que es solucion del sistema Ax=b.
    m=len(b)
    x = np.zeros((1, m))
    for i in range(m):
        aux=0
        for j in range(i):
            aux+=A[i,j]*x[0,j]
        aux= (1/A[i,i])*(b[i] - aux)
        x[0,i]=aux
    return x[0]


def diagonal_dominante(A):
    # Esta funcion verifica si una matriz A es diagonal dominante
    # Parametros de entrada: A=matriz que se verificara
    # Parametros de salida: True si es diagonal dominante y False de lo contrario
    for i in range(0,len(A)):
        d=np.abs(A[i,i])
        f=np.abs(A[i])
        if d <= np.sum(f)-d:
            return False

    return True


def is_squared(A):
    # Esta funcion verifica si una matriz A es cuadrada
    # Parametros de entrada: A=matriz que se verificara
    # Parametros de salida: True si la matriz A es cuadrada y False de lo contrario
    if A.shape[0] == A.shape[1]:
        return True
    print("la matriz A no es cuadrada")
    return False

def gauss_seidel(A,b,x0,tol,iterMax):
    # Esta funcion resuelve un sistema de ecuaciones mediante el metodo iterativo de Guass-Seidel
    # Parametros de entrada: A=matriz de coeficientes, b=vector de terminos independientes, x0=valor inicial,
    #                           tol=tolerancia, iterMax=iteraciones maximas
    # Parametros de salida: (aproximacion, error)

    if not is_squared(A):
        print("La matriz NO es cuadrada")
        return

    if not diagonal_dominante(A):
        print("La matriz NO es diagonalmente dominante")
        return


    L = np.tril(A,-1)
    D = np.diag(np.diag(A))
    U = A-L-D
    LpD = L+D

    y=sustitucion_adelante(LpD,b) #  (L+D)^-1*b    ---->   (L+D)y=b

    xk=x0
    error = []


    for k in range(0,iterMax):

        z=sustitucion_adelante(LpD,U*xk) #  (L+D)^-1*U*xk     ---->  (L+D)z=U*xk
        z=np.transpose(z)

        xk=-z+y
        xk=xk.reshape((-1,1))

        err =np.linalg.norm(A*xk-b)
        error.append(err)

        if err<tol:
            print("El error del metodo es " + str(err))
            print("La cantidad de iteraciones es " + str(k))
            fig, graf = plt.subplots()
            ejex = np.arange(0, k+1, 1)
            graf.plot(ejex, error)
            graf.set_xlabel('Iteraciones ($k$)')
            graf.set_ylabel('Error del sistema')
            graf.set_title('Metodo de Gauss-Seidel (Iteraciones vrs Error)');
            graf.grid(True)
            plt.show()

            return (xk, err)

    return (xk, err)


A=np.matrix(([4,-2,1],
             [2,-17,12],
             [-1,13,17]))

b=np.matrix(([11],
             [70],
             [17]))
iterMax=1000
tol=1e-10
x0=np.ones((len(A),1))
xk=gauss_seidel(A,b,x0,tol,iterMax)
print("El valor aproximado de xk:")
print(xk[0])
