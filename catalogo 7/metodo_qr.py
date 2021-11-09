from sympy import *
from scipy import linalg
import numpy as np
import matplotlib.pyplot as plt

# Ejemplo
error = []
def ejemplo():
    # Parametros.
    A = np.matrix(([0, 11, -5], [-2, 17, -7], [-4, 26, -10]))
    iterMax=25
    tol=10**-10

    # Aproximacion.
    valores, vectores, k, err = metodo_qr(A, iterMax, tol)
    print("\n\n ➤   METODO QR\n")

    for i in range(len(vectores)):
        string = "λ{0} : {1}   y   ν{0} : [ {2}".format(i+1, valores[i], vectores[i][0])
        for j in range(1, len(vectores)):
            string += ", {0}".format(vectores[i][j])
        string += " ]'\n"
        print(string)
    print("Iteraciones : ", k)
    print("\nError : ", err)
    print("\n")

    fig, graf = plt.subplots()
    ejex = np.arange(1, k+1)
    graf.plot(ejex, error)
    graf.set_xlabel('Iteraciones ($k$)')
    graf.set_ylabel('Error')
    graf.set_title('Metodo QR (Iteraciones vs Error)');
    graf.grid(True)
    plt.show()

def metodo_qr(A, iterMax, tol):    
    """
    Esta funcion utiliza el metodo QR para calcular los
    valores propios de una matrix de entrada simetrica
    asi como los vectores propios respectivos a estos.
    
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


    i = iterMax
    for i in range(iterMax):
        [Q, R] = np.linalg.qr(A)
        Ak = R*Q
        U = U*Q

        Vpk = np.diag(Ak)
        Vp = np.diag(A)
        A = Ak;
        
        err = np.linalg.norm(Vpk-Vp)/np.linalg.norm(Vpk)  # Error.
        error.append(err) 

        # Error.
        if err < tol:
            break

    # Salida.
    valores  = np.diag(A)
    vectores = np.transpose(U)
    err = error[-1]
    
    return valores, vectores, i+1, err


ejemplo()