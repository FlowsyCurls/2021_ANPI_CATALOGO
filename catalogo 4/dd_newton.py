from sympy import *

# Ejemplo
def ejemplo():
    xk=[-2, 0, 1] ## Vector de valores X de los pares ordenados.
    yk=[0, 1, -1] ## Vector de valores Y de los pares ordenados.
    Dd = dd_newton(xk,yk)
    p = calc_polynomial(Dd, xk)

    print("\n\n➤  La matriz de diferencias divididas de Newton es:\n")
    pprint(Dd)
    
    print("\n➤  El polinomio de interpolacion resultante es:\n")
    pprint(p)
    print("\n")

# Metodo de Lagrange
def dd_newton(xk, yk):
    n = len(xk)-1          ## Grado maximo del polinomio de interpolacion.
    m = n+1
    Dd = zeros(m, m)  ## Generar matriz de ceros para almacenar diferencias divididas.
    Dd[:,0] = yk           ## La columna 0 corresponde a los valores de y. Por la regla f(xi) =yi
    ## Se inicia la recursion en la columna 2, puesto que  ya fue asignada la 1.
    for j in range(1, m):
        for i in range(j, m):
            ## Por la regla:
            ## f( xi, x[i+1], ..., x[i+j] ) =   
            ## (  f( x[i+1], ..., x[i+j] ) - f( xi, x[i+1], ..., x[i+j-1])  ) / (x[i] - x[i-j])
            #Viendo la matriz es mas sencillo.
            fx1 = Dd[i,j-1]
            fx0 = Dd[i-1,j-1]
            x1 = xk[i]
            x0 = xk[i-j]
            Dd[i,j] = ( fx1 - fx0 ) / ( x1 - x0 )
    return Dd

def calc_polynomial(Dd, X):
    x = symbols('x')   ## Variable simbolica.
    n = len(X)-1       ## Grado maximo del polinomio de interpolacion.

    p = Dd[0,0];
    poly = 1;          ## Valor de 1 para que NO afecte resultados.
    for i in range(n):
        poly *= (x-X[i])
        p += poly * Dd[i+1,i+1]
    p = expand(p)
    return p

ejemplo()