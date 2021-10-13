from sympy import *
init_printing(use_latex=true)


# Ejemplo
def ejemplo():
    xk=[-2, 0, 1] ## Vector de valores X de los pares ordenados.
    yk=[0, 1, -1] ## Vector de valores Y de los pares ordenados.
    p = lagrange(xk,yk)
    print("\n➤  El polinomio de interpolacion resultante es:\n")
    pprint(p)

# Metodo de Lagrange
def lagrange(xk, yk):
    x = symbols('x')   ## Variable simbolica.
    
    n = len(xk)-1      ## Grado maximo del polinomio de interpolacion.
    p = 0              ## Para el calculo del primer polinomio.
    for k in range(n+1):
        # Formula. 
        p=p+yk[k]*do_Lk(xk,k);
    p = expand(p)
    return p

# Funcion que calcula el Lk para un k especifico.
def do_Lk(xk,k):
    x = symbols('x')        ## Variable simbolica.
    # k= 0,1,2,3,...,n
    n = len(xk)-1           ## Grado maximo del polinomio de interpolacion.
    Lk=1                    ## Valor de 1 para que NO afecte resultados.
    for j in range(0, n+1): ## Es excluyente.
        if j!=k:
            # Formula. 
            Lk=Lk*(x-xk[j])/(xk[k]-xk[j])
    Lk = expand(Lk)
    return Lk

# CALL
ejemplo()
print()
