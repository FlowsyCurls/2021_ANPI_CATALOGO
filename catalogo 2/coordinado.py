import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from math import sqrt

# Input: variables = texto con las variables que deben ser simbolicas.
#    NOTA: Las variables deben estar en una posicion impar. Ej: (x|y|z)
# Output: vector con las variables simbolicas.
def variables_simbolicas(variables):
    vrbs=[]
    for i in range(0, len(variables), 2):
        vrbs.append(sp.Symbol(variables[i]))
    return vrbs

# Input: f = funcion, vrbs = vector de n variables, 
# Output: gradiente de la funcion ingresada respecto a las variables
def gradiente(f,vrbs):
    g=[]
    for var in vrbs:
        g.append(sp.diff(f,var))
    return g

# Input: gs = funcion gradiente simbolica, vrbs = vector de n variables.
# Output: cada funcion del gradiente evaluada con los valores de constante recibidos.
def eval_gradiente(gs, vrbs, vals):
    gn = []
    for i in range(len(vrbs)):
        gn.append(gs[i].subs(zip(vrbs,vals)))
    return gn





# Input: gn = gradiente numerico.
# Output: norma del gradiente.
def norm(gn):
    suma=0
    for comp in gn:
        suma+=comp**2
    return sqrt(suma)

# Funcion que encuentra la preimagen del minimo de una funcion.
# Input: fs - funcion simbolica, pos_var: posicion de la variables
# vrbs -  variables simbolicas, vals - valores iniciales.
def fmin_preimg(fs, pos_var, vrbs, vals):
    ## Sustituir TODA variable que sea la analizada por su valor constante.
    for i in range(0, len(vrbs)):
        if i!=pos_var:
            fs = fs.subs(vrbs[i], vals[i])
    fn=sp.lambdify(vrbs[pos_var], fs)  # Evaluar la variable de forma numerica.
    minimum = optimize.fmin(fn, 0, disp=False)   # Encontrar el minimo.
    return minimum[0]

# Grafica el error.
def graph(k, err):
    fig, graf = plt.subplots()
    graf.plot(k, err)
    graf.set_xlabel('Iteraciones ($k$)')
    graf.set_ylabel('Error $|f(x_(k+1))|$')
    graf.set_title('Metodo de Descenso Coordinado');
    graf.grid(True)
    plt.show()

def coordinado(f, variables, iniciales, tol, iterMax):
    """
    Esta funcion de optimizacion utiliza un algoritmo que minimiza 
    sucesivamente a lo largo de las direcciones de coordenadas para 
    encontrar el minimo de una funcion de entrada.
        En cada iteracion, el algoritmo determina una coordenada a 
    traves de una regla de seleccion (en este caso se utiliza la 
    regla de Gauss-Seidel) de coordenadas y luego minimiza de 
    manera exacta mientras fija todas las demas coordenadas.
    
    Sintaxis:  coordinado(f,variables,iniciales,tol,iterMax)
    
    Parametros Iniciales: 
        f : una cadena de caracteres (string) que representa
            a la funcion f
        variables : es un string con las varibles simbolicas
            a tomar en cuenta (solo se toman en cuenta las
            posiciones impares).
        iniciales : es un string con los valores iniciales 
            de las variables.
        tol : un numero que representa la tolerancia para 
            el criterio |gradiente(f(x_k))|<tol
        iterMax :  numero maximo de iteraciones permitidas.

                
    Parametros Salida: 
        [xk, error, k] donde:                
        xk : vector de coordenadas aproximadas.
        k : numero de iteraciones realizadas
        error :  |gradiente(f(x_k))|
    """

    # Funcion numerica.
    fs=sp.sympify(f)

    # Variables
    vrbs = variables_simbolicas(variables)

    # Aproximaciones y error
    gs = gradiente(f, vrbs)     # Gradiente simbolico.
    error = 0
    xk, it, err = [], [], []
    tmp = iniciales.split(" ")
    for i in range(len(tmp)):
        xk.append(float(tmp[i]))

    # Metodo
    for k in range(iterMax):
        for j in range(len(vrbs)):
            # Encontrar minimo.
            xk[j] = fmin_preimg(fs, j, vrbs, xk)   

        # Calculo del Error
        gn = eval_gradiente(gs, vrbs, xk) # Gradiente numerico.
        error = norm(gn)                  # Error: Norma del gradiente.
        err.append(error)                 # Agregar el error.
        it.append(k)                      # Agregar la iteracion.

        # Condicion Parada.
        if error <= tol:
            break

    graph(it, err)
    return (xk, error, len(it))

#f = '(x-2)^2 + (y+3)^2 + (x*y)'
f = 'x^3 + y^3 + z^3 - (2*x*y) - (2*x*z) - (2*y*z)'

variables = 'x y z'
vals = "1 1 1"
tol = 1e-8
iterMax = 8

valores = coordinado(f, variables, vals, tol, iterMax)
txt= "\n\nAprox: {0}\nError: {1}\nItera: {2}"
print(txt.format(valores[0], valores[1], valores[2]))

