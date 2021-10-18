from sympy import *
from scipy import optimize


# Ejemplo
def ejemplo():
    f = 'sin(x*pi/2)';
    ptos=[-1, 0, 1, 2];  ## Valores de soporte en Intervalo [-1, 2]
    val = 0.54
    cota = cota_poly_inter(f, ptos)
    ## Notacion cientifica.
    # cota = np.format_float_scientific(cota, precision=6)

    string = ("\n ➤  La cota de error del polinomio de interpolación de esta\
        \n    funcion en el intervalo de {0} a {1} es de: \n\n\t {2} \n\n\
        ").format(ptos[0], ptos[-1], cota);
    print(string)

# Cota de Error Polinomio de Interpolación
def cota_poly_inter(f, ptos):
    """
    Esta funcion calcula el error que se genera al realizar una aproximacion 
    con un polinomio de interpolacion en el intervalo [a,b] y utilizando un 
    vector soporte que contiene puntos dentro del intervalo mencionado.
    
    Sintaxis:  cota_poly_inter(f,ptos)
    
    Parametros Iniciales: 
        f : una cadena de caracteres (string) que representa
            a la funcion f
        ptos : es un vector que contiene los puntos de soporte
            dentro de un intervalo determinado.

    Parametros Salida: 
        cota :  corresponde a la cota del error
    """
    n=len(ptos)-1;                  ## Grado maximo del polinomio.
    fs = sympify(f)                 ## Funcion simbolica.

    ## Extremos del intervalo.
    a = ptos[0]
    b = ptos[-1]

    alpha = alpha_max(fs, n+1, a, b)
    p_max = polynomial_max(ptos, a, b)

    # Calculo del valor del polinomio con el valor maximo.
    suma = 1
    for k in ptos:
        suma *=abs(p_max-k)
    suma = simplify(suma)
    cota = alpha / factorial(n+1) * suma         ## Cota
    return cota

def alpha_max(fs, m, a, b):
    """
        Esta funcion calcula el alpha_max para el metodo de
        cota interpolacion.
        Calcula la derivada de d(n+1)f y obtiene el punto maximo.
        Retorna la funcion evaluada en este punto.
    """
    x = Symbol('x')                 ## Variable simbolica.
    dfm_s =  fs.diff(x, m)          ## Derivada m de f simbolica.
    dfm_n = lambdify(x, dfm_s)      ## Derivada m de f numerica.
    ## Calculo del maximo utilizando funcion del minimo sobre -f(x).
    ## min{ -f } en [a,b] -> max{ f } en [a,b].
    fs_aux = -1*abs(dfm_s)
    fn_aux = lambdify(x, fs_aux)
    x_max = optimize.fminbound(fn_aux, a, b)    ## x_max.
    alpha = abs(dfm_n(x_max))                   ## alpha.
    return alpha

def polynomial_max(ptos,a,b):
    """
        Esta funcion calcula el punto maximo del polinomio 
        de la funcion generado en el metodo.
        Retorna la preimagen del punto maximo.
    """
    x = Symbol('x')                 ## Variable simbolica.
    ## Creacion de polinomio simbolico.
    expr = 1
    for i in ptos:
        expr *= (x - i)

    expr_s = expand(expr)
    expr_n = lambdify(x, expr_s)
    deriv_expr_s = diff(expr_s)   ## Derivada simbolica del polinomio.
    
    # Obtener los puntos criticos dentro del intervalo [a, b] y evaluar
    # la expresion en los extremos del intervalo y en los puntos criticos.
    crit = solveset( Eq( deriv_expr_s, 0 ), domain=Interval(a, b))

    Xv = [a,b]
    Yv = [abs(simplify(expr_n(a))), abs(simplify(expr_n(b)))]
    for i in range(len(crit)-1, -1,-1):
        punto = crit.args[i]
        if a < punto < b:
            Xv.append(simplify(punto))
            Yv.append(abs(simplify(expr_n(punto))))

    # Obtener el maximo de los resultados.
    indice_max = Yv.index(max(Yv))
    p_max = Xv[indice_max]
    return p_max

ejemplo()