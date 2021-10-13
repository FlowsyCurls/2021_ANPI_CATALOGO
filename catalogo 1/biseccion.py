from sympy import  sympify
import numpy as np
import matplotlib.pyplot as plt

def biseccion(f,a,b,tol):

    """
    Esta funcion aproxima la solucion de la ecuacion f(x)=0,
    utilizando el metodo de la biseccion
    
    Sintaxis:  biseccion(f,a,b,tol)
    
    Parametros Iniciales: 
        f : una cadena de caracteres (string) que representa a la funcion f
        a,b : son los extremos del intervalo [a,b]
        tol : un numero que representa la tolerancia para el criterio |f(x_k)|<tol
                
    Parametros Salida: 
        [x, k] donde:                
        x : aproximacion del cero de la funcion f
        k : numero de iteraciones realizados
        error :  |f(x_k)|
    """
        
    funcion = sympify(f)
    fa=funcion.subs({'x': a})
    fb=funcion.subs({'x': b})
    k= 0;
    error= 0;
    if fa*fb>0:
        print('El teorema de Bolzano no se cumple, es decir, f(a)*f(b)>0')
    else:
        error=tol+1
        k=0
        it=[]
        er=[]
        while error>tol:
            k=k+1
            x=(a+b)/2
            fa=funcion.subs({'x': a})
            fx=funcion.subs({'x': x})
            error=abs(fx)
            it.append(k)
            er.append(error)
            if fa*fx<0:
                b=x                
            else:
                a=x   

    fig, graf = plt.subplots()
    graf.plot(it, er)
    graf.set_xlabel('Iteraciones ($k$)')
    graf.set_ylabel('Error $|f(x_k)|$')
    graf.set_title('Metodo de la Biseccion (Iteraciones vs Error)');
    graf.grid(True)
    plt.show()

    return [x,error]

f = 'exp(x) - x - 2'
a = 0
b = 2
tol = 0.0001
resultado = biseccion(f,a, b, tol)
print("xk:")
print(resultado[0])
print("error:")
print(resultado[1])

