# -*- coding: utf-8 -*-
"""
Ejemplo Trapecio Compuesto y su Cota de Error

@author: jusoto
"""
from sympy import sympify, Symbol, lambdify, Abs
from numpy import linspace
from scipy.optimize import fminbound

"Definir la función"
x=Symbol("x")
f="log(x)"
fs=sympify(f)

"Aproximación Utilizando el Método del Trapecio Compuesto"
num_punt=800
a=2
b=5
h=(b-a)/(num_punt-1)
xv=linspace(a,b,num_punt)
I=0
for i in range(0,num_punt-1):
    ai=xv[i]
    bi=xv[i+1]
    fai=fs.subs(x,ai)
    fbi=fs.subs(x,bi)
    I+=(bi-ai)*(fai+fbi)/2
print("\nAprox:\t",I)

"Cota del error"

"1. Calculo de la segunda derivada"
fs2d=Abs(fs.diff(x,2))
fn2d=lambdify(x,fs2d)
"2. Calculo de las funciones auxiliares: max f = min -f"
fs2d_aux=-Abs(fs.diff(x,2))
fn2d_aux=lambdify(x,fs2d_aux)
"3. Calculo de alpha_max"
xmax=fminbound(fn2d_aux,a,b)
alpha_max=fn2d(xmax)

"4. Calculo de la cota de error"
cota_error=((b-a)*h**2/12)*alpha_max

print("Cota:\t",cota_error)