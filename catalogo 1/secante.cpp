#include <iostream> 
#include <cmath>
#include <ginac/ginac.h>
#include <cln/dfloat_io.h>
#include "matplotlibcpp.h"
using namespace std;
using namespace GiNaC;
using namespace cln;
namespace plt = matplotlibcpp;

/** Metodo que realiza la evaluacion de una funcion en un valor especifico dado 
* utilizando la libreria GiNaC.
* @param funct La funcion a evaluar en formato de string
* @param value El valor en el cual se va a evaluar la funcion
*/
ex f(string funct, ex value){
  symbol x;
  symtab table;
  table["x"]=x;
  parser reader(table);
  ex f = reader(funct);
  return evalf(f.subs(x==value));
}

/** Funcion que calcula la aproximacion de un cero de una funcion dada mediante el
* metodo de secante
* @param x0 First initial value
* @param x1 Second initial value
* @param maxIt Max Iterations for the function
* @param tol Tolerancia de la funcion
* @param funct Funcion a evaluar
*/
ex * secante(double x0, double x1, int maxIt, string funct, ex tol){
  vector<double> graphX(maxIt), graphY(maxIt);
  ex xk = x1;
  ex xkm1 = x0;
  ex xkp1;
  int i=0;
  ex err=tol+1;
  static ex r[2]; // Array de resultado
  while(i<maxIt && err>tol){
    // Nuevo resultado
    xkp1 = xk-f(funct,xk)*(xk-xkm1)/(f(funct,xk)-f(funct,xkm1));
    // Reasignando valores para la siguiente iteracion
    xkm1 = xk;
    xk = xkp1;
    err = abs(f(funct, xk));
    // Para la grafica
    graphX.push_back(ex_to<numeric>(i).to_double());
    graphY.push_back(ex_to<numeric>(abs(f(funct,xk))).to_double());
    // Contador de iteracion
    i+=1;
  }
  plt::plot(graphX,graphY);
  plt::show();
  r[0] = xk; // Aproximacion
  r[1] = abs(f(funct,xk)); // Error
  return r;
}

int main(void)
{
  ex *r;
  r = secante(0.75,0.437193, 5, "(cos(2*x)^(2))-((x)^(2))",1e-10);
  cout << "Aproximacion: "<<*r<<endl;
  cout << "Error: " << *(r+1) << endl;
  return 0;
}
