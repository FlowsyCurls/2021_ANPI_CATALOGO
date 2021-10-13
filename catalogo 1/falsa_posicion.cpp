#include <cmath>
#include <iostream>
#include <ginac/ginac.h>
//#include "matplotlibcpp.h"
using namespace std;
using namespace GiNaC;
//namespace plt = matplotlibcpp;

/** Metodo encargado de evaluar una funcion representada en string en un valor
* especifico utilizando la libreria GiNaC
* @param funct La funcion que se quiere evaluar en formato de string
* @param value El valor en el cual se va a evaluar la funcion
*/
ex f(string funct, ex value){
  symbol x;
  symtab table;
  table["x"] = x;
  parser reader(table);
  ex f = reader(funct);
  return evalf(f.subs(x==value));
}

/** Metodo encargado de calcular una solucion aproximada a una funcion utilizando el
* metodo de falsa posicion
* @param funct String de la funcion a la cual se le va a calcular el cero aproximado
* @param startInt El inicio del intervalo
* @param endInt El final del intervalo
* @param maxIt Las iteraciones maximas que va a realizar la funcion
* @param tol La tolerancia maxima permitida
*/
ex * falsa_posicion(string funct,ex startInt, ex endInt, int maxIt, ex tol){
  ex a = startInt;
  ex b = endInt;
  ex xk = a-((a-b)/(f(funct, a)-f(funct, b)))*f(funct, a);
  ex functXk;
  vector<double> graphX(maxIt), graphY(maxIt);
  int i=0;
  static ex r[2];
  ex err = tol+1;
  while(i<maxIt && err>tol){
    functXk = f(funct, xk);
    if(f(funct, a)*functXk<0){
      xk = xk-functXk*(xk-a)/(functXk-f(funct, a));
      // Para la siguiente iteracion
      b = xk;
    }else if(functXk*f(funct, b)<0){
      xk = xk-functXk*(xk-b)/(functXk-f(funct, b));
      // Para la siguiente iteracion
      a = xk;
    }else{
      cout << "No se puede resolver por falsa posicion" << endl;
      return 0;
    }
    // Grafica
    graphX.push_back(i);
    graphY.push_back(ex_to<numeric>(abs(f(funct,xk))).to_double());
    i+=1;
  }
  //plt::plot(graphX, graphY);
  //plt::show();
  r[0] = evalf(xk);
  r[1] = abs(f(funct, xk));
  return r;
}

int main(void){
  ex *r;
  r = falsa_posicion("cos(x)-x", 1/2, Pi/4, 4, 1e-10); 
  cout << "Aproximacion: " << *r << endl;
  cout << "Error: " << *(r+1) << endl;
}
