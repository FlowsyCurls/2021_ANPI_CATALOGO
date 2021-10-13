#include <iostream>
#include <armadillo>
#include "matplotlibcpp.h"

using namespace std;
using namespace arma;
namespace plt = matplotlibcpp;

// vector para graficar
vector<double> error;
tuple<mat, int, double> newton_schultz(mat A, mat b, double tol, int iterMax)
{
    /*
    Metodo de Newton_Schultz para aproximar la pseudoinversa A+ de A.
    */
    // Declaracion de variables
    mat Xk, xk, xk_n;
    double err = tol + 1;

    Xk = (1 / pow(norm(A, 2), 2)) * A.t();  // Valor inicial.
    xk = Xk * b; // Solucion con el valor actual.

    int k = 1, m = A.n_rows;  
    mat I = eye(m, m); // Matriz identidad I - # de filas de A.

    for (k; k < iterMax; k++) // Iteracion maximas.
    {
        Xk = Xk * (2 * I - A * Xk);
        xk_n = Xk * b; // Solucion con el valor actual.
        err = norm(xk_n - xk) / norm(xk_n); // Error de sucesion de elementos.
        xk = xk_n;
        error.push_back(err);  

        // Condicion de Parada
        if (err < tol) break;
    }

    return {xk, k, err};
}


tuple<mat, int, double> pseudoinversa(mat A, mat b, double tol, int iterMax)
{
    /*
        Descripcion:
            Metodo de la pseudoinversa A+.
              Para encontrar el vector x que mejor aproxima la solucion Ax=b 
            de un sistema de ecuaciones puede tener infinitas soluciones o 
            ninguna solucion.
              Se utiliza el algoritmo iterativo para realizar el calculo de la 
            pseudoinversa conocido como metodo de Schulz-Ben-Israel, el cual se
            basa en el metodo de  Newton-Rapshon visto en las clases anteriores
            para aproximar el valor de A+.

        Parametros Iniciales: 
            A       : Matriz A de coeficientes.
            b       : Vector b de terminos independientes.
            tol     : tolerancia al error.
            iterMax : cantidad de iteraciones maximas.

        Parametros Salida: 

            xk  : aproximacion de la solucion que mejor cumple Ax=b.
            k   : cantidad de iteraciones.
            err : error relacionado a la aproximacion.
    */
    mat xk; int k; double err;
    tie(xk, k, err) = newton_schultz(A, b, tol, iterMax);

    // Error relacionado a la aproximacion.
    double err_xk = norm(A * xk - b);
    return {xk, k, err_xk};
}




int main()
{

    // Entradas.
    mat A = {{1, 2, 3}, {1, 8, 9}, {1, 4, 1}, {44, 0, 1}, {1, 4, 5}};
    mat b = {1, 1, 1, 1, 1};
    b=b.t();
    
    // mat b = tmp.t();
    double tol = 10e-5;
    int iterMax = 1000;

    // Salidas
    mat xk;
    int k;
    double err;
    tie(xk, k, err) = pseudoinversa(A, b, tol, iterMax);

    cout << "▶ A\n"
         << A << endl;
    cout << "▶ b\n"
         << b << endl;
    cout << "▶ Solucion, iteraciones, error\n"
         << "\n • xk  : \n"
         << xk
         << " • k   : " << k
         << "\n • err : " << err
         << "\n"
         << endl;

    plt::plot(error);
    plt::title("Metodo de la pseudoinversa (iteraciones vs error)");
    plt::xlabel("Iteraciones k");
    plt::ylabel("Error");
    plt::show();

    return 0;
}