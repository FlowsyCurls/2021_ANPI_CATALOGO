#define ARMA_USE_LAPACK
#include <armadillo>
#include <iostream>
#include <tuple>
#include "matplotlibcpp.h"

using namespace std;
using namespace arma;
namespace plt = matplotlibcpp;

vector<double> verror;
tuple<double, mat, int, double> potencia_inversa(mat A, mat x0, int iterMax, double tol){
    /*
    Esta funcion utiliza el metodo de Potencia Inversa para calcular el
    valor propio de menor magnitud una matrix de entrada 
    simetrica definida positiva asi como el vector propio 
    respectivo a este. 
    
    Sintaxis:   potencia_inversa(A, x0, iterMax, tol)
    
    Parametros Iniciales: 
        A :  matriz simetrica definida positiva a la cual se le deben 
        calcular el vector propio.
        x0 : vector inicial aleatorio no nulo.
        tol = un numero positivo que representa a la tolerancia 
        para el criterio |xkn - xk|<tol
        iterMax = cantidad de iteraciones maximas
        
    Parametros de Salida: 
        valor : valor propio de menor magnitud de la matriz.
        vector : vector propio respectivo al valor propio calculado.
        error : el error obtenido bajo el criterio de |xkn - xk|.
    */
   
   
    int n = A.n_rows;
    mat yk, xk = x0, xk_n;
    double ck, error;

    int i = 0;
    for (i; i < iterMax; i++) // Iteraraciones
    {
        yk = arma::solve(A, xk);
        ck = norm(yk, "inf");
        xk_n = (1 / ck) * yk;
        error = norm(xk_n - xk);
        verror.push_back(error);
        xk = xk_n;

        if (error < tol)
            break;
    }
    return {ck, xk.t(), i, error};
}

int main()
{

    // Entradas.
    mat A = {{3, -1, 0}, 
             {-1, 2, -1}, 
             {0, -1, 3}};

    mat x0 = {1, 1, 1};
    x0=x0.t();

    int iterMax = 9;
    double tol = 10e-10;


    // Salidas
    mat vector;
    int k;
    double valor, err;
    tie(valor, vector, k, err) = potencia_inversa(A, x0, iterMax, tol);


    cout << "\n • valor propio  : " << valor
         << "\n • vector propio :" << vector
           << " • iteraciones   : " << k
           << "\n • error aprox.  : " << err
         << endl << endl;


    plt::plot(verror);
    plt::title("Metodo Potencia Inversa (iteraciones vs error)");
    plt::xlabel("Iteraciones k");
    plt::ylabel("Error");
    plt::show();

    return 0;
}