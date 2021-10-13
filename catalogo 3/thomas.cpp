#include <iostream>
#include<armadillo>

using namespace std;
using namespace arma;



// Input: Matriz A de tamano mxn.
// Output: true si la matriz corresponde a una cuadrada, 
//         false en caso contrario.
bool is_cuadrada(mat A){
    /// Funcion que verifica si una matriz de entrada es cuadrada.
    return A.n_rows == A.n_cols;
}

// Input: Matriz A de tamano mxn.
// Output: true si la matriz corresponde a una tridiagonal, 
//         false en caso contrario.
int is_tridiagonal(mat A)
{
    /// Funcion que verifica si una matriz de entrada es tridiagonal.  
    
    // Es matriz cuadrada?
    if (is_cuadrada(A) != true) return false;

    int m = A.n_rows; 
    for(int i=0; i<m; i++) // Iteraracion en las filas.
    {        
        for (int j = 0; j < m; j++) // Iteracion en las columnas.
        {                             
            if (j>i+1 or j<i-1) // Si la parte sup o inf no es 0.
                if (A(i,j)!=0) 
                    return false;
        }
    }
    return true;
}


mat thomas(mat A, mat b)
{
    /*
        Esta funcion encuentra las soluciones de un sistema de
        ecuaciones utilizando el metodo de Thomas.

        "En cada iteracion atraves de las filas de la matriz, el 
        algoritmo va calculando valores pi y qi y rellenando dos 
        vectores p, de tamano m-1 y q, de tamano m, a traves de una
        formula dependiendo de la fila en la que se esta iterando.
        
        Finalmente da la solucion iterando desde m hacia atras y 
        asignando a las posiciones de un vector x las soluciones,
        utilizando otra formula dependiendo de la fila y los 
        vector p y q. "

    Sintaxis:  thomas(A,b);
    
    Parametros Iniciales: 
        A: Matriz A de tamano m (cuadrada, triangular superior e invertible).
        b: Vector b de tamano m.

    Parametros Salida: 
        x : Vector x de tamano m, que es solucion del sistema Ax=b.

    */

    // Es matriz tridiagonal?
    if (is_tridiagonal(A) != true) {
        cout << "Matriz no es tridiagonal." << endl;
        return NULL;
    }
    
    int m=A.n_rows;
    double aux;
    b = b.t();
    mat x = zeros<mat>(m, 1), p = zeros<mat>(m - 1, 1), q = zeros<mat>(m, 1);

    // Obtener los vectores p y q.
    for (int i = 0; i < m; i++) // Iteracion en las filas.
    {
        if (i==0) {    // si i es la primera columna;
            p(i) = A(i, i+1) / A(i, i);  // ci / bi
            q(i) = b(i) / A(i, i);      // di / bi
        }
        else {   
            aux = A(i, i) - p(i-1) * A(i, i - 1);         // (bi-p(i-1)*ai)
            q(i) = (b(i) - q(i-1) * A(i, i - 1)) / aux;  // (di - q(i-1)*ai) / aux

            if ( i<m-1 ){
                p(i) = A(i, i + 1) / aux;                     // ci / aux
            }
        }
    }

    int n = m-1;
    x(n) = q(n);            // xn = qn
    for (int i = n-1; i>=0; i--)  {  // Iteracion hacia atras.
        x(i) = q(i) - p(i) * x(i + 1);
    }
    return x;
}


int main()
{
    mat A={
        {-4,1,0,0,0},
        {1,-4,1,0,0},
        {0,1,-4,1,0},
        {0,0,1,-4,1},
        {0,0,0,1,-4}};
    mat b = {1,1,1,1,1};

    cout << "▶ A\n"
         << A << endl;
    cout << "▶ b\n"
         << b.t() << endl;

    mat x = thomas(A, b);
    cout << "▶ Soluciones del sistema: \n • x:\n" << x << endl;    
    return 0;
}