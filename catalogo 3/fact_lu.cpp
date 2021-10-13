#include <iostream>
#include<armadillo>
#include <tuple>

using namespace std;
using namespace arma;

// Input: Matriz A de tamano mxn.
// Output: true si la matriz corresponde a una cuadrada, 
//         false en caso contrario.
bool is_cuadrada(mat A){
    /// Función que verifica si una matriz de entrada es cuadrada.
    return A.n_rows == A.n_cols;
}

/// Función Sustitución Hacia Atrás 
// Input: Matriz A (cuadrada, triangular superior e invertible).
// Output: Vector x de tamano m, que es solucion del sistema Ax=b.
vec sust_atras(mat A, mat b){
  int m=A.n_rows;
  mat x=zeros<mat>(m,1);
  for (int i = m-1; i>=0;i--){
    double aux=0;
    for (int j=i+1; j<=m-1;j++){
      aux+=A(i,j)*x(j);
    }
    x(i)=(1/A(i,i))*(b(i)-aux);
  }
  return x;
}

/// Función Sustitución Hacia Adelante 
// Input: Matriz A (cuadrada, triangular superior e invertible).
// Output: Vector x de tamano m, que es solucion del sistema Ax=b.
mat sust_adelante(mat A, mat b){
  int m=A.n_rows;
  mat x=zeros<mat>(m,1);
  for (int i=0; i<m; i++)
  {
      double aux = 0;
      for (int j=0; j<i; j++)
      {
          aux += A(i, j) * x(j);
      }
      x(i) = (1 / A(i, i)) * (b(i) - aux);
  }
  return x;
}


tuple<mat,mat> lu(mat A){
    /*
        Esta funcion encuentra la matriz L y la matriz U de la
        factorizacion LU de la matriz A.

        "Para la obtencion de la matriz U se realiza la operacion
        de restar y/o sumar el multiplo de otra fila en la matriz A,
        convirtiendola en una matriz triangular superior."

        "Para la obtencion de la matriz L se guardan los multiplicadores 
        (con signo opuesto) utilizados para convertir en cero las entradas
        debajo de la diagonal principal."


    Sintaxis:  lu(A);
    
    Parametros Iniciales: 
        A: Matriz A de tamano m (cuadrada e invertible).

    Parametros Salida: 
        (L,U) : Tupla con matriz L y matriz U de la factorizacion LU
                de la matriz A.

    */

    // Es matriz cuadrada?
    if (is_cuadrada(A) != true) {
        cout << "No es matriz cuadrada.\n" << endl;
        return {{0}, {0}};
    }

    int m=A.n_rows;
    for (int i=0; i<m; i++) {

        // Obtener submatriz de tamano i.
        mat submat = A.submat(0, 0, i, i);

        // Es invertible?
        if(det(submat) == 0){
            cout << "Una submatriz no es invertible.\n"  
                << "No existe una unica solucion.\n" << endl;
            return {{0}, {0}};
        }
    }

    mat L = eye(m, m); // Matriz identidad L.
    for (int k = 0; k < m; k++) {
        // Elementos por debajo de la diagonal.
        for (int i = k+1; i < m; i++) {
            L(i, k) = A(i, k) / A(k, k); // Obtener factor multiplicador.
            A(i, k) = 0;                 // Hacer cero.

            // Multiplicar los elementos restantes de la fila.
            for (int j = k+1; j < m; j++) {
                A(i, j) = A(i, j) - L(i, k) * A(k, j);
            }
        }
    }
    return {L, A};
}


mat fact_lu(mat A, mat b)
{
    /*
        Esta funcion utiliza factorizacion LU de la matriz A, obtenida
        llamando al metodo lu(A), para encontrar la solucion al sistema
        de ecuaciones lineal  Ax = b, con x como vector solucion.

        "Pasos de la factorizacion: Ax = b
            Se sabe que A = LU, por lo tanto LUx = b;

            Tomando y=Ux, se resuelve el sistema Ly=b con 'y' como vector 
            solucion, utilizando la sustitucion hacia adelante (pues L es
            triangular inferior).

            Una vez encontrados los valores de 'y', se resuelve el sistema 
            Ux=y con 'x' como vector solucion, utilizando sustitucion hacia 
            atras (pues U es triangular superior.)

    Sintaxis:  fact_lu(A,b);
    
    Parametros Iniciales: 
        A: Matriz A de tamano m (cuadrada e invertible).
        b: Vector b de tamano m.

    Parametros Salida: 
        x : Vector x de tamano m, que es solucion del sistema Ax=b.

    */

    mat L, U;
    tie(L, U) = lu(A); // Tupla con L y U.
    cout << "▶ L\n" << L << endl;
    cout << "▶ U\n" << U << endl;


    // // Si no se logra la factorizacion.
    if ((L(0) == 0) & U(0) == 0)
        return NULL;

    mat y = sust_adelante(L, b.t());   // Ly = b
    mat x = sust_atras(U, y.t());      // Ux = y
    return x;
}

int main()
{
    mat A = {
        {2, 3, 0, 1},
        {4, 5, 3, 3},
        {-2, -6, 7, 7},
        {8, 9, 5, 21}};
    mat b = {1, 1, 1, 1};

    cout << "▶ A\n"
         << A << endl;
    cout << "▶ b\n"
         << b.t() << endl;

    mat x = fact_lu(A,b);
    cout << "▶ Soluciones del sistema: \n • x:\n" << x << endl;
    return 0;
}