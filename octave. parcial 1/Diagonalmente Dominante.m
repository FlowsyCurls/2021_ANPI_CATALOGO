%Comparar si la matriz es diag. domin.clc; clear;A=[10 -2 1 1;   1 -20 1 1;   -1 3 -3 -8;   0 4 12 -6]m = rows(A);for f=1:m     %Valor de la diagonal de dicha fila     d=A(f,f)     %Otros valores de la fila, sin incluir la diagonal     f1=A(f,:);     f1(f)=0;     s = sum(abs(f1))     if abs(d)<=sum(abs(f1))          display('No es diagonalmente dominante.')          display('')     else       display('Si es diagonalmente dominante.')        display('')     end  end  