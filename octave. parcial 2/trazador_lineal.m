function ejemplo()
     clc; clear;
     pkg load symbolic;

     warning('off','all');
     xv = [ 1, 2, 4];
     yv = [1, 0.5, 0.25];
     s=traz_lineal(xv, yv);
  
end

function s= traz_lineal(xv, yv)
     syms x;
     
     n=length(xv);
     for i=1:(n-1)
          x_curr = xv(i);        
          x_next = xv(i+1);
          
          y_curr = yv(i);        
          y_next = yv(i+1);
      
          s(i) =  (y_curr+((y_next-y_curr)/(x_next-x_curr))*(x-x_curr));

          if i==1
               display("S(0)\n");
               disp(s(i));
               continue;
          end
          printf ("\nS(%d)\n\n", i-1);
          disp(s(i));
     end
    
end
  