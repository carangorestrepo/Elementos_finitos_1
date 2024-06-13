clear
clc
syms xi eta 
L1 = 1 - xi - eta;
L2 = xi;
L3 = eta;

Nforma=  simplify([(L1*(9*L1^2 - 9*L1 + 2))/2;
         (L2*(9*L2^2 - 9*L2 + 2))/2;
         (L3*(9*L3^2 - 9*L3 + 2))/2;
         (9*L1*L2*(3*L1 - 1))/2;
         (9*L1*L2*(3*L2 - 1))/2;
         (9*L2*L3*(3*L2 - 1))/2;
         (9*L2*L3*(3*L3 - 1))/2;
         (9*L1*L3*(3*L3 - 1))/2;
         (9*L1*L3*(3*L1 - 1))/2;
         27*L1*L2*L3])


dN_dxi = simplify(diff(Nforma,xi))
     
dN_deta = simplify(diff(Nforma,eta))


Nforma =[ -((eta + xi - 1)*(9*eta + 9*xi + 9*(eta + xi - 1)^2 - 7))/2;
                                  (xi*(9*xi^2 - 9*xi + 2))/2;
                               (eta*(9*eta^2 - 9*eta + 2))/2;
                (xi*(3*eta + 3*xi - 2)*(9*eta + 9*xi - 9))/2;
                       -(xi*(3*xi - 1)*(9*eta + 9*xi - 9))/2;
                                     (9*eta*xi*(3*xi - 1))/2;
                                    (9*eta*xi*(3*eta - 1))/2;
                     -(eta*(3*eta - 1)*(9*eta + 9*xi - 9))/2;
               (eta*(3*eta + 3*xi - 2)*(9*eta + 9*xi - 9))/2;
                               -eta*xi*(27*eta + 27*xi - 27)];
 
 
dN_dxi =[18*eta + 18*xi - 27*eta*xi - (27*eta^2)/2 - (27*xi^2)/2 - 11/2;
                                          (27*xi^2)/2 - 9*xi + 1;
                                                               0;
 (27*eta^2)/2 + 54*eta*xi - (45*eta)/2 + (81*xi^2)/2 - 45*xi + 9;
               (9*eta)/2 + 36*xi - 27*eta*xi - (81*xi^2)/2 - 9/2;
                                            (9*eta*(6*xi - 1))/2;
                                           (9*eta*(3*eta - 1))/2;
                                          -(9*eta*(3*eta - 1))/2;
                                    (9*eta*(6*eta + 6*xi - 5))/2;
                                        -27*eta*(eta + 2*xi - 1)];
 
 
dN_deta =[18*eta + 18*xi - 27*eta*xi - (27*eta^2)/2 - (27*xi^2)/2 - 11/2;
                                                               0;
                                        (27*eta^2)/2 - 9*eta + 1;
                                     (9*xi*(6*eta + 6*xi - 5))/2;
                                            -(9*xi*(3*xi - 1))/2;
                                             (9*xi*(3*xi - 1))/2;
                                            (9*xi*(6*eta - 1))/2;
              36*eta + (9*xi)/2 - 27*eta*xi - (81*eta^2)/2 - 9/2;
 (81*eta^2)/2 + 54*eta*xi - 45*eta + (27*xi^2)/2 - (45*xi)/2 + 9;
                                         -27*xi*(2*eta + xi - 1)];
 

Nforma =[ -((eta  - 1)*(9*eta + 9*xi + 9*(eta + xi - 1)^2 - 7))/2;
                                  (xi*(9*xi^2 - 9*xi + 2))/2;
                               (eta*(9*eta^2 - 9*eta + 2))/2;
                (xi*(3*eta + 3*xi - 2)*(9*eta + 9*xi - 9))/2;


                                     
                                     
Nforma = @(xi,eta) [ ...
    (eta + xi - 1)*(2*eta + 2*xi - 1);
                     xi*(2*xi - 1);
                   eta*(2*eta - 1);
            -xi*(4*eta + 4*xi - 4);
                          4*eta*xi;
           -eta*(4*eta + 4*xi - 4)];
 
 
dN_dxi = @(xi,eta) [ ...
 4*eta + 4*xi - 3;
         4*xi - 1;
                0;
 4 - 8*xi - 4*eta;
            4*eta;
           -4*eta];
 
 
dN_deta = @(xi,eta) [ ... 
 4*eta + 4*xi - 3;
                0;
        4*eta - 1;
            -4*xi;
             4*xi;
 4 - 4*xi - 8*eta];