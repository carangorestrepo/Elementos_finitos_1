clc
clear
syms xi eta 
L1 = 1 - xi - eta;
L2 = xi;
L3 = eta;

Nforma = [L1;
         L2;
         L3];
       
dN_dxi = simplify(diff(Nforma,xi));
     
dN_deta = simplify(diff(Nforma,eta));