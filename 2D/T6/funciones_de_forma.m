syms xi eta 
L1 = 1 - xi - eta;
L2 = xi;
L3 = eta;

N= [L1*(2*L1 - 1);
 L2*(2*L2 - 1);
 L3*(2*L3 - 1);
 4*L1*L2;
 4*L2*L3;
 4*L1*L3];

dNxi=diff(N,xi,1);
dNeta=diff(N,eta,1);


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