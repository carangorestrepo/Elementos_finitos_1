syms xi eta
Nforma = [( 1.0 - 9.0 * xi * eta ) * ( 1.0 - xi - eta );
           xi * ( 1.0 - 9.0 * ( 1.0 - xi - eta ) * eta );
           eta * ( 1.0 - 9.0 * ( 1.0 - xi - eta ) * xi );
           27.0 * ( 1.0 - xi - eta ) * xi * eta];
    
dN_dxi=diff(Nforma,xi)
dN_deta=diff(Nforma,eta)

dN_dxi =
 
     9*eta*(eta + xi - 1) + 9*eta*xi - 1
   9*eta*xi + eta*(9*eta + 9*xi - 9) + 1
                 eta*(9*eta + 18*xi - 9)
 - 27*eta*xi - eta*(27*eta + 27*xi - 27)
 
 
dN_deta =
 
     9*xi*(eta + xi - 1) + 9*eta*xi - 1
                 xi*(18*eta + 9*xi - 9)
   9*eta*xi + xi*(9*eta + 9*xi - 9) + 1
 - 27*eta*xi - xi*(27*eta + 27*xi - 27)