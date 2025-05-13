%Funciones de forma serendipitas del elemento rectangular de 8 nodos:
Nforma=@(xi,eta)[-((eta - 1)*(xi - 1)*(eta + xi + 1))/4;
                 1 - xi^2;
                 ((eta - 1)*(xi + 1)*(eta - xi + 1))/4;
                 1 - eta^2;
                 ((eta + 1)*(xi + 1)*(eta + xi - 1))/4;
                 1 - xi^2;
                 ((eta + 1)*(xi - 1)*(xi - eta + 1))/4;
                 1 - eta^2];
%Derivadas con respecto a xi:
dN_dxi=@(xi,eta)[-((eta + 2*xi)*(eta - 1))/4;
                -2*xi;
                 ((eta - 2*xi)*(eta - 1))/4;
                 0;
                ((eta + 2*xi)*(eta + 1))/4;
                 -2*xi;
                 -((eta - 2*xi)*(eta + 1))/4;
                  0];
%Derivadas con respecto a eta:
dN_deta=@(xi,eta)[-((2*eta + xi)*(xi - 1))/4;
                 0;
                 ((xi + 1)*(2*eta - xi))/4;
                 -2*eta;
                 ((2*eta + xi)*(xi + 1))/4;
                 0;
                 -((xi - 1)*(2*eta - xi))/4;
                 -2*eta];
NN = @(xi)[...
          xi*(xi-1)/2;      % N1
          -2*xi;
          xi*(1+xi)/2 ];   % N3
dNN_dxi = @(xi) [ ...
                   xi - 1/2         % dN1_dxi
                   -2*xi            % dN2_dxi
                   xi + 1/2 ];      % dN3_dxi