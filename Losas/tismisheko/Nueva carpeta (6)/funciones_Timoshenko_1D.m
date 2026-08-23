function [Nv,Nt,dNv,dNt] = funciones_Timoshenko_1D(x,L,EI,Ks)
%==========================================================================
% FUNCIONES EXACTAS DE VIGA TIMOSHENKO SIMPLE
% Estado: [v t M V]
% Ecuaciones: V'=0, M'=V, t'=M/EI, v'=t-V/Ks
% GDL: [v1 t1 v2 t2]
%==========================================================================

A = [1, x, x^2/(2*EI), x^3/(6*EI)-x/Ks;
     0, 1, x/EI,       x^2/(2*EI);
     0, 0, 1,          x;
     0, 0, 0,          1];

dA = [0, 1, x/EI, x^2/(2*EI)-1/Ks;
      0, 0, 1/EI, x/EI;
      0, 0, 0,    1;
      0, 0, 0,    0];

D = [1 0 0 0;
     0 1 0 0];

A0 = eye(4);
AL = [1, L, L^2/(2*EI), L^3/(6*EI)-L/Ks;
      0, 1, L/EI,       L^2/(2*EI);
      0, 0, 1,          L;
      0, 0, 0,          1];

B = [D*A0;
     D*AL];

H = A/B;
dH = dA/B;

Nv = H(1,:);
Nt = H(2,:);
dNv = dH(1,:);
dNt = dH(2,:);
end
