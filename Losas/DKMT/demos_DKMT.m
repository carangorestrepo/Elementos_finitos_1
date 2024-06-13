clc
clear
syms w1 w2 w3 w4 L5 L6 L7 L8 c5 c6 c7 c8 s5 s6 s7 s8 bx1 bx2 bx3 bx4 by1 by2 by3 by4 db5 db6 db7 db8 phi5 phi6 phi7 phi8 dBs5 dBs6 dBs7 dBs8
syms x43 x32 x21 x14 y43 y32 y21 y14 v h xi eta
c5 = x21/L5;      s5 = y21/L5;
c6 = x32/L6;      s6 = y32/L6;
c7 = x43/L7;      s7 = y43/L7;
c8 = x14/L8;      s8 = y14/L8;
dbs51 = (w2 - w1 + L5*(c5*bx1 + s5*by1)/2 + L5*(c5*bx2 + s5*by2)/2)/(-2*L5*(1+phi5)/3);
dbs61 = (w3 - w2 + L6*(c6*bx2 + s6*by2)/2 + L6*(c6*bx3 + s6*by3)/2)/(-2*L6*(1+phi6)/3);
dbs71 = (w4 - w3 + L7*(c7*bx3 + s7*by3)/2 + L7*(c7*bx4 + s7*by4)/2)/(-2*L7*(1+phi7)/3);
dbs81 = (w1 - w4 + L8*(c8*bx4 + s8*by4)/2 + L8*(c8*bx1 + s8*by1)/2)/(-2*L8*(1+phi8)/3);

dbs5 = -(w2 - w1 + L5*(c5*bx1 + s5*by1)/2 + L5*(c5*bx2 + s5*by2)/2);
dbs6 = -(w3 - w2 + L6*(c6*bx2 + s6*by2)/2 + L6*(c6*bx3 + s6*by3)/2);
dbs7 = -(w4 - w3 + L7*(c7*bx3 + s7*by3)/2 + L7*(c7*bx4 + s7*by4)/2);
dbs8 = -(w1 - w4 + L8*(c8*bx4 + s8*by4)/2 + L8*(c8*bx1 + s8*by1)/2);

Un = [ w1 bx1 by1 w2 bx2 by2 w3 bx3 by3 w4 bx4 by4 ].';

An = equationsToMatrix([ dbs5; dbs6; dbs7; dbs8 ], Un);


An1 = equationsToMatrix([ dbs51; dbs61; dbs71; dbs81 ], Un);

Adb = diag((2/3)*[ L5*(1+phi5) L6*(1+phi6) L7*(1+phi7) L8*(1+phi8) ]);

An=Adb^-1*An

N5=1/2*(1-eta);
N6=1/2*(1+xi);
N7=1/2*(1+eta);
N8=1/2*(1-xi);


syms xi eta x1 x2 x3 x4 y3 y1 y2 y3 y4
N1 = (1 - xi  )*(1 - eta  )/4;
N2 = (1 + xi  )*(1 - eta  )/4;
N3 = (1 + xi  )*(1 + eta  )/4;
N4 = (1 - xi  )*(1 + eta  )/4;

N5 = (1 - xi^2)*(1 - eta  )/2;
N6 = (1 + xi  )*(1 - eta^2)/2;
N7 = (1 - xi^2)*(1 + eta  )/2;
N8 = (1 - xi  )*(1 - eta^2)/2;

%% Jacobian and inverse Jacobian
% isoparametric interpolation
x = N1*x1 + N2*x2 + N3*x3 + N4*x4;
y = N1*y1 + N2*y2 + N3*y3 + N4*y4;

J = simplify([ diff(x,xi)   diff(y,xi)
               diff(x,eta)  diff(y,eta) ]);
  
gsz5=-2/3*(2/(5/6*(1-v)))*(h^2/L5^2)*dbs5;
gsz6=-2/3*(2/(5/6*(1-v)))*(h^2/L6^2)*dbs6;
gsz7=-2/3*(2/(5/6*(1-v)))*(h^2/L7^2)*dbs7;
gsz8=-2/3*(2/(5/6*(1-v)))*(h^2/L8^2)*dbs8;

  
ec = equationsToMatrix([ gsz5+ gsz7; gsz6+ gsz8 ], Un);
  
gbar=simplify(J*ec)

