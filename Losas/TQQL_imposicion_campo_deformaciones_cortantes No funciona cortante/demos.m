
syms S4 S5 S6 
syms C4 C5 C6  
syms Bx1 Bx2 Bx3 Bx4
syms By1 By2 By3 By4
syms dphi4 dphi5 dphi6 

L12xi=1;
L13xi=1;
L23xi=(2)^(1/2);

C4 = x21/L4;      S4 = y21/L4;
C5 = x32/L5;      S5 = y32/L5;
C6 = x13/L6;      S6 = y13/L6;

Un1 = [ w1 Bx1 By1 w2 Bx2 By2 w3 Bx3 By3].';

g12xi=(w1 - w2)/L4 -1/2*((Bx1*C4)+(By1*S4)+(By2*S4)+(Bx2*C4))-dphi4*2/3*L4/L12xi;
g23xi=(w2 - w3)/L5 -1/2*((Bx2*C5)+(By2*S5)+(By3*S5)+(Bx3*C5))-dphi5*2/3*L5/L23xi;
g31xi=(w3 - w1)/L6 -1/2*((Bx3*C6)+(By3*S6)+(By1*S6)+(Bx1*C6))-dphi6*2/3*L6/L13xi;


Au = equationsToMatrix([ g12xi; g23xi; g31xi ], Un);