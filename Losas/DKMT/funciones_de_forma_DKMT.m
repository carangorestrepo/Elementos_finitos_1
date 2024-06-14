syms x21 y21 x32 y32 x31 y31 eta_gl xi_gl phi4 phi5 phi6

L12=(x21^2+y21^2)^(1/2);
L23=(x32^2+y32^2)^(1/2);
L31=(x31^2+y31^2)^(1/2);

C4=x21/L12; 
C5=x32/L23;
C6=x31/L31;

S4=y21/L12;
S5=y32/L23;
S6=y31/L31;

A1=C4*S6-C6*S4;
A2=C5*S4-C4*S5;
A3=C6*S5-C5*S6;

landa=1-eta_gl-xi_gl;
A_delta=-2/3*[(1+phi4),0       ,0;
                     0,(1+phi5),0;
                     0,       0,(1+phi6)];
A_phi=-2/3*[(phi4),0     ,0;
                 0,(phi5),0;
                 0,     0,(phi6)];                 
Au=1/2*[-2/L12,C4,S4, 2/L12,C4,S4, 0    ,0  ,0;
         0,  0,  0,-2/L23,C5,S5, 2/L23,C5,S5;
     2/L31,C6,S6,     0,  0,  0,-2/L31,C6,S6];
 
Bsgama=[(S6*landa/A1-S5*xi_gl/A2),(S4*xi_gl/A2-S6*eta_gl/A3),(S5*xi_gl/A3-S4*landa/A1);
         (C5*xi_gl/A2-C6*landa/A1),(C6*eta_gl/A3-C4*xi_gl/A2),(C4*landa/A1-C5*eta_gl/A3)];
 
         
Bs=simplify(Bsgama*A_phi*A_delta^(-1)*Au);       