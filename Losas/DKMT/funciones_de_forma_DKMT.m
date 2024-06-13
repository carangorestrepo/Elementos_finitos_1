syms x21 y21 x32 y32 x31 y31 eta_gl xi_gl phi4 phi5 phi6

L12=(x21^2+y21^2)^(1/2);
L23=(x32^2+y32^2)^(1/2);
L31=(x31^2+y31^2)^(1/2);

C12=x21/L12; 
C23=x32/L23;
C31=x31/L31;

S12=y21/L12;
S23=y32/L23;
S31=y31/L31;

A1=C12*S31-C31*S12;
A2=C23*S12-C12*S23;
A3=C31*S23-C23*S31;

landa=1-eta_gl-xi_gl;
A_delta=-2/3*[(1+phi4),0       ,0;
                     0,(1+phi5),0;
                     0,       0,(1+phi6)];
A_phi=-2/3*[(phi4),0     ,0;
                 0,(phi5),0;
                 0,     0,(phi6)];                 
Au=1/2*[-2/L12,C12,S12, 2/L12,C12,S12, 0    ,0  ,0;
         0,  0,  0,-2/L23,C23,S23, 2/L23,C23,S23;
     2/L31,C31,S31,     0,  0,  0,-2/L31,C31,S31];
 
 Bsgama=[(S31*landa/A1-S23*xi_gl/A2),(S12*xi_gl/A2-S31*eta_gl/A3),(S23*xi_gl/A3-S12*landa/A1);
         (C23*xi_gl/A2-C31*landa/A1),(C31*eta_gl/A3-C12*xi_gl/A2),(C12*landa/A1-C23*eta_gl/A3)];
 
         
Bs=simplify(Bsgama*A_phi*A_delta^(-1)*Au);       