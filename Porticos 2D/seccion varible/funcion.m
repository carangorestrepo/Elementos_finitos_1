function [MEE]=funcion(q1,q2,nq,sec,b1,b2,nb,h1,h2,nh,G,E,tf,tw,L,g,v,a,b)

syms x
%V=C1 + q1*x - (x*x^nq*(q1 - q2))/(L^nq + L^nq*nq);
%M=C2 + C1*x + (q1*x^2)/2 - (x^nq*x^2*(q1 - q2))/(2*L^nq + L^nq*nq^2 + 3*L^nq*nq);
V= q1*x - (x*x^nq*(q1 - q2))/(L^nq + L^nq*nq);
M=(q1*x^2)/2 - (x^nq*x^2*(q1 - q2))/(2*L^nq + L^nq*nq^2 + 3*L^nq*nq);

[~,EIx,Acx]=secciones_1(sec,b1,b2,nb,h1,h2,nh,G,E,tf,tw,a,b);

MC2=area_cuadraturas(a,b,matlabFunction(1/EIx,"Vars",{x}));
MC1=area_cuadraturas(a,b,matlabFunction(x/EIx,"Vars",{x}));
Mq=area_cuadraturas(a,b,matlabFunction(M/EIx,"Vars",{x}));
MC2x=area_cuadraturas(a,b,matlabFunction(1/EIx*x,"Vars",{x}));
MC1x=area_cuadraturas(a,b,matlabFunction(x/EIx*x,"Vars",{x}));
Mqx=area_cuadraturas(a,b,matlabFunction(M/EIx*x,"Vars",{x}));
Vqx=area_cuadraturas(a,b,matlabFunction(V/Acx,"Vars",{x}));
VC1=area_cuadraturas(a,b,matlabFunction(1/Acx,"Vars",{x}));
            
A=[abs(MC1),MC2;
   MC1x+VC1,abs(MC2x)];
B=[Mq+g;
   Mqx+Vqx+v]; 
VM=A\B;

V2=-VM(1) + q1*L - (L*L^nq*(q1 - q2))/(L^nq + L^nq*nq);
M2=-VM(2) + -VM(1)*L + (q1*L^2)/2 - (L^nq*L^2*(q1 - q2))/(2*L^nq + L^nq*nq^2 + 3*L^nq*nq);

MEE=[VM(1);-VM(2);V2;-M2]; 
