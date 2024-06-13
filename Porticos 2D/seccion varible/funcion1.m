function [R]=funcion1(tf,tw,b1,b2,nb,h1,h2,nh,q1,q2,nq,L,E,G,v,xi,xf,sec)
%syms xf xi x b2 b1 h2 h1 q1 q2 C1 C2 nb nh nq 
%bx=((b2-b1)/(xf-xi)^nb*(x-xi)^nb+b1);
%hx=((h2-h1)/(xf-xi)^nh*(x-xi)^nh+h1);
%qx=((q2-q1)/(xf-xi)^nq*(x-xi)^nq+q1);

%Ax=bx*hx
%Ax=((h1 - ((h1 - h2)*(x - xi)^nh)/(xf - xi)^nh)*(b1 - ((b1 - b2)*(x - xi)^nb)/(xf - xi)^nb))
%A=qx+C1
%def=expand(A*x/Ax+C2)
syms x
V= q1*x - (x*x^nq*(q1 - q2))/(L^nq + L^nq*nq);

%[Ax,Ix,As2x]=secciones_1(sec,b1,b2,nb,h1,h2,nh,G,E,tf,tw,x1,x2)
[Ae,~,~]=    secciones_1(sec,b1,b2,nb,h1,h2,nh,G,E,tf,tw,0,L);
Axial  =area_cuadraturas(0,L,matlabFunction((V/Ae),"Vars",{x}));
C1xfAE =area_cuadraturas(0,L,matlabFunction((1/Ae),"Vars",{x}));
A=C1xfAE;
B=Axial+v;
RA=A\B;
RB=-RA + q1*L - (L*L^nq*(q1 - q2))/(L^nq + L^nq*nq); 
R=[RA;RB];