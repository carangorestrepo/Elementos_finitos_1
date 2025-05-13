function [MEE]=Ax_Ix_As2(b1,b2,nb,h1,h2,nh,q1,q2,nq,bq1,bq2,nbq,a,b,E,G,L,tw,tf,sec,g,v)
%{
MC2=simplify(1/Ix)
MC1=simplify(x/Ix)
Mq=simplify(M/Ix)
MC2x=simplify(1/Ix*x)
MC1x=simplify(x/Ix*x)
Mqx=simplify(M/Ix*x)
Vqx=simplify(V/As2x)
VC1=simplify(1/As2x)
Axial  =simplify(VA/Ax)
C1xfAE =simplify(1/Ax)
%}
if sec==1 
%seccion C
    MC2 =area_cuadraturas(a,b,@(x)(12/(E*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3))));
    MC1 =area_cuadraturas(a,b,@(x)((12*x)/(E*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3))));
    Mq =area_cuadraturas(a,b,@(x)((6*x^2*(2*q2*x^nq - 2*q1*x^nq + 2*L^nq*q1 + 3*L^nq*nq*q1 + L^nq*nq^2*q1))/(E*L^nq*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3)*(nq^2 + 3*nq + 2))));
    MC2x =area_cuadraturas(a,b,@(x)((12*x)/(E*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3))));
    MC1x =area_cuadraturas(a,b,@(x)((12*x^2)/(E*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3))));
    Mqx =area_cuadraturas(a,b,@(x)((6*x^3*(2*q2*x^nq - 2*q1*x^nq + 2*L^nq*q1 + 3*L^nq*nq*q1 + L^nq*nq^2*q1))/(E*L^nq*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3)*(nq^2 + 3*nq + 2))));
    Vqx =area_cuadraturas(a,b,@(x)((q1*x - (x*x^nq*(q1 - q2))/(L^nq + L^nq*nq))/(G*tw*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh))));
    VC1 =area_cuadraturas(a,b,@(x)(1/(G*tw*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh))));
    Axial =area_cuadraturas(a,b,@(x)((bq1*x - (x*x^nbq*(bq1 - bq2))/(L^nbq + L^nbq*nbq))/(E*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - (tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))));
    C1xfAE =area_cuadraturas(a,b,@(x)(1/(E*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - (tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))));
elseif sec==2  
%Seccion rectangular 
    MC2 =area_cuadraturas(a,b,@(x)(12/(E*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3)));
    MC1 =area_cuadraturas(a,b,@(x)((12*x)/(E*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3)));
    Mq =area_cuadraturas(a,b,@(x)((6*x^2*(2*q2*x^nq - 2*q1*x^nq + 2*L^nq*q1 + 3*L^nq*nq*q1 + L^nq*nq^2*q1))/(E*L^nq*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3*(nq^2 + 3*nq + 2))));
    MC2x =area_cuadraturas(a,b,@(x)((12*x)/(E*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3)));
    MC1x =area_cuadraturas(a,b,@(x)((12*x^2)/(E*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3)));
    Mqx =area_cuadraturas(a,b,@(x)((6*x^3*(2*q2*x^nq - 2*q1*x^nq + 2*L^nq*q1 + 3*L^nq*nq*q1 + L^nq*nq^2*q1))/(E*L^nq*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3*(nq^2 + 3*nq + 2)))); 
    Vqx =area_cuadraturas(a,b,@(x)((6*(q1*x - (x*x^nq*(q1 - q2))/(L^nq + L^nq*nq)))/(5*G*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh))));
    VC1 =area_cuadraturas(a,b,@(x)(6/(5*G*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh))));
    Axial =area_cuadraturas(a,b,@(x)((bq1*x - (x*x^nbq*(bq1 - bq2))/(L^nbq + L^nbq*nbq))/(E*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh))));
    C1xfAE =area_cuadraturas(a,b,@(x)(1/(E*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh))));
elseif sec==3
%Circular solida
    MC2 =area_cuadraturas(a,b,@(x)(64/(E*pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4)));
    MC1 =area_cuadraturas(a,b,@(x)((64*x)/(E*pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4)));
    Mq =area_cuadraturas(a,b,@(x)((32*x^2*(2*q2*x^nq - 2*q1*x^nq + 2*L^nq*q1 + 3*L^nq*nq*q1 + L^nq*nq^2*q1))/(E*L^nq*pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4*(nq^2 + 3*nq + 2))));
    MC2x =area_cuadraturas(a,b,@(x)((64*x)/(E*pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4)));
    MC1x =area_cuadraturas(a,b,@(x)((64*x^2)/(E*pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4)));
    Mqx =area_cuadraturas(a,b,@(x)((32*x^3*(2*q2*x^nq - 2*q1*x^nq + 2*L^nq*q1 + 3*L^nq*nq*q1 + L^nq*nq^2*q1))/(E*L^nq*pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4*(nq^2 + 3*nq + 2))));
    Vqx =area_cuadraturas(a,b,@(x)((40*(q1*x - (x*x^nq*(q1 - q2))/(L^nq + L^nq*nq)))/(9*G*pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)));
    VC1 =area_cuadraturas(a,b,@(x)(40/(9*G*pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)));
    Axial =area_cuadraturas(a,b,@(x)((4*(bq1*x - (x*x^nbq*(bq1 - bq2))/(L^nbq + L^nbq*nbq)))/(E*pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)));
    C1xfAE =area_cuadraturas(a,b,@(x)(4/(E*pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)));
elseif sec==4
%%Seccion I
    MC2 =area_cuadraturas(a,b,@(x)(12/(E*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3))));
    MC1 =area_cuadraturas(a,b,@(x)((12*x)/(E*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3))));
    Mq =area_cuadraturas(a,b,@(x)((6*x^2*(2*q2*x^nq - 2*q1*x^nq + 2*L^nq*q1 + 3*L^nq*nq*q1 + L^nq*nq^2*q1))/(E*L^nq*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3)*(nq^2 + 3*nq + 2))));
    MC2x =area_cuadraturas(a,b,@(x)((12*x)/(E*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3))));
    MC1x =area_cuadraturas(a,b,@(x)((12*x^2)/(E*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3))));
    Mqx =area_cuadraturas(a,b,@(x)((6*x^3*(2*q2*x^nq - 2*q1*x^nq + 2*L^nq*q1 + 3*L^nq*nq*q1 + L^nq*nq^2*q1))/(E*L^nq*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3)*(nq^2 + 3*nq + 2))));
    Vqx =area_cuadraturas(a,b,@(x)((q1*x - (x*x^nq*(q1 - q2))/(L^nq + L^nq*nq))/(G*tw*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh))));
    VC1 =area_cuadraturas(a,b,@(x)(1/(G*tw*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh))));
    Axial =area_cuadraturas(a,b,@(x)((bq1*x - (x*x^nbq*(bq1 - bq2))/(L^nbq + L^nbq*nbq))/(E*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - (tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))));
    C1xfAE =area_cuadraturas(a,b,@(x)(1/(E*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - (tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))));
elseif sec==5
%Seccion circular hueca
    MC2 =area_cuadraturas(a,b,@(x)(-64/(E*(pi*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4 - pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4)))); 
    MC1 =area_cuadraturas(a,b,@(x)(-(64*x)/(E*(pi*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4 - pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4)))); 
    Mq =area_cuadraturas(a,b,@(x)(-(32*x^2*(2*q2*x^nq - 2*q1*x^nq + 2*L^nq*q1 + 3*L^nq*nq*q1 + L^nq*nq^2*q1))/(E*L^nq*(pi*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4 - pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4)*(nq^2 + 3*nq + 2)))); 
    MC2x =area_cuadraturas(a,b,@(x)(-(64*x)/(E*(pi*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4 - pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4))));
    MC1x =area_cuadraturas(a,b,@(x)(-(64*x^2)/(E*(pi*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4 - pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4))));
    Mqx =area_cuadraturas(a,b,@(x)(-(32*x^3*(2*q2*x^nq - 2*q1*x^nq + 2*L^nq*q1 + 3*L^nq*nq*q1 + L^nq*nq^2*q1))/(E*L^nq*(pi*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4 - pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4)*(nq^2 + 3*nq + 2))));
    Vqx =area_cuadraturas(a,b,@(x)(-(x*(10*h1*(b - a)^nh - 10*h1*(x - a)^nh + 10*h2*(x - a)^nh)*(q2*x^nq - q1*x^nq + L^nq*q1 + L^nq*nq*q1))/(G*L^nq*tf*pi*(nq + 1)*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)*(5*h1*(b - a)^nh + 8*tf*(b - a)^nh - 5*h1*(x - a)^nh + 5*h2*(x - a)^nh))));
    VC1 =area_cuadraturas(a,b,@(x)(-(10*h1*(b - a)^nh - 10*h1*(x - a)^nh + 10*h2*(x - a)^nh)/(G*tf*pi*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)*(5*h1*(b - a)^nh + 8*tf*(b - a)^nh - 5*h1*(x - a)^nh + 5*h2*(x - a)^nh))));
    Axial =area_cuadraturas(a,b,@(x)(-(bq1*x - (x*x^nbq*(bq1 - bq2))/(L^nbq + L^nbq*nbq))/(E*((pi*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)/4 - (pi*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)/4))));
    C1xfAE =area_cuadraturas(a,b,@(x)((b - a)^nh/(E*tf*pi*(h1*(b - a)^nh - tf*(b - a)^nh - h1*(x - a)^nh + h2*(x - a)^nh))));
elseif sec==6
%Tubular hueca 
    MC2 =area_cuadraturas(a,b,@(x)(12/(E*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (2*tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3)))); 
    MC1 =area_cuadraturas(a,b,@(x)((12*x)/(E*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (2*tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3))));
    Mq =area_cuadraturas(a,b,@(x)((6*x^2*(2*q2*x^nq - 2*q1*x^nq + 2*L^nq*q1 + 3*L^nq*nq*q1 + L^nq*nq^2*q1))/(E*L^nq*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (2*tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3)*(nq^2 + 3*nq + 2)))); 
    MC2x =area_cuadraturas(a,b,@(x)((12*x)/(E*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (2*tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3)))); 
    MC1x =area_cuadraturas(a,b,@(x)((12*x^2)/(E*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (2*tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3))));
    Mqx =area_cuadraturas(a,b,@(x)((6*x^3*(2*q2*x^nq - 2*q1*x^nq + 2*L^nq*q1 + 3*L^nq*nq*q1 + L^nq*nq^2*q1))/(E*L^nq*((b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - (2*tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3)*(nq^2 + 3*nq + 2))));
    Vqx =area_cuadraturas(a,b,@(x)((q1*x - (x*x^nq*(q1 - q2))/(L^nq + L^nq*nq))/(2*G*tw*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh))));
    VC1 =area_cuadraturas(a,b,@(x)(1/(2*G*tw*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh))));
    Axial =area_cuadraturas(a,b,@(x)(-(bq1*x - (x*x^nbq*(bq1 - bq2))/(L^nbq + L^nbq*nbq))/(E*((2*tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - (b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))));
    C1xfAE =area_cuadraturas(a,b,@(x)(-1/(E*((2*tw - b1 + ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - (b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))));
elseif sec==7
%Perfil T
    MC2 =area_cuadraturas(a,b,@(x)(-12/(E*(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - tf^3*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb) + 3*tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)*(h1 - tf + (2*((tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)/2 - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(tf/2 - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))/(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)) - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2 - 3*tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*h1 - tf + (2*((tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)/2 - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(tf/2 - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))/(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)) - (2*(h1 - h2)*(x - a)^nh)/(b - a)^nh)^2))));
    MC1 =area_cuadraturas(a,b,@(x)(-(12*x)/(E*(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - tf^3*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb) + 3*tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)*(h1 - tf + (2*((tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)/2 - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(tf/2 - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))/(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)) - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2 - 3*tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*h1 - tf + (2*((tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)/2 - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(tf/2 - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))/(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)) - (2*(h1 - h2)*(x - a)^nh)/(b - a)^nh)^2))));
    Mq =area_cuadraturas(a,b,@(x)(-(6*x^2*(2*q2*x^nq - 2*q1*x^nq + 2*L^nq*q1 + 3*L^nq*nq*q1 + L^nq*nq^2*q1))/(E*L^nq*(nq^2 + 3*nq + 2)*(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - tf^3*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb) + 3*tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)*(h1 - tf + (2*((tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)/2 - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(tf/2 - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))/(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)) - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2 - 3*tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*h1 - tf + (2*((tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)/2 - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(tf/2 - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))/(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)) - (2*(h1 - h2)*(x - a)^nh)/(b - a)^nh)^2))));
    MC2x =area_cuadraturas(a,b,@(x)(-(12*x)/(E*(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - tf^3*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb) + 3*tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)*(h1 - tf + (2*((tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)/2 - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(tf/2 - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))/(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)) - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2 - 3*tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*h1 - tf + (2*((tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)/2 - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(tf/2 - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))/(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)) - (2*(h1 - h2)*(x - a)^nh)/(b - a)^nh)^2))));
    MC1x =area_cuadraturas(a,b,@(x)(-(12*x^2)/(E*(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - tf^3*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb) + 3*tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)*(h1 - tf + (2*((tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)/2 - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(tf/2 - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))/(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)) - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2 - 3*tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*h1 - tf + (2*((tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)/2 - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(tf/2 - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))/(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)) - (2*(h1 - h2)*(x - a)^nh)/(b - a)^nh)^2))));
    Mqx =area_cuadraturas(a,b,@(x)(-(6*x^3*(2*q2*x^nq - 2*q1*x^nq + 2*L^nq*q1 + 3*L^nq*nq*q1 + L^nq*nq^2*q1))/(E*L^nq*(nq^2 + 3*nq + 2)*(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3 - tf^3*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb) + 3*tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)*(h1 - tf + (2*((tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)/2 - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(tf/2 - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))/(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)) - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2 - 3*tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(2*h1 - tf + (2*((tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2)/2 - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)*(tf/2 - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh)))/(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)) - (2*(h1 - h2)*(x - a)^nh)/(b - a)^nh)^2))));
    Vqx =area_cuadraturas(a,b,@(x)((q1*x - (x*x^nq*(q1 - q2))/(L^nq + L^nq*nq))/(2*G*tw*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh))));
    VC1 =area_cuadraturas(a,b,@(x)(1/(2*G*tw*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh))));
    Axial =area_cuadraturas(a,b,@(x)(-(bq1*x - (x*x^nbq*(bq1 - bq2))/(L^nbq + L^nbq*nbq))/(E*(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)))));
    C1xfAE =area_cuadraturas(a,b,@(x)(-1/(E*(tw*(tf - h1 + ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - tf*(b1 - ((b1 - b2)*(x - a)^nb)/(b - a)^nb)))));
elseif sec==8
% perfil L 
    MC2 =area_cuadraturas(a,b,@(x)(-(12*tf - 24*h1 + (24*(h1 - h2)*(x - a)^nh)/(b - a)^nh)/(E*tf*(2*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4 - tf^3*(6*h1 - (6*(h1 - h2)*(x - a)^nh)/(b - a)^nh) + 8*tf^2*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2 + tf^4 - 4*tf*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3))));
    MC1 =area_cuadraturas(a,b,@(x)(-(x*(12*tf - 24*h1 + (24*(h1 - h2)*(x - a)^nh)/(b - a)^nh))/(E*tf*(2*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4 - tf^3*(6*h1 - (6*(h1 - h2)*(x - a)^nh)/(b - a)^nh) + 8*tf^2*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2 + tf^4 - 4*tf*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3))));
    Mq =area_cuadraturas(a,b,@(x)(-(x^2*(12*tf - 24*h1 + (24*(h1 - h2)*(x - a)^nh)/(b - a)^nh)*(2*q2*x^nq - 2*q1*x^nq + 2*L^nq*q1 + 3*L^nq*nq*q1 + L^nq*nq^2*q1))/(2*E*L^nq*tf*(nq^2 + 3*nq + 2)*(2*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4 - tf^3*(6*h1 - (6*(h1 - h2)*(x - a)^nh)/(b - a)^nh) + 8*tf^2*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2 + tf^4 - 4*tf*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3))));
    MC2x =area_cuadraturas(a,b,@(x)(-(x*(12*tf - 24*h1 + (24*(h1 - h2)*(x - a)^nh)/(b - a)^nh))/(E*tf*(2*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4 - tf^3*(6*h1 - (6*(h1 - h2)*(x - a)^nh)/(b - a)^nh) + 8*tf^2*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2 + tf^4 - 4*tf*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3))));
    MC1x =area_cuadraturas(a,b,@(x)(-(x^2*(12*tf - 24*h1 + (24*(h1 - h2)*(x - a)^nh)/(b - a)^nh))/(E*tf*(2*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4 - tf^3*(6*h1 - (6*(h1 - h2)*(x - a)^nh)/(b - a)^nh) + 8*tf^2*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2 + tf^4 - 4*tf*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3))));
    Mqx =area_cuadraturas(a,b,@(x)(-(x^3*(12*tf - 24*h1 + (24*(h1 - h2)*(x - a)^nh)/(b - a)^nh)*(2*q2*x^nq - 2*q1*x^nq + 2*L^nq*q1 + 3*L^nq*nq*q1 + L^nq*nq^2*q1))/(2*E*L^nq*tf*(nq^2 + 3*nq + 2)*(2*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^4 - tf^3*(6*h1 - (6*(h1 - h2)*(x - a)^nh)/(b - a)^nh) + 8*tf^2*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^2 + tf^4 - 4*tf*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh)^3))));
    Vqx =area_cuadraturas(a,b,@(x)((q1*x - (x*x^nq*(q1 - q2))/(L^nq + L^nq*nq))/(G*tf*((5*h1)/6 - (5*(h1 - h2)*(x - a)^nh)/(6*(b - a)^nh)))));
    VC1 =area_cuadraturas(a,b,@(x)(6/(5*G*tf*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh))));
    Axial =area_cuadraturas(a,b,@(x)((bq1*x - (x*x^nbq*(bq1 - bq2))/(L^nbq + L^nbq*nbq))/(E*(2*tf*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - tf^2))));
    C1xfAE =area_cuadraturas(a,b,@(x)(1/(E*(2*tf*(h1 - ((h1 - h2)*(x - a)^nh)/(b - a)^nh) - tf^2))));
end

A=[abs(MC1),MC2;
   MC1x+VC1,abs(MC2x)];
B=[Mq+g;
   Mqx+Vqx+v]; 
VM=A\B;
V2=-VM(1) + q1*L - (L*L^nq*(q1 - q2))/(L^nq + L^nq*nq);
M2=-VM(2) + -VM(1)*L + (q1*L^2)/2 - (L^nq*L^2*(q1 - q2))/(2*L^nq + L^nq*nq^2 + 3*L^nq*nq);

A=C1xfAE;
B=Axial+v;
RA=A\B;
RB=-RA + bq1*L - (L*L^nbq*(bq1 - bq2))/(L^nbq + L^nbq*nbq); 
MEE=[RA;VM(1);-VM(2);RB;V2;-M2]; 
 