clc
clear
syms x1 x2 x3 y1 y2 y3 eta_gl xi_gl phi4 phi5 phi6
%% An evaluation on the performance of two simple triangular
%% bending plate elements
nu = 0.3;         %         coeficiente de Poisson
h  = 0.05;        % [m]     espesor de la losa
xi_gl=0.166666666666670;
eta_gl=0.166666666666670;

xe=[0.142857142857143;0;0];
ye=[0;0.142857142857143;0];

 
x21=xe(2)-xe(1);
x32=xe(3)-xe(2);
x13=xe(1)-xe(3);

y21=ye(2)-ye(1);
y32=ye(3)-ye(2);
y13=ye(1)-ye(3);

L12=(x21^2+y21^2)^(1/2);
L23=(x32^2+y32^2)^(1/2);
L31=(x13^2+y13^2)^(1/2);

C12=x21/L12;
C23=x32/L23;
C31=x13/L31;
S12=y21/L12;
S23=y32/L23;
S31=y13/L31;

L4=hypot(x21,y21);
L5=hypot(x32,y32);
L6=hypot(x13,y13);
phi4=(2/((5/6)*(1 - nu))) .* (h./L4).^2;
phi5=(2/((5/6)*(1 - nu))) .* (h./L5).^2;
phi6=(2/((5/6)*(1 - nu))) .* (h./L6).^2; 

A1=C12*S31-C31*S12;
A2=C23*S12-C12*S23;
A3=C31*S23-C23*S31;
landa=1-xi_gl-eta_gl;

Bsgma=[(S31*landa/A1-S23*xi_gl/A2),(S12*xi_gl/A2-S31*eta_gl/A3),(S23*xi_gl/A3-S12*landa/A1);
      (C23*xi_gl/A2-C31*landa/A1),(C31*eta_gl/A3-C12*xi_gl/A2),(C12*landa/A1-C23*eta_gl/A3)];
  
Au=1/2*[-2/L12,C12,S12, 2/L12,C12,S12,     0,  0,  0;
             0,  0,  0,-2/L23,C23,S23, 2/L23,C23,S23;
         2/L31,C31,S31,     0,  0,  0,-2/L31,C31,S31];
 
Adelta=-2/3*[1+phi4,0     ,     0;
                  0,1+phi5,     0;
                  0,     0,1+phi6];
              
Aphi=-2/3*[phi4,   0,   0;
              0,phi5,   0;
              0,   0,phi6];
          
Bs=(Bsgma*Aphi*Adelta^(-1)*Au);
