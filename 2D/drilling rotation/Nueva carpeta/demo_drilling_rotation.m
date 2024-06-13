clc
clear
syms a1 a2 a3 a4 a5 a6 x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6 
syms b1 b2 b3 b4 b5 b6
syms d1 d2 d3 d4 d5 d6
x1=0;
y1=0;
y2=0;
%x2=1;
%y3=1;
%x3=0;
x4=x2/2;
y4=y2/2;
x5=(x2+x3)/2;
y5=(y2+y3)/2;
x6=x3/2;
y6=y3/2;

%% coordenadas elemetntp
xy=[x1,y1;
    x4,y4;
    x2,y2;
    x5,y5;
    x3,y3;
    x6,y6;
    x1,y1];
%plot(xy(:,1),xy(:,2)),hold on
%plot(xy(:,1),xy(:,2),'*r'),hold on
%{
L12=((x1-x2)^2+(y1-y2)^2)^(1/2);
L23=((x2-x3)^2+(y2-y3)^2)^(1/2);
L31=((x3-x3)^2+(y2-y3)^2)^(1/2);
C12=(x1-x2)/L12;
C23=(x2-x3)/L23;
C31=(x3-x1)/L31;
S12=(y1-y2)/L12;
S23=(y2-y3)/L23;
S31=(y3-y1)/L31;
%}
%% ec 14 haciendo 1 el lado 4
ec1a=a1+a2*x1+a3*y1+a4*x1^2+a5*x1*y1+a6*y1^2==0;
ec4a=a1+a2*x4+a3*y4+a4*x4^2+a5*x4*y4+a6*y4^2==1;
ec2a=a1+a2*x2+a3*y2+a4*x2^2+a5*x2*y2+a6*y2^2==0;
ec5a=a1+a2*x5+a3*y5+a4*x5^2+a5*x5*y5+a6*y5^2==0;
ec3a=a1+a2*x3+a3*y3+a4*x3^2+a5*x3*y3+a6*y3^2==0;
ec6a=a1+a2*x6+a3*y6+a4*x6^2+a5*x6*y6+a6*y6^2==0;

[a1,a2,a3,a4,a5,a6]=solve(ec1a,ec2a,ec3a,ec4a,ec5a,ec6a,[a1,a2,a3,a4,a5,a6]);

%% ec 18 haciendo 1 el lado 5
ec1b=b1+b2*x1+b3*y1+b4*x1^2+b5*x1*y1+b6*y1^2==0;
ec4b=b1+b2*x4+b3*y4+b4*x4^2+b5*x4*y4+b6*y4^2==0;
ec2b=b1+b2*x2+b3*y2+b4*x2^2+b5*x2*y2+b6*y2^2==0;
ec5b=b1+b2*x5+b3*y5+b4*x5^2+b5*x5*y5+b6*y5^2==1;
ec3b=b1+b2*x3+b3*y3+b4*x3^2+b5*x3*y3+b6*y3^2==0;
ec6b=b1+b2*x6+b3*y6+b4*x6^2+b5*x6*y6+b6*y6^2==0;

[b1,b2,b3,b4,b5,b6]=solve(ec1b,ec2b,ec3b,ec4b,ec5b,ec6b,[b1,b2,b3,b4,b5,b6]);

%% ec 19 haciendo 1 el lado 6

ec1d=d1+d2*x1+d3*y1+d4*x1^2+d5*x1*y1+d6*y1^2==0;
ec4d=d1+d2*x4+d3*y4+d4*x4^2+d5*x4*y4+d6*y4^2==0;
ec2d=d1+d2*x2+d3*y2+d4*x2^2+d5*x2*y2+d6*y2^2==0;
ec5d=d1+d2*x5+d3*y5+d4*x5^2+d5*x5*y5+d6*y5^2==0;
ec3d=d1+d2*x3+d3*y3+d4*x3^2+d5*x3*y3+d6*y3^2==0;
ec6d=d1+d2*x6+d3*y6+d4*x6^2+d5*x6*y6+d6*y6^2==1;

[d1,d2,d3,d4,d5,d6]=solve(ec1d,ec2d,ec3d,ec4d,ec5d,ec6d,[d1,d2,d3,d4,d5,d6]);

syms C12 C23 C31 S12 S23 S31 P4 P5 P6
syms alfa1 alfa2 alfa3 alfa4 alfa5 alfa6 alfa7 alfa8 alfa9 x y
%syms a1 a2 a3 a4 a5 a6 x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6
%syms b1 b2 b3 b4 b5 b6
%syms d1 d2 d3 d4 d5 d6
syms a2 a3 a4 a5 a6  b5  b6 d3 d5 d6
%% ec13
dp4=(a1+a2*x+a3*y+a4*x^2+a5*x*y+a6*y^2);
dp5=(b1+b2*x+b3*y+b4*x^2+b5*x*y+b6*y^2);
dp6=(d1+d2*x+d3*y+d4*x^2+d5*x*y+d6*y^2);
%%ec20
up4=dp4*C12;
up5=dp5*C23;
up6=dp6*C31;

vp4=dp4*S12;
vp5=dp5*S23;
vp6=dp6*S31;

P4=alfa7;
P5=alfa8;
P6=alfa9;
%ec20-ec21
u=alfa1+alfa2*x+alfa3*y+up4*P4+up5*P5+up6*P6;
v=alfa4+alfa5*x+alfa6*y+vp4*P4+vp5*P5+vp6*P6;
%% ec7
gama=1/2*(diff(v,x)-diff(u,y));
%% la ecuaacion 26 fue revisada y ccomprobada con documento
ex=diff(u,x);
ey=diff(v,y);
gamaxy=diff(u,y)+diff(v,x);

ALAFAS = [alfa1 alfa2 alfa3 alfa4 alfa5 alfa6 alfa7 alfa8 alfa9].';

%% ec25
UVGAMA=simplify(equationsToMatrix([ u; v; gama], ALAFAS));
%ec27
B=simplify(equationsToMatrix([ ex; ey; gamaxy], ALAFAS));

%L12=sqrt((x2-x1)^2+(y2-y1)^2);
%L23=sqrt((x3-x1)^2+(y3-y2)^2);
%L31=sqrt((x1-x3)^2+(y1-y3)^2);
%C12=(x2-x1)/L12;
%C23=(x3-x1)/L23;
%C31=(x0-x3)/L31;

%S12=(y2-y1)/L12;
%S23=(y3-y1)/L23;
%S31=(y1-y3)/L31;

A=[subs(UVGAMA,{x,y,C12,S12},{x1,x2,0,-1});
    subs(UVGAMA,{x,y,x1,y2,C12,-S12},{x2,y2,0,0,0,-1});  
    subs(UVGAMA,{x,y,x1,y2,C12,-S12},{x3,y3,0,0,0,-1})]; 
%K=(A^-1)'*B*D*B*(A^-1);

syms xi yi
%Notas la ecuacion 30 esta acorde al documento 
Ai=[[1 xi yi 0 0 0 0 0 0];
    [0 0 0 1 xi yi 0 0 0];
    [0,0,-1/2,0,1/2,0,1/2*S12*a2-1/2*C12*a3-1/2*C12*a5*xi-C12*a6*yi+S12*a4*xi+1/2*S12*a5*yi,1/2*S23*b5*yi-C23*b6*yi-1/2*C23*b5*xi,1/2*S31*d5*yi-1/2*C31*d5*xi-C31*d6*yi-1/2*C31*d3]];

As=[subs(Ai,{xi,yi},{x1,y1});
    subs(Ai,{xi,yi},{x2,y2});
    subs(Ai,{xi,yi},{x3,y3})];
%Notas la ecuacion 27 dle documento tiene errores
bs= [0,1,0,0,0,0,C12*a2 + 2*C12*a4*x + C12*a5*y                                 ,y*b5*C23*b5                   ,C31*d5*y;
     0,0,0,0,0,1,S12*a3 + S12*a5*x + 2*S12*a6*y                                 ,(x*b5+ 2*y*b6)*S23            ,S31*d3 + S31*d5*x + 2*S31*d6*y;
     0,0,1 0,1,0,C12*a3 + S12*a2 + C12*a5*x + 2*C12*a6*y + 2*S12*a4*x + S12*a5*y,x*C23*b5+2*y*C23*b6+y*S23*b5  ,C31*d3 + C31*d5*x + 2*C31*d6*y + S31*d5*y];
 

 
 
 