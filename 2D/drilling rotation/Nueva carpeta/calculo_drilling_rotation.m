E=30000;
nu=0.25;
t=0.2;
xe=[0;1.5;0];
ye=[0;0;1.5];
x21 = xe(2) - xe(1);         y21 = ye(2) - ye(1); 
x32 = xe(3) - xe(2);         y32 = ye(3) - ye(2);
x31 = xe(3) - xe(1);         y31 = ye(3) - ye(1);
x1=xe(1);
x2=xe(2);
x3=xe(3);
y1=ye(1);
y2=ye(2);
y3=ye(3);

plot(xe,ye)

x4=mean([x1, x2]);
y4=mean([y1, y2]);

x5=mean([x3, x2]);
y5=mean([y3, y2]);

x6=mean([x3, x1]);
y6=mean([y3, y1]);

xji = [ x21 x32 x31];   yji = [y21 y32 y31]; 

Lk = hypot(xji, yji);      Sk =-xji./Lk;     Ck = yji./Lk; %% figure 4

L5=Lk(1);
L6=Lk(2);
L7=Lk(3);
C12=Ck(1);
C23=Ck(2);
C31=Ck(3);
S12=Sk(1);
S23=Sk(2);
S31=Sk(3);

a1 =0;
a2 =4/x2;
a3 =-(4*x3)/(x2*y3);
a4 =-4/x2^2;
a5 =-(4*(x2 - 2*x3))/(x2^2*y3);
a6 =(4*x3*(x2 - x3))/(x2^2*y3^2);
b1 =0;
b2 =0;
b3 =0;
b4 =0;
b5 =4/(x2*y3);
b6 =-(4*x3)/(x2*y3^2);
d1 =0;
d2 =0;
d3 =4/y3;
d4 =0;
d5 =-4/(x2*y3);
d6 =-(4*(x2 - x3))/(x2*y3^2);
x=0.75;
y=0;
abd=[a2,a3,a4,a5,a6,b5,b6,d3,d5,d6];
bm =[0,1,0,0,0,0,C12*a2+2*C12*a4*x+C12*a5*y                           ,C23*b5*y                    ,C31*d5*y;
     0,0,0,0,0,1,S12*a3+S12*a5*x+2*S12*a6*y                           ,S23*b5*x+2*S23*b6*y         ,S31*d3+S31*d5*x+2*S31*d6*y;
     0,0,1,0,1,0,C12*a3+S12*a2+C12*a5*x+2*C12*a6*y+2*S12*a4*x+S12*a5*y,C23*b5*x+2*C23*b6*y+S23*b5*y,C31*d3+C31*d5*x+2*C31*d6*y+S31*d5*y];


syms alfa1 alfa2 alfa3 alfa4 alfa5 alfa6 alfa7 alfa8 alfa9 x y 
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

ex=diff(u,x);
ey=diff(v,y);
gamaxy=diff(u,y)+diff(v,x);

ALAFAS = [alfa1 alfa2 alfa3 alfa4 alfa5 alfa6 alfa7 alfa8 alfa9].';

%% ec25
UVGAMA=simplify(equationsToMatrix([ u; v; gama], ALAFAS));
%ec27
B=simplify(equationsToMatrix([ ex; ey; gamaxy], ALAFAS));

A=[subs(UVGAMA,{x,y},{x1,y1});
   subs(UVGAMA,{x,y},{x2,y2});  
   subs(UVGAMA,{x,y},{x3,y3})]; 
De = (E/(1-nu^2)) * [ 1      nu  0
                               nu  1      0
                               0      0      (1-nu)/2 ];
K=double((A^(-1))'*int(int((B'*De*B),y,0,-x+1),x,0,1)*A^(-1))*t;


 