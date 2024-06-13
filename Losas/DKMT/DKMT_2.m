clc
clear
syms Bss Bnn b Db Ds Bsn v Bns dBsk Lk s Bni Bnj Bsi Bsj wi wj phik Ck Sk Bxi Bxj Byi Byj

syms w1 w2 w3 w4 phik Ck Sk Bxi Bxj Byi Byj
syms phi5 phi6 phi7 phi8 
syms S5 S6 S7 S8 
syms C5 C6 C7 C8 
syms Bx1 Bx2 Bx3 Bx4
syms By1 By2 By3 By4
syms L5 L6 L7 L8 

syms x43 x32 x21 x14 y43 y32 y21 y14
C5 = x21/L5;      S5 = y21/L5;
C6 = x32/L6;      S6 = y32/L6;
C7 = x43/L7;      S7 = y43/L7;
C8 = x14/L8;      S8 = y14/L8;

Mss=Db*(Bss+v*Bnn);

Mns=Db*(1-v)/2*(Bsn+Bns);

Ts=Mss+Mns;
gamaaz=Ts/Ds;  %16
Bsi=Ck*Bxi+Sk*Byi;
Bsj=Ck*Bxj+Sk*Byj;
Bn=(1-s/Lk)*Bni+(s/Lk)*Bnj;
Bs=(1-s/Lk)*Bsi+(s/Lk)*Bsj+4*s/Lk*(1-s/Lk)*dBsk;
ws=(1-s/Lk)*wi +(s/Lk)*wj;
gamazk=-2/3*phik*dBsk;%%22a
gamask=Bs+diff(ws,s);%%%32
ec=int(gamask-gamazk,s,0,Lk)==0;%% 35b
dBsk_r=solve(ec,dBsk);%% 35b despejo deltaBsk

% eq 40b
% dbsk = (wj - wi + Lk*(ck*bxi + sk*byi)/2 + Lk*(ck*bxj + sk*byj)/2)/(-2*Lk*(1+phik)/3);
dbs5=subs(dBsk_r,{wj,wi,Bxj,Bxi,Byj,Byi,Lk,Ck,Sk,phik},{w2,w1,Bx2,Bx1,By2,By1,L5,C5,S5,phi5});
dbs6=subs(dBsk_r,{wj,wi,Bxj,Bxi,Byj,Byi,Lk,Ck,Sk,phik},{w3,w2,Bx3,Bx2,By3,By2,L6,C6,S6,phi6});
dbs7=subs(dBsk_r,{wj,wi,Bxj,Bxi,Byj,Byi,Lk,Ck,Sk,phik},{w4,w3,Bx4,Bx3,By4,By3,L7,C7,S7,phi7});
dbs8=subs(dBsk_r,{wj,wi,Bxj,Bxi,Byj,Byi,Lk,Ck,Sk,phik},{w1,w4,Bx1,Bx4,By1,By4,L8,C8,S8,phi8});

Un = [ w1 Bx1 By1 w2 Bx2 By2 w3 Bx3 By3 w4 Bx4 By4 ].';
% Equation 41
An = simplify(equationsToMatrix([ dbs5; dbs6; dbs7; dbs8 ], Un));
% Equation 43
Adb = diag((2/3)*[ L5*(1+phi5) L6*(1+phi6) L7*(1+phi7) L8*(1+phi8) ]);
% From equation 42, we obtain equation 44
Aw =simplify(Adb*An);




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

J = [ diff(x,xi)   diff(y,xi)
      diff(x,eta)  diff(y,eta) ];