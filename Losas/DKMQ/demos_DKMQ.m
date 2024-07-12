clc
clear
syms Bss Bnn b Db Ds Bsn v Bns dBsk Lk s Bni Bnj Bsi Bsj wi wj phik Ck Sk Bxi Bxj Byi Byj

syms w1 w2 w3 w4 phik Ck Sk Bxi Bxj Byi Byj
syms phi4 phi5 phi6 phi7 phi8
syms S5 S6 S7 S8 
syms C5 C6 C7 C8 
syms Bx1 Bx2 Bx3 Bx4
syms By1 By2 By3 By4
syms L4 L5 L6 L7 L8  

syms x21 x32 x13 x43  x14 y21  y32 y13 y43 y14
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

Un = [ w1 Bx1 By1 w2 Bx2 By2 w3 Bx3 By3 w4 Bx4 By4].';
% Equation 41
An = simplify(equationsToMatrix([ dbs5; dbs6; dbs7 ;dbs8  ], Un));
% Equation 43
Adb = diag((2/3)*[ L5*(1+phi5) L6*(1+phi6) L7*(1+phi7) L8*(1+phi8)]);
% From equation 42, we obtain equation 44
Aw =simplify(Adb*An);

%%
clear
syms gxz1 gyz1 gxz2 gyz2 gxz3 gyz3 gsz4 gsz5 gsz7 gsz6 gxz4 gyz4 gsz8
syms C4 C5 C6 C7 C8
syms S4 S5 S6 S7 S8
syms A1 A2 A3 A4
syms L4 L5 L6 L7 L8 
syms x21 x32 x13  y21 y32 y13
syms xi_gl eta_gl x1 x2 x3 x4 y3 y1 y2 y3 y4 landa
N1 = 1/4*(1-xi_gl)*(1-eta_gl);
N2 = 1/4*(1+xi_gl)*(1-eta_gl);
N3 = 1/4*(1+xi_gl)*(1+eta_gl);
N4 = 1/4*(1-xi_gl)*(1+eta_gl);

A1=C5*S8-C8*S5;
A2=C6*S5-C5*S6;
A3=C7*S6-C6*S7;
A4=C8*S7-C7*S8;

gxzbar = N1*gxz1 + N2*gxz2 + N3*gxz3 + N4*gxz4 ;
gyzbar = N1*gyz1 + N2*gyz2 + N3*gyz3 + N4*gyz4;
M6 = equationsToMatrix([ gxzbar; gyzbar ], [ gxz1; gyz1; gxz2; gyz2; gxz3; gyz3; gxz4; gyz4]);

gxz1 = ( S8*gsz5 - S5*gsz8)/A1;
gyz1 = (-C8*gsz5 + C5*gsz8)/A1;

gxz2 = ( S5*gsz6 - S6*gsz5)/A2;
gyz2 = (-C5*gsz6 + C6*gsz5)/A2;

gxz3 = ( S6*gsz7 - S7*gsz6)/A3;
gyz3 = (-C6*gsz7 + C7*gsz6)/A3;

gxz4 = ( S7*gsz8 - S8*gsz7)/A4;
gyz4 = (-C7*gsz8 + C8*gsz7)/A4;

M5 = equationsToMatrix([ gxz1; gyz1; gxz2; gyz2; gxz3; gyz3; gxz4; gyz4 ], [ gsz5; gsz6; gsz7; gsz8 ]);

syms phi5 phi6 phi6 phi7 phi8
M8 = diag(-(2/3)*[ phi5 phi6 phi7 phi8]);

Ngamma = M6*M5;
Bsdb = Ngamma*M8;
