clc
clear
syms Bss Bnn b Db Ds Bsn v Bns dBsk Lk s Bni Bnj Bsi Bsj wi wj phik Ck Sk Bxi Bxj Byi Byj

syms w1 w2 w3 w4 phik Ck Sk Bxi Bxj Byi Byj
syms phi4 phi5 phi6
syms S5 S6 S7 S8 
syms C5 C6 C7 C8 
syms Bx1 Bx2 Bx3 Bx4
syms By1 By2 By3 By4
syms L4 L5 L6  

syms x21 x32 x13  y21 y32 y13
C4 = x21/L4;      S4 = y21/L4;
C5 = x32/L5;      S5 = y32/L5;
C6 = x13/L6;      S6 = y13/L6;

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
dbs4=subs(dBsk_r,{wj,wi,Bxj,Bxi,Byj,Byi,Lk,Ck,Sk,phik},{w2,w1,Bx2,Bx1,By2,By1,L4,C4,S4,phi4});
dbs5=subs(dBsk_r,{wj,wi,Bxj,Bxi,Byj,Byi,Lk,Ck,Sk,phik},{w3,w2,Bx3,Bx2,By3,By2,L5,C5,S5,phi5});
dbs6=subs(dBsk_r,{wj,wi,Bxj,Bxi,Byj,Byi,Lk,Ck,Sk,phik},{w1,w3,Bx1,Bx3,By1,By3,L6,C6,S6,phi6});

Un = [ w1 Bx1 By1 w2 Bx2 By2 w3 Bx3 By3].';
% Equation 41
An = simplify(equationsToMatrix([ dbs4; dbs5; dbs6 ], Un));
% Equation 43
Adb = diag((2/3)*[ L4*(1+phi4) L5*(1+phi5) L6*(1+phi6)]);
% From equation 42, we obtain equation 44
Aw =simplify(Adb*An);


%%
clear
syms gxz1 gyz1 gxz2 gyz2 gxz3 gyz3 gsz4 gsz5 gsz7 gsz6
syms C4 C5 C6 S4 S5 S6 
syms A1 A2 A3
syms L4 L5 L6  
syms x21 x32 x13  y21 y32 y13
syms xi_gl eta_gl x1 x2 x3 x4 y3 y1 y2 y3 y4 landa
N1 = 1-xi_gl-eta_gl;%landa;
N2 = xi_gl;
N3 = eta_gl;

%C4 = x21/L4;      S4 = y21/L4;
%C5 = x32/L5;      S5 = y32/L5;
%C6 = x13/L6;      S6 = y13/L6;

A1=C4*S6-C6*S4;
A2=C5*S4-C4*S5;
A3=C6*S5-C5*S6;

gxzbar = N1*gxz1 + N2*gxz2 + N3*gxz3 ;
gyzbar = N1*gyz1 + N2*gyz2 + N3*gyz3 ;
M6 = equationsToMatrix([ gxzbar; gyzbar ], [ gxz1; gyz1; gxz2; gyz2; gxz3; gyz3]);

gxz1 = ( S6*gsz4 - S4*gsz6)/A1;
gyz1 = (-C6*gsz4 + C4*gsz6)/A1;

gxz2 = ( S4*gsz5 - S5*gsz4)/A2;
gyz2 = (-C4*gsz5 + C5*gsz4)/A2;

gxz3 = ( S5*gsz6 - S6*gsz5)/A3;
gyz3 = (-C5*gsz6 + C6*gsz5)/A3;


M5 = equationsToMatrix([ gxz1; gyz1; gxz2; gyz2; gxz3; gyz3 ], [ gsz4; gsz5; gsz6 ]);

syms phi4 phi5 phi6
M8 = diag(-(2/3)*[ phi4 phi5 phi6]);

Ngamma = M6*M5;
Bsdb = Ngamma*M8
