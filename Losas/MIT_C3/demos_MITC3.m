%% Demonstration of the equations relative to the MITC3 finite element in:
%

%% On the formulation and evaluation of old and new efficient low order
%% triangular plate bending elements with shear effects
%% https://doi.org/10.1007/s00466-021-02020-6
clear, clc
disp('EQUATIONS Bs KATILI ET. AL.')
disp('*** Calculation of the inverse Jacobian j = inv(J) ***')

% Jacobian and inverse Jacobian
% isoparametric interpolation

syms xi eta x1 x2 x3 x4 y1 y2 y3 y4
N1 = 1-xi-eta; % = ec 6
N2 = xi;
N3 = eta;


x = N1*x1 + N2*x2 + N3*x3 ;
y = N1*y1 + N2*y2 + N3*y3 ;

disp('J = '); 
J = [ diff(x,xi)   diff(y,xi)
      diff(x,eta)  diff(y,eta) ]

disp('det(J) = ');
detJ = det(J);

invJ = simplify(inv(J));
% invJ = [ diff(xi,x)  diff(eta,x)
%          diff(xi,y)  diff(eta,y) ];

disp('inv(J)*det(J) = '); collect(invJ*detJ,{'r','s','xi','eta'})

%% Equation , equation 10
disp('*** DEMO, EQUATION 10 ***')

syms gs4 gs5 gs6 
g_xieta = [ gs4 gs5 gs6].';

% Bathe, equation 6
% grz = (1/2)*(1-s)*grzC + (1/2)*(1+s)*grzA;
% gsz = (1/2)*(1-r)*gszB + (1/2)*(1+r)*gszD;

gxiz  = (1-eta)*gs4 -(2)^(1/2)*eta*gs5+eta*gs6;
getaz =      xi*gs4 +(2)^(1/2) *xi*gs5+(1-xi)*gs6;

gbar = [ gxiz; getaz ];

Ng = equationsToMatrix(gbar, g_xieta);
disp('Ngamma = '); disp(Ng)

%%
disp('*** DEMO KATILI, EQUATION 24 ***')
syms Bss Bnn b Db Ds Bsn v Bns dBsk Lk s Bni Bnj Bsi Bsj wi wj phik Ck Sk Bxi Bxj Byi Byj
syms w1 bx1 by1 w2 bx2 by2 w3 bx3 by3 w4 bx4 by4
syms x43 x32 x21 x13 y43 y32 y21 y13 L4 L5 L6 c4 c5 c6 s4 s5 s6
syms w1 w2 w3 w4 phik Ck Sk Bxi Bxj Byi Byj
syms phi4 phi5 phi6
syms S4 S5 S6 
syms C4 C5 C6  
syms Bx1 Bx2 Bx3 Bx4
syms By1 By2 By3 By4
syms L4 L5 L6  

%c4 = x21/L4;      s4 = y21/L4;
%c5 = x32/L5;      s5 = y32/L5;
%c6 = x13/L6;      s6 = y13/L6;

dBsk=0;
Un = [ w1 bx1 by1 w2 bx2 by2 w3 bx3 by3].';
Un1 = [ w1 Bx1 By1 w2 Bx2 By2 w3 Bx3 By3].';

Mss=Db*(Bss+v*Bnn);

Mns=Db*(1-v)/2*(Bsn+Bns);

Ts=Mss+Mns;
gamaaz=Ts/Ds;  %16
Bsi=Ck*Bxi+Sk*Byi;
Bsj=Ck*Bxj+Sk*Byj;
Bn=(1-s/Lk)*Bni+(s/Lk)*Bnj;
Bs=((1-s/Lk)*Bsi+(s/Lk)*Bsj)/Lk+4*s/Lk*(1-s/Lk)*dBsk;
ws=(-(1-s/Lk)*wi +(s/Lk)*wj)/Lk^2*2;
gamazk=-2/3*phik*dBsk;%%22a
gamask=Bs+ws;%%%32
ec=int(gamask-gamazk,s,0,Lk);%% 35b

gs4=subs(ec,{wj,wi,Bxj,Bxi,Byj,Byi,Lk,Ck,Sk},{w2,w1,bx2,bx1,By2,By1,L4,c4,s4});
gs5=subs(ec,{wj,wi,Bxj,Bxi,Byj,Byi,Lk,Ck,Sk},{w3,w2,bx3,bx2,By3,By2,L5,c5,s5});
gs6=subs(ec,{wj,wi,Bxj,Bxi,Byj,Byi,Lk,Ck,Sk},{w1,w3,bx1,bx3,By1,By3,L6,c6,s6});

%An = simplify(equationsToMatrix([ dbs4; dbs5; dbs6 ], Un1));

% eq 49
% gsk = (wj - wi)/Lk + (ck*bxi + sk*byi)/2 + (ck*bxj + sk*byj)/2;
%gs4   = (w2 - w1)/L4 + (c4*bx1 + s4*by1)/2 + (c4*bx2 + s4*by2)/2;
%gs5   = (w3 - w2)/L5 + (c5*bx2 + s5*by2)/2 + (c5*bx3 + s5*by3)/2;
%gs6   = (w1 - w3)/L6 + (c6*bx3 + s6*by3)/2 + (c6*bx1 + s6*by1)/2;

% Equation 26
Au = equationsToMatrix([ gs4; gs5; gs6 ], Un);
disp('Au = (1/2)*'); pretty(2*Au);

%%
disp('*** DEMO KATILI, EQUATION 18 ***')

% Equation 18
Ag = diag([ L4 L5/(2)^(1/2) -L6])

% Equation 52
Ng_Ag_Au = simplify(expand(Ng*Ag*Au));
disp('(Ng*Ag*Au)^T = '); Ng_Ag_Au.'

Bs = invJ*Ng_Ag_Au


%%
disp('************************************************************************')
disp('EQUATIONS Bs BATHE-DVORKIN')
clear
syms r s eta xi 
syms w1 tx1 ty1 w2 tx2 ty2 w3 tx3 ty3 w4 tx4 ty4 
syms x14 x21 x32 x43 y14 y21 y32 y43 x13 y13

%{
eq_grz = (1+s)*((w1 - w2)/2 + ((x1 - x2)/4)*(ty1 + ty2) - ((y1 - y2)/4)*(tx1 + tx2)) + ...
         (1-s)*((w4 - w3)/2 + ((x4 - x3)/4)*(ty4 + ty3) - ((y4 - y3)/4)*(tx4 + tx3));
eq_gsz = (1+r)*((w1 - w4)/2 + ((x1 - x4)/4)*(ty1 + ty4) - ((y1 - y4)/4)*(tx1 + tx4)) + ...
         (1-r)*((w2 - w3)/2 + ((x2 - x3)/4)*(ty2 + ty3) - ((y2 - y3)/4)*(tx2 + tx3));
%}
eq_grz = (1+s)*((w1 - w2)/2 - (x21/4)*(ty1 + ty2) + (y21/4)*(tx1 + tx2)) + ...
         (1-s)*((w4 - w3)/2 + (x43/4)*(ty4 + ty3) - (y43/4)*(tx4 + tx3));
eq_gsz = (1+r)*((w1 - w4)/2 + (x14/4)*(ty1 + ty4) - (y14/4)*(tx1 + tx4)) + ...
         (1-r)*((w2 - w3)/2 - (x32/4)*(ty2 + ty3) + (y32/4)*(tx2 + tx3));

ae = [ w1 tx1 ty1 w2 tx2 ty2 w3 tx3 ty3 ].';
Bgrz = simplify(expand(equationsToMatrix(eq_grz, ae))).'
Bgsz = simplify(expand(equationsToMatrix(eq_gsz, ae))).'


B1=[-1,1/2*(x21+x32*eta), 1/2*(y21+y32*eta);
    -1,1/2*(x13+x32*xi) ,-1/2*(y13+y32*xi)];  

B2=[ 1,1/2*(x21+x13*eta),1/2*(y21+y13*eta);
     0,-1/2*x13*xi      ,-1/2*y13*xi];  

B3=[ 0,1/2*x21*eta      ,1/2*y21*eta;
     1,-1/2*(x13+x21*xi),-1/2*(y13+y21*xi)]; 
j=1/(2*A)*[-y13,-y21;
            x13,x21];
gama=[B1,B3,B3];
Bs1=j*gama 




 

 