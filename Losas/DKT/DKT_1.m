clc
clear
x=[0,1,0];
y=[0,0,1];

plot(x,y)

syms x23 x31 x12 y23 y31 y12 xi eta v E t;
%{
E  = 4.0*4880.0;       % [Pa]    modulo de elasticidad = 210GPa
v = 1/3;         %         coeficiente de Poisson
t  = 1;        % [m]     espesor de la losa
x23=x(2)-x(3);x31=x(3)-x(1);x12=x(1)-x(2);
y23=y(2)-y(3);y31=y(3)-y(1);y12=y(1)-y(2);
Db=E*t^3/(12*(1-v))*[1,v,0;
                    v,1,0;
                    0,0,(1-v)/2];
%}
                    
N=[2*(1-xi-eta)*(1/2-xi-eta);
   xi*(2*xi-1);
   eta*(2*eta-1);
   4*xi*eta;
   4*eta*(1-xi-eta);
   4*xi *(1-xi-eta)];

xij=[x23;x31;x12];
 
yij=[y23;y31;y12];
 
lij2=(xij.^2+yij.^2);
ak=-xij./lij2;
bk=-3/4*xij.*yij./lij2;
dk=-yij./lij2;
ek=(1/4*yij.^2-1/2*xij.^2)./lij2;
ck=(1/4*xij.^2-1/2*yij.^2)./lij2;
ek=[0;0;0;ek];
ak=[0;0;0;ak];
bk=[0;0;0;bk];
ck=[0;0;0;ck];
dk=[0;0;0;dk];

Hx=[1.5*(ak(6)*N(6)-ak(5)*N(5));
         bk(5)*N(5)+bk(6)*N(6);
    N(1)-ck(5)*N(5)-ck(6)*N(6);
    1.5*(ak(4)*N(4)-ak(6)*N(6));
         bk(6)*N(6)+bk(4)*N(4);
    N(2)-ck(6)*N(6)-ck(4)*N(4);
    1.5*(ak(5)*N(5)-ak(4)*N(4));
         bk(4)*N(4)+bk(5)*N(5);  
    N(3)-ck(4)*N(4)-ck(5)*N(5)];


Hx1=[1.5*(ak(6)*N(6) - ak(5)*N(5))
    bk(5)*N(5) + bk(6)*N(6)
    N(1) - ck(5)*N(5) - ck(6)*N(6)
    1.5*(ak(4)*N(4) - ak(6)*N(6))
    bk(4)*N(4) + bk(6)*N(6)
    N(2) - ck(4)*N(4) - ck(6)*N(6)
    1.5*(ak(5)*N(5) - ak(4)*N(4))
    bk(4)*N(4) + bk(5)*N(5)
    N(3) - ck(4)*N(4) - ck(5)*N(5)];
  
Hy=[1.5*(dk(6)*N(6)-dk(5)*N(5))
   -N(1)+ek(5)*N(5)-ek(6)*N(6)
        -bk(5)*N(5)-bk(6)*N(6)
    1.5*(dk(4)*N(4)-dk(6)*N(6))
   -N(2)+ek(6)*N(6)+ek(4)*N(4)    
        -bk(6)*N(6)-bk(4)*N(4)
    1.5*(dk(5)*N(5)-dk(4)*N(4))
   -N(3)+ek(4)*N(4)+ek(5)*N(5)
        -bk(4)*N(4)-bk(5)*N(5)];
    
    
Hy1 = [1.5*(dk(6)*N(6) - dk(5)*N(5));
       -N(1) + ek(5)*N(5) + ek(6)*N(6);
       -bk(5)*N(5) - bk(6)*N(6);
       1.5*(dk(4)*N(4) - dk(6)*N(6));
       -N(2) + ek(4)*N(4) + ek(6)*N(6);
       -bk(4)*N(4) - bk(6)*N(6);
       1.5*(dk(5)*N(5) - dk(4)*N(4));
       -N(3) + ek(4)*N(4) + ek(5)*N(5);
       -bk(4)*N(4) - bk(5)*N(5)];
dHxdxi=diff(Hx1,xi);
dHydxi=diff(Hy1,xi);

dHxdeta=diff(Hx1,eta);
dHydeta=diff(Hy1,eta);      
A2=x31*y12-x12*y31;
B = 1.0/(x31*y12 - x12*y31) * [y31*dHxdxi' + y12*dHxdeta';
                               -x31*dHydeta' - x12*dHydeta';
    -x31*dHxdxi'- x12*dHxdeta' + y31*dHydxi' + y12*dHydeta'];
K=double(A2*int(int(B'*Db*B,xi,0,1-eta),eta,0,1));
     
     
     