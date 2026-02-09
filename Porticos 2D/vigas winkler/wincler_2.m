clc
clear
syms V D k v q M EI Ac  t C1 C2 C3 C4 x C5 C6 xi AE L
%{
E = 24870062.3;                 % Modulo de elasticidad del concreto [kPa]
G = 0.4*E;                      % Modulo de cortante [kPa]
L = 4;                          % Longitud [m]
Ae = 0.4^2;                     % Area [m^2]
I = 0.4^4/12;                   % Inercia [m^4]
Ac = Ae*G*5/6;                  % Coeficiente de corrección por cortante
k = 500;                        % Coeficiente de balasto
EI = E*I;                       % Rigidez a flexión
AE = Ae*E;                      % Rigidez axial
rho = 2.4;                      % Densidad del concreto [Mg/m^3]
g = 9.8066502;                  % Aceleración de gravedad
P = 1;                          % Carga axial de pandeo [kN/m]
q = 0;                          % Carga vertical  o para matriz de rigidez[kN/m]
b = 0;                          % Carga axial o para matriz de rigidez[kN/m]
EA=E*Ae;
%}
%despejo todo en funcion de M
ec1=V*D == -k*v+q; %   (1) %Equilibrio de fuerzas cortantes
ec2=M*D == V      ;% (2)  %%Relación momento-cortante
ec3=t*D == M/EI   ;% (3)  %Relación giro-momento
ec4=v*D == t - V/Ac;%(4)  %Relación desplazamiento-giro

[V1,M1,t1,v1]=solve(ec1,ec2,ec3,ec4,[V,M,t,v]);

V1 =      (Ac*D^3*EI*q)/(Ac*EI*D^4 - EI*k*D^2 + Ac*k);
M1 =      (Ac*D^2*EI*q)/(Ac*EI*D^4 - EI*k*D^2 + Ac*k);
t1 =           (Ac*D*q)/(Ac*EI*D^4 - EI*k*D^2 + Ac*k);
v1 =(q*(- EI*D^2 + Ac))/(Ac*EI*D^4 - EI*k*D^2 + Ac*k);

s=solve((D^4 - k*D^2/Ac + k/EI),D);
%ss=double(s);

a=1;
b=-k/Ac;
c=k/EI;

u1=(-b + (b^2-4*a*c)^(1/2)) / 2;
u2=(-b - (b^2-4*a*c)^(1/2)) / 2;

r1=+(u1)^(1/2);
r2=-(u1)^(1/2);
r3=+(u2)^(1/2);
r4=-(u2)^(1/2);

%M=C1*exp(r1*x)+C2*exp(r2*x)+C3*exp(r3*x)+C4*exp(r4*x);

La1=(-u1)^(1/2);
La2=(-u2)^(1/2);

syms La1 La2


M=C1*cos(La1*x)+C2*sin(La1*x)+C3*cos(La2*x)+C4*sin(La2*x);

V=diff(M,x,1);
v=-diff(V,x,1)/k;
t=V/Ac+diff(v,x,1);
bp=0;
%se definen las ecuaciones diferenciales a carga axial
A=int(bp,x)+C5;
u=int(A/AE,x)+C6;

%# Se calcula la matrix de rigidez
K_TE2 = sym(zeros(6,6));
N_u2 = sym(zeros(1,6));
N_w2 = sym(zeros(1,6));
N_t2 = sym(zeros(1,6));

for i = 1:6
    [c1,c2,c3,c4,c5,c6]=solve(subs(u,x,0)==(i==1),...
                              subs(v,x,0)==(i==2),...% con sus respectivas condiciones de frontera
                              subs(t,x,0)==(i==3),...
                              subs(u,x,L)==(i==4),...
                              subs(v,x,L)==(i==5),...
                              subs(t,x,L)==(i==6),...
                              [C1,C2,C3,C4,C5,C6]);
    % # se evaluan las reacciones horizontales y verticales y los momentos en los apoyos
	K_TE2(:,i)=[-subs(A,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,0}); % X1
                 subs(V,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,0}); % Y2
                -subs(M,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,0}); % M2
                 subs(A,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L}); % X2 
                -subs(V,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L}); % Y2
                 subs(M,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L})];% M2
	N_t2(i) = subs(t,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L*(1+xi)/2});
	N_w2(i) = subs(v,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L*(1+xi)/2});
	N_u2(i) = subs(u,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L*(1+xi)/2}); 
end
%K_TE2=	double(real(K_TE2));



K_TE2=	simplify((K_TE2));



                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        (k*(La1^2 - La2^2)*(Ac*La1^2 + k)*(Ac*La2^2 + k)*(cos(L*La1) - cos(L*La2)))/(2*Ac^2*La1^4*La2^4*cos(L*La1)*cos(L*La2) - 2*La1^2*La2^2*k^2 - 2*Ac*La1^2*La2^4*k - 2*Ac*La1^4*La2^2*k - 2*Ac^2*La1^4*La2^4 + Ac^2*La1^3*La2^5*sin(L*La1)*sin(L*La2) + Ac^2*La1^5*La2^3*sin(L*La1)*sin(L*La2) + 2*La1^2*La2^2*k^2*cos(L*La1)*cos(L*La2) + La1*La2^3*k^2*sin(L*La1)*sin(L*La2) + La1^3*La2*k^2*sin(L*La1)*sin(L*La2) + 2*Ac*La1^2*La2^4*k*cos(L*La1)*cos(L*La2) + 2*Ac*La1^4*La2^2*k*cos(L*La1)*cos(L*La2) + 4*Ac*La1^3*La2^3*k*sin(L*La1)*sin(L*La2)),                                                                                                                                  (Ac*k*(La1^2 - La2^2)*(La1*k*sin(L*La1) - La2*k*sin(L*La2) + Ac*La1*La2^2*sin(L*La1) - Ac*La1^2*La2*sin(L*La2)))/(2*Ac^2*La1^4*La2^4*cos(L*La1)*cos(L*La2) - 2*La1^2*La2^2*k^2 - 2*Ac*La1^2*La2^4*k - 2*Ac*La1^4*La2^2*k - 2*Ac^2*La1^4*La2^4 + Ac^2*La1^3*La2^5*sin(L*La1)*sin(L*La2) + Ac^2*La1^5*La2^3*sin(L*La1)*sin(L*La2) + 2*La1^2*La2^2*k^2*cos(L*La1)*cos(L*La2) + La1*La2^3*k^2*sin(L*La1)*sin(L*La2) + La1^3*La2*k^2*sin(L*La1)*sin(L*La2) + 2*Ac*La1^2*La2^4*k*cos(L*La1)*cos(L*La2) + 2*Ac*La1^4*La2^2*k*cos(L*La1)*cos(L*La2) + 4*Ac*La1^3*La2^3*k*sin(L*La1)*sin(L*La2)),     0, (La1^2*k^3*cos(L*La1)*cos(L*La2) - La2^2*k^3 - Ac*La1^4*k^2 - Ac*La2^4*k^2 - La1^2*k^3 + La2^2*k^3*cos(L*La1)*cos(L*La2) - 2*Ac*La1^2*La2^2*k^2 - Ac^2*La1^2*La2^4*k - Ac^2*La1^4*La2^2*k + 2*La1*La2*k^3*sin(L*La1)*sin(L*La2) + Ac*La1^4*k^2*cos(L*La1)*cos(L*La2) + Ac*La2^4*k^2*cos(L*La1)*cos(L*La2) + 2*Ac*La1*La2^3*k^2*sin(L*La1)*sin(L*La2) + 2*Ac*La1^3*La2*k^2*sin(L*La1)*sin(L*La2) + Ac^2*La1*La2^5*k*sin(L*La1)*sin(L*La2) + Ac^2*La1^5*La2*k*sin(L*La1)*sin(L*La2) + 2*Ac*La1^2*La2^2*k^2*cos(L*La1)*cos(L*La2) + Ac^2*La1^2*La2^4*k*cos(L*La1)*cos(L*La2) + Ac^2*La1^4*La2^2*k*cos(L*La1)*cos(L*La2))/(2*Ac^2*La1^4*La2^4*cos(L*La1)*cos(L*La2) - 2*La1^2*La2^2*k^2 - 2*Ac*La1^2*La2^4*k - 2*Ac*La1^4*La2^2*k - 2*Ac^2*La1^4*La2^4 + Ac^2*La1^3*La2^5*sin(L*La1)*sin(L*La2) + Ac^2*La1^5*La2^3*sin(L*La1)*sin(L*La2) + 2*La1^2*La2^2*k^2*cos(L*La1)*cos(L*La2) + La1*La2^3*k^2*sin(L*La1)*sin(L*La2) + La1^3*La2*k^2*sin(L*La1)*sin(L*La2) + 2*Ac*La1^2*La2^4*k*cos(L*La1)*cos(L*La2) + 2*Ac*La1^4*La2^2*k*cos(L*La1)*cos(L*La2) + 4*Ac*La1^3*La2^3*k*sin(L*La1)*sin(L*La2)),                                                                                      -(Ac*k*(La1^2 - La2^2)*(La1*k*cos(L*La2)*sin(L*La1) - La2*k*cos(L*La1)*sin(L*La2) + Ac*La1*La2^2*cos(L*La2)*sin(L*La1) - Ac*La1^2*La2*cos(L*La1)*sin(L*La2)))/(2*Ac^2*La1^4*La2^4*cos(L*La1)*cos(L*La2) - 2*La1^2*La2^2*k^2 - 2*Ac*La1^2*La2^4*k - 2*Ac*La1^4*La2^2*k - 2*Ac^2*La1^4*La2^4 + Ac^2*La1^3*La2^5*sin(L*La1)*sin(L*La2) + Ac^2*La1^5*La2^3*sin(L*La1)*sin(L*La2) + 2*La1^2*La2^2*k^2*cos(L*La1)*cos(L*La2) + La1*La2^3*k^2*sin(L*La1)*sin(L*La2) + La1^3*La2*k^2*sin(L*La1)*sin(L*La2) + 2*Ac*La1^2*La2^4*k*cos(L*La1)*cos(L*La2) + 2*Ac*La1^4*La2^2*k*cos(L*La1)*cos(L*La2) + 4*Ac*La1^3*La2^3*k*sin(L*La1)*sin(L*La2))]
 