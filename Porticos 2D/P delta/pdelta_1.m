clc
clear
syms D x C1 C2 C3 C4 w Q V M v EI Ac P t


ec1=V == Q - P*D*v;          % (1) Relación fuerza cortante
ec2=D*M == Q;                % (2) Equilibrio de momentos
ec3=D*V == -w;               % (3) Equilibrio de fuerzas verticales
ec4=M == EI*D*t;             % (4) Relación momento-curvatura
ec5=Q ==Ac*(t + D*v);    % (5) Relación fuerza transversal-deformación

[Q1,V1,M1,t1,v1]=solve(ec1,ec2,ec3,ec4,ec5,[Q,V,M,t,v]);
 
Q1 =-(Ac*D*EI*w)/(Ac*P + Ac*D^2*EI - D^2*EI*P)
V1 =-w/D
M1 =-(Ac*EI*w)/(Ac*P + Ac*D^2*EI - D^2*EI*P)
t1 =-(Ac*w)/(D*(Ac*P + Ac*D^2*EI - D^2*EI*P))
v1 =(w*(- EI*D^2 + Ac))/(D*(Ac*D*P + Ac*D^3*EI - D^3*EI*P))
 

ec5=D*M ==Ac*(t + D*v);    % (5) Relación fuerza transversal-deformación
ec4=D*t == M/EI;             % (4) Relación momento-curvatura
ec1=D*v == (D*M - V)/P;          % (1) Relación fuerza cortante
ec3=D*V == -w;               % (3) Equilibrio de fuerzas verticales


ec5=D*v ==t + D*M/Ac;    % (5) Relación fuerza transversal-deformación
ec4=D*t == M/EI;             % (4) Relación momento-curvatura
ec1=D*M == P*D*v - V;          % (1) Relación fuerza cortante
ec3=D*V == -w;               % (3) Equilibrio de fuerzas verticales


[Q1,V1,M1,t1,v1]=solve(ec1,ec3,ec4,ec5,[Q,V,M,t,v]);
 
Q1 =-(Ac*D*EI*w)/(Ac*P + Ac*D^2*EI - D^2*EI*P)

V1 =-w/D
M1 =-(Ac*EI*w)/(Ac*P + Ac*D^2*EI - D^2*EI*P)
t1 =-(Ac*w)/(D*(Ac*P + Ac*D^2*EI - D^2*EI*P))
v1 =(w*(- EI*D^2 + Ac))/(D*(Ac*D*P + Ac*D^3*EI - D^3*EI*P))

  dydx(v_) = y(t_)-(dydx(M_)+(P)*dydx(v_))/Ac;      % = theta
  dydx(t_) = y(M_)/(EI);          % = M/(EI)
  dydx(M_) = y(V_) - (P)*dydx(v_);% = V-P*diff(v)
  dydx(V_) = qyloc(x)-k*y(v_);    % = qyloc
  dydx(u_) = y(fax_)/(AE);        % = u
  dydx(fax_) = -qxloc(x);         % = faxial


 