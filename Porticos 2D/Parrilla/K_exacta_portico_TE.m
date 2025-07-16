clc
clear
% # Definción de variables
syms x L A(x) V(x) M(x) t(x) v(x) u(x) EI EA GJ Ac b  xi k b2 b1 wb wa Ae I g qa qb;
syms va vb ta tb ua ub
%{
%% Se describen las propiedades de los materiales
g=9.8066502;      % aceleracion de la gravedad
E = 24870062.3;       % [kPa] modulo de elasticidad del concreto (NSR 10  - art C.8.5)
G = 0.4*E;                    % [kPa] modulo de cortante (NSR 10  - art CR8.8.2)
L = 4;                        % [m]   longitud
Ae = 0.4^2;                   % [m^2] area  
I = 0.4^4/12;                 % [m^4] inercia
Ac=Ae*G*5/6;                  % coeficiente de correccion del cortante para seccion rectangular
%}
%% # Se define la carga: hace la ecuación diferencial homogénea
q =0;  
m = 0;
%EA=Ae*E;
%EI=E*I;
%# Se calcula la matrix de rigidez

%EI=E*I;
%EA=E*Ae;

K_TE2 = sym(zeros(6));
N_u2 = sym(zeros(1,6));
N_w2 = sym(zeros(1,6));
N_t2 = sym(zeros(1,6));
for i = 1:6
    sol = dsolve(...       
           diff(V)== q,        ... % diff(V,x)+k*v==q se definen las ecuaciones diferenciales
           diff(M) == V - m,...%-P*diff(v,x),    ... %-P*diff(v,x),    ...
           diff(t) == M/(EI),   ... 
           diff(v) == t- V/Ac , ...- V/Ac
           diff(A) == 0,        ...
           diff(u) == A/(GJ),     ...  
           t(0) == (i==1),      ...
           u(0) == (i==2),      ... % con sus respectivas condiciones de frontera
           v(0) == (i==3),      ...
           t(L) == (i==4),      ...
           u(L) == (i==5),      ...   
           v(L) == (i==6)); 
   N_t2(i) = sol.t;
   N_w2(i) = sol.v;
   N_u2(i) = sol.u; 
   % # se evaluan las reacciones horizontales y verticales y los momentos en los apoyos
       K_TE2(:,i)=[-subs(sol.M,x,0); % M2
                   -subs(sol.A,x,0); % X1
                    subs(sol.V,x,0); % Y2
                    subs(sol.M,x,L);% M2
                    subs(sol.A,x,L); % X2 
                   -subs(sol.V,x,L)]; % Y2
end
K_TE2=(K_TE2);
%# Se convierten las funciones de forma a coordenadas naturales,

F_TE2 = sym(zeros(6,1));
q=qa;
%q=(q2-q1)/L*x+q1;
b=b1;
%b=(b2-b1)/L*x+b1;
sol = dsolve(...       
           diff(V)== q,        ... % diff(V,x)+k*v==q se definen las ecuaciones diferenciales
           diff(M) == V - m,...%-P*diff(v,x),    ... %-P*diff(v,x),    ...
           diff(t) == M/(EI),   ... 
           diff(v) == t- V/Ac, ...- V/Ac
           diff(A) == b,        ...
           diff(u) == A/(GJ),     ...  
           t(0) == ta,      ...
           u(0) == ua,      ... % con sus respectivas condiciones de frontera
           v(0) == va,      ...
           t(L) == tb,      ...
           u(L) == ub,      ...   
           v(L) == vb); 

% # se evaluan las reacciones horizontales y verticales y los momentos en los apoyos
F_TE2(:,1)=[subs(sol.M,x,0); % M2
           subs(sol.A,x,0); % X1
            subs(sol.V,x,0); % Y2
            subs(sol.M,x,L);% M2
            subs(sol.A,x,L); % X2 
           -subs(sol.V,x,L)]; % Y2
       
a=1