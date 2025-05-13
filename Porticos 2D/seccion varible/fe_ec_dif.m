function [X1,Y1,M1,X2,Y2,M2] = fe_ec_dif(tf,tw,nb1loc,nq1loc,nb,nh,L,E,G,b2,b1,h2,h1,sec, b2loc,b1loc, q2loc,q1loc,k,P,puntos_graficas)

%{
syms x L q1 q2 b1 b2 m b
x1 = 0;
x2 = L;
sol = solve(q1 == m*x1+b, q2 == m*x2+b, m, b);
q = simplify(sol.m*x + sol.b)
sol = solve(b1 == m*x1+b, b2 == m*x2+b, m, b);
b = simplify(sol.m*x + sol.b)
%}


qxloc = @(x) b1loc - (x^nb1loc*(b1loc - b2loc))/L^nb1loc;
qyloc = @(x) q1loc - (x^nq1loc*(q1loc - q2loc))/L^nq1loc;
%I =  @(x)       ((b1 - (x^nb*(b1 - b2))/L^nb)*(h1 - (x^nh*(h1 - h2))/L^nh)^3)/12;
%Ae = @(x)       ((b1 - (x^nb*(b1 - b2))/L^nb)*(h1 - (x^nh*(h1 - h2))/L^nh));
%Ac = @(x)(As2*G* (b1 - (x^nb*(b1 - b2))/L^nb)*(h1 - (x^nh*(h1 - h2))/L^nh));

[Ae,I,Ac]=secciones(sec,b1,b2,nb,h1,h2,nh,G,tf,tw,L);


%% se definen algunas constantes
v_ = 1; t_ = 2; M_ = 3; V_ = 4; u_ = 5; fax_ = 6;

npuntos = puntos_graficas;
xinit   = linspace(0, L, npuntos);
solini  = bvpinit(xinit, zeros(6,1));
sol     = bvp5c(@ecuacion_diferencial, @condiciones_de_apoyo, solini);

%% Calculos intermedios
y = deval(sol, [0 L]);
faxial = y(6,:);           % Fuerza axial [kN]
V      = y(4,:);           % Fuerza cortante [kN]
M      = y(3,:);           % Momento flector [kN/m]
%u     = y(5,:);           % Desplazamiento horizontal de la viga [m]
%v     = y(1,:);           % Desplazamiento vertical de la viga [m]
%theta = atan(sol.y(2,:)); % Angulo de giro  [rad]
X1 = +faxial(1);   Y1 = -V(1);   M1 = +M(1);   % 1   => en x=0
X2 = -faxial(end); Y2 = +V(end); M2 = -M(end); % end => en x=L

fe = [ X1; Y1; M1; X2; Y2; M2 ];

%% -----------------------------------------------------------------------
   function dydx = ecuacion_diferencial(x,y)
      % aqui se implementa la ecuacion diferencial para vigas de material
      % homogeneo y seccion transversal constante (A, E, I, qx, qy las 
      % provee la funcion exterior)
      %      d^4 v(x)
      % E I ---------- = q(x)
      %        dx^4
      %
      %      d^2 u(x)
      % A E ---------- = -b(x)
      %        dx^2

      dydx = zeros(6,1);
      %         y(1)          = v
      dydx(v_) = y(t_)-y(V_)/Ac(x);       % = theta
      dydx(t_) = y(M_)/(E*I(x)); % = M/(EI)
      dydx(M_) = y(V_)-(P)*dydx(v_);       % = V
      dydx(V_) = qyloc(x)-k*y(v_);   % = qyloc
      dydx(u_) = y(fax_)/(Ae(x)*E); % = u
      dydx(fax_) = -qxloc(x);  % = faxial
   end

%% ------------------------------------------------------------------------

   function res = condiciones_de_apoyo(YL,YR)
      % condiciones de apoyo (cond. de frontera de la ecuacion diferencial)
      v_  = 1; t_ = 2; M_ = 3; V_ = 4; u_ = 5; fax_ = 6;
      res = [ % YL: apoyo izquierdo (LEFT), YR: apoyo derecho (RIGHT)
              YL(u_) - 0          % uloc(0)     = 0
              YL(v_) - 0          % vloc(0)     = 0
              YL(t_) - 0          % thetaloc(0) = 0
              YR(u_) - 0          % uloc(L)     = 0
              YR(v_) - 0          % vloc(L)     = 0
              YR(t_) - 0 ];       % thetaloc(L) = 0
   end
end
%% bye, bye!