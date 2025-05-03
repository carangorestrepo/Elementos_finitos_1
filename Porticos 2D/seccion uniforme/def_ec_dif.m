function [faxial,V,M,u,v] = def_ec_dif(L, EI, Ac, AE, b1loc,b2loc, q1loc,q2loc,k,P,v1,v2,t1,t2,u1,u2,puntos_graficas)

%{
syms x L q1 q2 b1 b2 m b
x1 = 0;
x2 = L;
sol = solve(q1 == m*x1+b, q2 == m*x2+b, m, b);
q = simplify(sol.m*x + sol.b)
sol = solve(b1 == m*x1+b, b2 == m*x2+b, m, b);
b = simplify(sol.m*x + sol.b)
%}
qxloc = @(x) b1loc - (x*(b1loc - b2loc))/L;
qyloc = @(x) q1loc - (x*(q1loc - q2loc))/L;

%% se definen algunas constantes
v_ = 1; t_ = 2; M_ = 3; V_ = 4; u_ = 5; fax_ = 6;

npuntos = puntos_graficas;
xinit   = linspace(0, L, npuntos);
solini  = bvpinit(xinit, zeros(6,1));
sol     = bvp5c(@ecuacion_diferencial, @condiciones_de_apoyo, solini);

%% Calculos intermedios
y = deval(sol, xinit);
faxial = y(6,:);           % Fuerza axial [kN]
V      = y(4,:);           % Fuerza cortante [kN]
M      = y(3,:);           % Momento flector [kN/m]
u     = y(5,:);           % Desplazamiento horizontal de la viga [m]
v     = y(1,:);           % Desplazamiento vertical de la viga [m]
%theta = atan(sol.y(2,:)); % Angulo de giro  [rad]
%X1 = +faxial(1);   Y1 = -V(1);   M1 = +M(1);   % 1   => en x=0
%X2 = -faxial(end); Y2 = +V(end); M2 = -M(end); % end => en x=L

%fe = [ X1; Y1; M1; X2; Y2; M2 ];

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
      dydx(v_) = y(t_)-y(V_)/Ac;      % = theta
      dydx(t_) = y(M_)/(EI);          % = M/(EI)
      dydx(M_) = y(V_) - (P)*dydx(v_);% = V-P*diff(v)
      dydx(V_) = qyloc(x)-k*y(v_);    % = qyloc
      dydx(u_) = y(fax_)/(AE);        % = u
      dydx(fax_) = -qxloc(x);         % = faxial
   end

%% ------------------------------------------------------------------------

   function res = condiciones_de_apoyo(YL,YR)
      % condiciones de apoyo (cond. de frontera de la ecuacion diferencial)
      v_  = 1; t_ = 2; M_ = 3; V_ = 4; u_ = 5; fax_ = 6;
      res = [ % YL: apoyo izquierdo (LEFT), YR: apoyo derecho (RIGHT)
              YL(u_) - u1          % uloc(0)     = 0
              YL(v_) - v1          % vloc(0)     = 0
              YL(t_) - t1          % thetaloc(0) = 0
              YR(u_) - u2          % uloc(L)     = 0
              YR(v_) - v2          % vloc(L)     = 0
              YR(t_) - t2 ];       % thetaloc(L) = 0
   end
end
%% bye, bye!