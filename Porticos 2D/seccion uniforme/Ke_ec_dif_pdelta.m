function Ke = Ke_ec_dif_pdelta(L,EI, AE,Ac,k,P,puntos_graficas)

%% se definen algunas constantes
v_ = 1; t_ = 2; M_ = 3; V_ = 4; u_ = 5; fax_ = 6;

npuntos = puntos_graficas;
xinit   = linspace(0, L, npuntos);
solini  = bvpinit(xinit, zeros(6,1));
Ke = zeros(6);
for i = 1:6
    ae = [i==1, i==2, i==3, i==4, i==5, i==6];
    sol = bvp5c(@ecuacion_diferencial, @condiciones_de_apoyo, solini);
    Ke(:,i) = [ -sol.y(fax_,   1)
                +sol.y(V_,     1)
                -sol.y(M_,     1)
                +sol.y(fax_, end)
                -sol.y(V_,   end)
                +sol.y(M_,   end) ];
end
%{
L2=L^2; 
L3=L^3;
AE = A*E;
EI = E*I;
Ke2 = [ AE/L   0         0        -AE/L    0          0       
        0     12*EI/L3   6*EI/L2   0     -12*EI/L3   6*EI/L2
        0      6*EI/L2   4*EI/L    0      -6*EI/L2   2*EI/L
       -AE/L   0         0         AE/L    0         0
        0    -12*EI/L3  -6*EI/L2   0      12*EI/L3  -6*EI/L2
        0      6*EI/L2   2*EI/L    0      -6*EI/L2   4*EI/L];
% error entre respuesta anal�tica y soluci�n por ecuaci�n diferencial:
max(max(abs(Ke - Ke2)))
%}

%% se define la ecuacion diferencial asociada
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
      %            y(v_)           = v
      dydx(v_)   = y(t_)-(dydx(M_)+(P)*dydx(v_))/Ac;      % = theta
      dydx(t_)   = y(M_)/(EI);          % = M/(EI)
      dydx(M_)   = y(V_) - (P)*dydx(v_);% = V-P*diff(v)
      dydx(V_)   = 0-k*y(v_);           % = qyloc
      dydx(u_)   = y(fax_)/(AE);        % = u
      dydx(fax_) = 0;                   % = faxial
   end

%% se definen las condiciones de frontera de la ecuacion diferencial
   function res = condiciones_de_apoyo(YL,YR)
      % condiciones de apoyo (cond. de frontera de la ecuacion diferencial)
      u1  = 1; v1 = 2; t1 = 3; u2 = 4; v2 = 5; t2   = 6;
      res = [ % YL: apoyo izquierdo (LEFT), YR: apoyo derecho (RIGHT)
              YL(u_) - ae(u1)          % uloc(0)     = u1
              YL(v_) - ae(v1)          % vloc(0)     = v1
              YL(t_) - ae(t1)          % thetaloc(0) = t1
              YR(u_) - ae(u2)          % uloc(L)     = u2
              YR(v_) - ae(v2)          % vloc(L)     = v2
              YR(t_) - ae(t2) ];       % thetaloc(L) = t2
   end
end
%% bye, bye!