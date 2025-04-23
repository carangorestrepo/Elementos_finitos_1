clc
clear
syms eta_gl xi_gl
syms x1 x2 x3 x4 y1 y2 y3 y4
syms C4 C5 C6 C7 C8
syms S4 S5 S6 S7 S8
syms A1 A2 A3 A4
syms L4 L5 L6 L7 L8 
syms x21 x32 x13  y21 y32 y13
syms x21 x32 x13 x43  x14 y21  y32 y13 y43 y14


%xi_gl=-0.577350269189626;
%eta_gl=-0.577350269189626;

%% Derivadas de N con respecto a xi   
dN_dxi = [-1/4*(1-eta_gl)             % dN1_dxi
          1/4*(1-eta_gl)              % dN2_dxi
          1/4*(1+eta_gl)              % dN3_dxi
         -1/4*(1+eta_gl)    ];        % dN4_dxi
%% Derivadas de N con respecto a eta    
dN_deta = [-1/4*(1-xi_gl)             % dN1_deta
          -1/4*(1+xi_gl)              % dN2_deta
           1/4*(1+xi_gl)              % dN3_deta
           1/4*(1-xi_gl)    ];        % dN4_deta 
       
%% Derivadas de P con respecto a xi
dP_dxi = [-xi_gl*(1-eta_gl)              % dP1_dxi
         1/2*(1-eta_gl^2)             % dP2_dxi
         -xi_gl*(1+eta_gl)               % dP3_dxi
         -1/2*(1-eta_gl^2)  ];        % dP4_dxi
                        
%% Derivadas de P con respecto a eta    
dP_deta =[-1/2*(1-xi_gl^2)           % dP1_deta
         -eta_gl*(1+xi_gl)              % dP2_deta
          1/2*(1-xi_gl^2)            % dP3_deta
         -eta_gl*(1-xi_gl)    ];        % dP4_deta     

%xe=[0;0.0500000000000000;0.0500000000000000;0];
%ye=[0;0;0.0500000000000000;0.0500000000000000];
%x21 = xe(2) - xe(1);         y21 = ye(2) - ye(1); 
%x32 = xe(3) - xe(2);         y32 = ye(3) - ye(2);
%x43 = xe(4) - xe(3);         y43 = ye(4) - ye(3);    
%x14 = xe(1) - xe(4);         y14 = ye(1) - ye(4);

xe=[x1;x2;x3;x4];
ye=[y1;y2;y3;y4];
xji = [ x21 x32 x43 x14 ];   yji = [ y21 y32 y43 y14 ];   
Lk = hypot(xji, yji);      Ck =xji./Lk;      Sk = yji./Lk; %% figure 4
%% Matriz jacobiana, su inversa y determinante
% Se ensambla la matriz jacobiana
dx_dxi  = sum(dN_dxi .*xe);   dy_dxi  = sum(dN_dxi .*ye);
dx_deta = sum(dN_deta.*xe);   dy_deta = sum(dN_deta.*ye);
% Se calcula su inversa
Je = [ dx_dxi    dy_dxi
       dx_deta   dy_deta ];
   
% Se calcula su inversa
inv_Je = inv(Je);
j11 = inv_Je(1,1);              j12 = inv_Je(1,2);
j21 = inv_Je(2,1);              j22 = inv_Je(2,2);                       
% y su determinante (el Jacobiano)
det_Je = det(Je);

%% Se calcula Bb_beta (ecuacion 12)
Bb_beta = sym(zeros(3,12));
%Bb_beta =zeros(3,4);
for i = 1:4                
    dNi_dx = j11*dN_dxi(i) + j12*dN_deta(i); % = ai
    dNi_dy = j21*dN_dxi(i) + j22*dN_deta(i); % = bi               
    Bb_beta(:,[3*i-2 3*i-1 3*i]) = [ 0    dNi_dx         0
                                     0         0    dNi_dy
                                     0    dNi_dy    dNi_dx ];
end
Bb_beta=simplify(Bb_beta);
%% Se calcula Bb_dbeta (ecuacion 13)
Bb_dbeta = sym(zeros(3,4));            
for k = 1:4            
    dPk_dx = j11*dP_dxi(k) + j12*dP_deta(k);
    dPk_dy = j21*dP_dxi(k) + j22*dP_deta(k);
    Bb_dbeta(:,k) = [ dPk_dx*Ck(k)
                      dPk_dy*Sk(k)
                      dPk_dy*Ck(k) + dPk_dx*Sk(k) ];
end   

Bb_dbeta=simplify(Bb_dbeta);
a=1