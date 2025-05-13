%% Limpieza inicial
%clear, clc, close all 

%% Definición de constantes para mejor legibilidad
X = 1; Y = 2; Z = 3;   % Coordenadas
ww = 1; tx = 2; ty = 3; % Grados de libertad por nodo

%% Definimos la geometria de la losa
nno = size(xnod,1);  % Número de nodos
ngdl = 3*nno;        % Grados de libertad totales
gdl = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % Asignación de GDL
nef = size(LaG,1);   % Número de elementos

%% Visualización de la malla
figure;
hold on;
for e = 1:nef
   plot(xnod(LaG(e,[1 2 3 4 1]), X), xnod(LaG(e,[1 2 3 4 1]), Y), 'b-');
   cgx = mean(xnod(LaG(e,:), X));
   cgy = mean(xnod(LaG(e,:), Y));
   text(cgx, cgy, num2str(e), 'Color', 'r');
end
plot(xnod(:,X), xnod(:,Y), 'r*');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
title('Malla de elementos finitos');

%% Matrices constitutivas
Db = (E*h^3/(12*(1-nu^2)));   % plate rigidity eq 7
H_b = Db * [ 1  nu 0          % matriz constitutiva de flexion generalizada
            nu 1  0           % (Dbe en la nomenclatura del curso) 
            0  0  (1-nu)/2 ]; 
D_s = k * h * G * eye(2);
D = blkdiag(H_b, D_s);
%T = rho * diag([h, h^3/12, h^3/12]);
T  = rho * diag([h, 0, 0]);

%% Funciones de forma y derivadas para cuadrilátero bilineal
Nforma = @(xi,eta) [1/4*(1-xi)*(1-eta)        % N1
                    1/4*(1+xi)*(1-eta)        % N2
                    1/4*(1+xi)*(1+eta)        % N3
                    1/4*(1-xi)*(1+eta)];      % N4

%% Derivadas de N con respecto a xi    
dN_dxi = @(xi,eta) [-1/4*(1-eta)              % dN1_dxi
                     1/4*(1-eta)              % dN2_dxi
                     1/4*(1+eta)              % dN3_dxi
                    -1/4*(1+eta)    ];        % dN4_dxi
                        
%% Derivadas de N con respecto a eta    
dN_deta = @(xi,eta) [-1/4*(1-xi)              % dN1_deta
                     -1/4*(1+xi)              % dN2_deta
                      1/4*(1+xi)              % dN3_deta
                      1/4*(1-xi)    ];        % dN4_deta     

%% Matrices [P] y [Q]
P = @(xi,eta) [...
                1,-xi, -eta,  -xi^2/2, -(xi^2*eta)/2 + eta^3/12, -eta^2/2,-(xi*eta^2)/2 + xi^3/12,-xi*eta/2, xi/2, xi*eta/2,eta/2, xi*eta/2;
                0,1  , 0   ,       xi,                   xi*eta,        0,       xi^2/4 + eta^2/2,    eta/2,  1/2,    eta/2,    0, -eta/2;
                0,0  , 1   ,        0,         eta^2/4 + xi^2/2,      eta,                 xi*eta,     xi/2,    0,    -xi/2,  1/2, xi/2];
Q = @(xi,eta) [...
                0,0,0,1,eta  ,0 ,xi/2 ,0,0,  0,0,0;
                0,0,0,0,eta/2,1 ,xi   ,0,0,  0,0,0;
                0,0,0,0,2*xi ,0 ,2*eta,1,0,  0,0,0;
                0,0,0,0,0    ,0 ,0    ,0,1,eta,0,0;
                0,0,0,0,0    ,0 ,0    ,0,0,  0,1,xi];

%% Integración numérica (Gauss 2x2)
n_gl = 2;
%[x_gl, w_gl] = gauss_points(n_gl);
[x_gl, w_gl] = gausslegendre_quad(n_gl);

%% Ensamblaje de matrices globales
K = zeros(ngdl, ngdl);
M = zeros(ngdl, ngdl);

for e = 1:nef
    xe = xnod(LaG(e,:), X);
    ye = xnod(LaG(e,:), Y);

    %% Construcción de la matriz C evaluando P en coordenadas naturales
    C = zeros(12,12);
    natural_coords = [-1 -1;
                       1 -1;
                       1 1;
                      -1 1]; % Coordenadas naturales de los nodos
    %for i = 1:4
        %C(3*(i-1)+1:3*i, :) = P(natural_coords(i,1), natural_coords(i,2));
    %    C(3*(i-1)+1:3*i, :) = P(xe(i,1), ye(i,1));
    %end
    for i = 1:4
        xi = natural_coords(i,1); yi =  natural_coords(i,2);
        %xi = xe(i); yi = ye(i);
        C(3*(i-1)+1:3*i, :) = [...
                                1,-xi, -yi,  -xi^2/2, -(xi^2*yi)/2 + yi^3/12, -yi^2/2,-(xi*yi^2)/2 + xi^3/12,-xi*yi/2, xi/2, xi*yi/2,yi/2, xi*yi/2;
                                0,1  , 0  ,       xi,                  xi*yi,       0,       xi^2/4 + yi^2/2,    yi/2,  1/2,    yi/2,   0, -yi/2;
                                0,0  , 1  ,        0,        yi^2/4 + xi^2/2,      yi,                 xi*yi,    xi/2,    0,   -xi/2, 1/2, xi/2];
    end
    invC = inv(C);
    K0 = zeros(12);
    M0 = zeros(12);
    det_Je = zeros(n_gl,n_gl); % almacenara los Jacobianos
    for pp = 1:n_gl
        for qq = 1:n_gl
            xi_gl  = x_gl(pp);            
            eta_gl = x_gl(qq);
            NN       = Nforma (xi_gl, eta_gl);
            %xm = sum(NN .* xe);
            %ym = sum(NN .* ye);
            ddN_dxi  = dN_dxi (xi_gl, eta_gl);       
            ddN_deta = dN_deta(xi_gl, eta_gl);  
            dx_dxi  = sum(ddN_dxi .*xe);   dy_dxi  = sum(ddN_dxi .*ye);
            dx_deta = sum(ddN_deta.*xe);   dy_deta = sum(ddN_deta.*ye);

            % Evaluación del Jacobiano
            Je = [ dx_dxi    dy_dxi
                   dx_deta   dy_deta ];
            det_Je(pp,qq) = det(Je);

            % Evaluar P y Q en punto de Gauss (coordenadas naturales)
            P_val = P(xi_gl, eta_gl);
            Q_val = Q(xi_gl, eta_gl);
            %B = Q_val * invC;
            %Bm = P_val * invC;
            %B = Q_val * invC; % Matriz deformación-desplazamiento
            %N = P_val * invC; % Funciones de forma de desplazamient
            K0 = K0 +  Q_val' * D * Q_val  * det_Je(pp,qq) * w_gl(pp) * w_gl(qq);
            M0 = M0 +  P_val' * T * P_val  * det_Je(pp,qq) * w_gl(pp) * w_gl(qq);
            %K0 = K0 +  invC'*Q_val' * D * Q_val *invC * det_Je(pp,qq) * w_gl(pp) * w_gl(qq);
            %M0 = M0 +  invC'*P_val' * T * P_val *invC * det_Je(pp,qq) * w_gl(pp) * w_gl(qq);
        end
    end
    %% se verifica que todos los determinantes sean positivos
    if any(det_Je(:) <= 0)
        error('Existen elementos con det(Je(xi,eta)) <= 0 %d.\n', e);
    end
    idx = [gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:) gdl(LaG(e,4),:)];
    K(idx,idx) = K(idx,idx) +  invC'*K0*invC ;
    M(idx,idx) = M(idx,idx) +  invC'*M0*invC ;
    %K(idx,idx) = K(idx,idx) +  K0 ;
    %M(idx,idx) = M(idx,idx) +  M0 ;
end
