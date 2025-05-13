clc
clear all
%%
% Numeracion local del EF serendipito rectangular de 8 nodos:

%  Reference Element T3:
%    |
%    1  3
%    |  |\
%    |  | \
%  eta  6  5
%    |  |   \
%    |  |    \
%    0  1--4--2
%    |
%    +--0--xi--1-->
% Coordenadas de los nodos

        %  xi   eta     
xnod =   [1    0    
          0    1    
          0    0    
          1/2  1/2  
          0    1/2  
          1/2  0    ];

nno = size(xnod, 1);
plot(xnod(:,1),xnod(:,2))
%% Se define la cuadratura de Gauss Legendre a utilizar
n=3;
xw=TriGaussPoints(n);
x_gl = xw(:,1);
e_gl = xw(:,2);
w_gl =  xw(:,3);
n_gl = size(x_gl,1);  %# Número de puntos de Gauss.

%[x_gl, w_gl]  = gausslegendre_quad(n_gl);
x_gl = sym(x_gl);
e_gl = sym(e_gl);
%% numero de terminos del polinomio interpolador
nterm = 4; % 1   xi_gl   eta_gl   xi_gl*eta_gl

%%  Se define la matriz A1
A1 = sym(zeros(nterm));
i = 0;
for p = 1:n_gl
      i = i+1;
      xi_gl  = x_gl(p);
      eta_gl = e_gl(p);
      A1(i,:) = [ 1....
            xi_gl, eta_gl,...
            xi_gl.* eta_gl ];

end

%% Se define la matriz A2
A2 = sym(zeros(nno, nterm));
for i = 1:nno
      xi_gl  = xnod(i,1);
      eta_gl = xnod(i,2);      
      A2(i,:) = [  1....
                xi_gl, eta_gl,...
                xi_gl.* eta_gl ];
end
        mensaje = 'EF lagrangiano rectangular de 16 nodos';   
%% Se reporta la matriz
fprintf('La matriz de interpolacion del %s es:\n', mensaje)
A = simplify(A2/A1); %= A2*inv(A1)
disp(A) 
disp('o en forma numerica:')
format long
disp(double(A))
