clc
clear all
%%
% Numeracion local del EF serendipito rectangular de 25 nodos:
%        ^ eta
%        |
%        |
% 13---12--11--10---9
%  |   |   |   |
% 14--23--22--21----8
%  |   |   |   |
% 15--24--25--20----7
%  |   |   |   |      ----> xi
% 16--17--18--19----6
%  |   |   |   |
%  1---2---3---4----5

% Coordenadas de los nodos
xnod = [ ...
%  xi   eta     % nodo   
   -1    -1     %  1
   -1/2  -1     %  2
    0    -1     %  3
    1/2  -1     %  4
    1    -1     %  5
    1   -1/2    %  6
    1    0      %  7
    1    1/2    %  8
    1    1      %  9
    1/2  1      % 10
    0    1      % 11
   -1/2  1      % 12
   -1    1      % 13
   -1    1/2    % 14
   -1    0      % 15
   -1   -1/2    % 16
   -1/2 -1/2    % 17
    0   -1/2    % 18
    1/2 -1/2    % 19
    1/2  0      % 20
    1/2  1/2    % 21
    0    1/2    % 22
   -1/2  1/2    % 23
   -1/2  0      % 24
    0    0];    % 25


nno = size(xnod, 1);
plot(xnod(:,1),xnod(:,2))
%% Se define la cuadratura de Gauss Legendre a utilizar
n_gl = 5;
[x_gl, w_gl]  = gausslegendre_quad(n_gl);
x_gl = sym(x_gl);
%% numero de terminos del polinomio interpolador
nterm = 25; % 1   xi_gl   eta_gl   xi_gl*eta_gl
%%  Se define la matriz A1
A1 = sym(zeros(n_gl*n_gl,nterm));
i = 0;
for p = 1:n_gl
   for q = 1:n_gl
      i = i+1;
      xi_gl  = x_gl(p);
      eta_gl = x_gl(q);
      %A1(i,:) = [ 1 xi_gl eta_gl xi_gl*eta_gl ];%% 4 terminios
      %A1(i,:) = [ 1 xi_gl eta_gl xi_gl^2 xi_gl*eta_gl eta_gl^2 xi_gl^2*eta_gl xi_gl*eta_gl^2 xi_gl^2*eta_gl^2];%% 9 terminios
      
     %A1(i,:) = [ 1 xi_gl eta_gl xi_gl^2 xi_gl*eta_gl eta_gl^2];%% 6 terminios 
     A1(i,:) = [       1,...
                 xi_gl, eta_gl,...
          xi_gl^2, xi_gl*eta_gl, eta_gl^2,...
        xi_gl^3, xi_gl^2*eta_gl, xi_gl*eta_gl^2,eta_gl^3,...
    xi_gl^4,xi_gl^3*eta_gl, xi_gl^2*eta_gl^2, xi_gl*eta_gl^3,eta_gl^4,...
    xi_gl^4*eta_gl, xi_gl^3*eta_gl^2, xi_gl^2*eta_gl^3, xi_gl*eta_gl^4,...
        xi_gl^4*eta_gl^2, xi_gl^3*eta_gl^3 xi_gl^2*eta_gl^4,...
                 xi_gl^4*eta_gl^3, xi_gl^3*eta_gl^4,...
                    xi_gl^4*eta_gl^4];%% 25 terminios
   end
end

%% Se define la matriz A2
A2 = sym(zeros(nno, nterm));
for i = 1:nno
      xi_gl  = xnod(i,1);
      eta_gl = xnod(i,2);      
      %A2(i,:) = [ 1 xi_gl eta_gl xi_gl*eta_gl ];
      %A2(i,:) =[ 1 xi_gl eta_gl xi_gl^2 xi_gl*eta_gl eta_gl^2];%% 6 terminios 
      %A2(i,:) =[ 1 xi_gl eta_gl xi_gl^2 xi_gl*eta_gl eta_gl^2 xi_gl^2*eta_gl xi_gl*eta_gl^2 xi_gl^2*eta_gl^2];
       A2(i,:) = [ 1,...
             xi_gl, eta_gl,...
             xi_gl^2, xi_gl*eta_gl, eta_gl^2,...
             xi_gl^3,  xi_gl^2*eta_gl, xi_gl*eta_gl^2,eta_gl^3,...
             xi_gl^4,  xi_gl^3*eta_gl, xi_gl^2*eta_gl^2, xi_gl*eta_gl^3,eta_gl^4,...
             xi_gl^4*eta_gl, xi_gl^3*eta_gl^2, xi_gl^2*eta_gl^3, xi_gl*eta_gl^4,...
             xi_gl^4*eta_gl^2, xi_gl^3*eta_gl^3, xi_gl^2*eta_gl^4,...
             xi_gl^4*eta_gl^3, xi_gl^3*eta_gl^4,...
             xi_gl^4*eta_gl^4];%% 25 terminios
      
end
        mensaje = 'EF lagrangiano rectangular de 25 nodos';   
%% Se reporta la matriz
fprintf('La matriz de interpolacion del %s es:\n', mensaje)
A = simplify(A2/A1); %= A2*inv(A1)
disp(A) 
disp('o en forma numerica:')
format long
disp(double(A))
