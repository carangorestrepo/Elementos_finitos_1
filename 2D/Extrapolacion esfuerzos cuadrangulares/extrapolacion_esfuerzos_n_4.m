clc
clear all
%%
% Numeracion local del EF serendipito rectangular de 25 nodos:
% Numeracion local:
%        ^ eta
%        |
%        |
% 10---9---8---7
%  |   |   |   |
% 11--16--15---6
%  |   |   |   |----> xi
% 12--13--14---5
%  |   |   |   |
%  1---2---3---4

% Coordenadas de los nodos
xnod = [ ...
%  xi   eta     % nodo   
   -1   -1      %  1
   -1/3 -1      %  2
    1/3 -1      %  3
    1   -1      %  4
    1   -1/3    %  5
    1    1/3    %  6
    1    1      %  7
    1/3  1      %  8
   -1/3  1      %  9
   -1    1      % 10
   -1    1/3    % 11
   -1   -1/3    % 12
   -1/3 -1/3    % 13
    1/3 -1/3    % 14
    1/3  1/3    % 15
   -1/3  1/3 ]; % 16


nno = size(xnod, 1);
plot(xnod(:,1),xnod(:,2))
%% Se define la cuadratura de Gauss Legendre a utilizar
n_gl = 4;
[x_gl, w_gl]  = gausslegendre_quad(n_gl);
x_gl = sym(x_gl);

%% numero de terminos del polinomio interpolador
nterm = 16; % 1   xi_gl   eta_gl   xi_gl*eta_gl

%%  Se define la matriz A1
A1 = sym(zeros(nterm));

i = 0;
for p = 1:n_gl
   for q = 1:n_gl
      i = i+1;
      xi_gl  = x_gl(p);
      eta_gl = x_gl(q);
         
      A1(i,:)= [ 1,...
        xi_gl, eta_gl,...
   xi_gl.^2, xi_gl.*eta_gl, eta_gl.^2....
 xi_gl.^3,xi_gl.^2.*eta_gl,xi_gl.*eta_gl.^2,eta_gl.^3,...
  xi_gl.^3.*eta_gl,xi_gl.^2.*eta_gl.^2,xi_gl.*eta_gl.^3,...
     xi_gl.^3.*eta_gl.^2,xi_gl.^2.*eta_gl.^3,...
          xi_gl.^3.*eta_gl.^3];
      
   end
end

%% Se define la matriz A2
A2 = sym(zeros(nno, nterm));
for i = 1:nno
      xi_gl  = xnod(i,1);
      eta_gl = xnod(i,2);      
               A2(i,:) = [1,...
                 xi_gl, eta_gl,...
            xi_gl.^2, xi_gl.*eta_gl, eta_gl.^2....
        xi_gl.^3,xi_gl.^2.*eta_gl,xi_gl.*eta_gl.^2,eta_gl.^3,...
    xi_gl.^3.*eta_gl,xi_gl.^2.*eta_gl.^2,xi_gl.*eta_gl.^3,...
     xi_gl.^3.*eta_gl.^2,xi_gl.^2.*eta_gl.^3,...
          xi_gl.^3.*eta_gl.^3];
end
        mensaje = 'EF lagrangiano rectangular de 16 nodos';   
%% Se reporta la matriz
fprintf('La matriz de interpolacion del %s es:\n', mensaje)
A = simplify(A2/A1); %= A2*inv(A1)
disp(A) 
disp('o en forma numerica:')
format long
disp(double(A))
