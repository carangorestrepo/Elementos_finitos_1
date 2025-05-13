%% Funciones de forma del elemento rectangular serendipito de 4 y 8 nodos 
clear, clc, close all

nno = 25; %  escoja entre {4, 8}.

X = 1; Y = 2;

%% Coordenadas de los nodos y numeracion local

% Numeracion local:
%        ^ eta
%        |
%        |
% 13---12--11--10---9
%  |   |   |   |    |
% 14--23--22--21----8
%  |   |   |   |    |
% 15--24--25--20----7
%  |   |   |   |    |  ----> xi
% 16--17--18--19----6
%  |   |   |   |    |
%  1---2---3---4----5

nod = [ ...
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
  
%% Se calculan las funciones de forma bidimensionales
xxi  = nod(:, X); 
eeta = nod(:, Y);

%A = [        ones(4,1),...
%    xxi.^2.*eeta, xxi.*eeta.^2,...
%           xxi.*eeta ]; 
A = [ ones(nno,1),...
         xxi, eeta,...
      xxi.^2, xxi.*eeta, eeta.^2,...
    xxi.^3,xxi.^2.*eeta,xxi.*eeta.^2,eeta.^3,...
 xxi.^4,xxi.^3.*eeta,xxi.^2.*eeta.^2,xxi.*eeta.^3,eeta.^4,...
 xxi.^4.*eeta,xxi.^3.*eeta.^2,xxi.^2.*eeta.^3,xxi.*eeta.^4,...
  xxi.^4.*eeta.^2,xxi.^3.*eeta.^3,xxi.^2.*eeta.^4,...
       xxi.^4.*eeta.^3,xxi.^3.*eeta.^4,...
        xxi.^4.*eeta.^4
          ];
N = cell(nno,1);
syms xi eta
for i = 1:nno
    % se arma el sistema de ecuaciones
    b = zeros(nno,1);   b(i) = 1;
    coef_alpha = A\b;

     N{i} = simplify([ 1,...
         xi, eta,...
      xi^2, xi*eta, eta^2,...
    xi^3,xi^2*eta,xi*eta^2,eta^3,...
 xi^4,xi^3*eta,xi^2*eta^2,xi*eta^3,eta^4,...
 xi^4*eta,xi^3*eta^2,xi^2*eta^3,xi*eta^4,...
  xi^4*eta^2,xi^3*eta^3,xi^2*eta^4,...
       xi^4*eta^3,xi^3*eta^4,...
        xi^4*eta^4]*coef_alpha);% 25 nodos
%    N{i} = simplify([ 1 eta*xi^2 xi*eta^2 xi*eta ]*coef_alpha);
end

%% Se imprimen las funciones de forma
fprintf('Funciones de forma serendipitas del elemento rectangular de %d nodos:\n', nno)
for i = 1:nno
   fprintf('N%d = ', i); disp(N{i})
end

%% Se calculan las derivadas de las funciones de forma con respecto a xi y
%% con respecto a eta y se imprimen (para referencias posteriores):
fprintf('\nDerivadas con respecto a xi:\n')
for i = 1:nno
   fprintf('dN%d_dxi = ',  i); disp(simplify(diff(N{i}, xi)))
end

fprintf('\nDerivadas con respecto a eta:\n')
for i = 1:nno
   fprintf('dN%d_deta = ', i); disp(simplify(diff(N{i}, eta)))
end

%% Se verifica la condicion de cuerpo rigido: sum(N) == 1
suma = 0;
for i = 1:nno
   suma = suma + N{i};
end
fprintf('\nSe verifica la condicion de cuerpo rigido: sum(N) == ');
disp(simplify(suma));

 
%% Se grafican las funciones de forma
XXI  = linspace(-1, 1, 50);
EETA = linspace(-1, 1, 50);
[XI,ETA] = meshgrid(XXI,EETA);

% calculo las esferitas
[xsp,ysp,zsp] = sphere;
xsp = 0.025*xsp;
ysp = 0.025*ysp;
zsp = 0.025*zsp;

for i = 1:nno
   figure                 % creo un lienzo
   grid on                % creo la rejilla
   hold on                % para que no se sobreescriban los graficos
   xlabel('\xi', 'FontSize',16)   % titulo eje X
   ylabel('\eta', 'FontSize',16)  % titulo eje Y
   title(sprintf('N_{%d}(\\xi,\\eta)',i), 'FontSize',20)

   % con este comando convierto la funcion de forma de tipo simbolico a
   % tipo funcion
   NN = matlabFunction(N{i}, 'Vars', [xi, eta]);   
   surf(XI, ETA, NN(XI,ETA))     % superficie
   
   % se grafican las esferitas en cada nodo
   for j=1:nno
      surf(xsp+nod(j,X), ysp+nod(j,Y), zsp+(i==j), 'facecolor', 'k')
   end 
   
   axis tight             % ejes apretados
   daspect([1 1 1])       % similar a axis equal pero en 3D
   view(3)                % vista tridimensional
end

%% Bye, bye!
return;