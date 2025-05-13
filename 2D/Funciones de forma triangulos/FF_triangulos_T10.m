clc
clear
%phi=a0+a1*x+a2*y+a3*x*y+a4*x2+a5*y5;
nno = 10; %  escoja entre {4, 8}.

X = 1; Y = 2;

     %  xi   eta     % nodo  
coord =  [1    0    0
          0    1    0
          0    0    1
          2/3  1/3  0
          1/3  2/3  0
          0    2/3  1/3
          0    1/3  2/3
          1/3  0    2/3
          2/3  0    1/3             
          1/3  1/3  1/3];
coord = round(coord*100)/100;
%% Se calculan las funciones de forma bidimensionales
xxi  = coord(:, X); 
eeta = coord(:, Y);
A = [ ones(nno,1) xxi eeta xxi.^2 xxi.*eeta eeta.^2 xxi.^3 xxi.^2.*eeta xxi.*eeta.^2 eeta.^3];
N = cell(nno,1);
syms xi eta
for i = 1:nno
   % se arma el sistema de ecuaciones
   b = zeros(nno,1);   b(i) = 1;
   coef_alpha = A\b;
   N{i} = simplify([ 1 xi eta xi^2 xi*eta eta^2 xi^3 xi.^2*eta xi.*eta.^2 eta.^3]*coef_alpha);
end
%% Se imprimen las funciones de forma
fprintf('Funciones de forma serendipitas del elemento rectangular de %d nodos:\n', nno)
for i = 1:nno
   fprintf('N%d = ', i); disp(N{i})
end

%% Se grafican las funciones de forma.
LL2 = 0:0.05:1;
LL3 = 0:0.05:1;
[L2,L3] = meshgrid(LL2,LL3);
L1 = 1 - L2 - L3;
L1 = round(100*L1)/100;   % redondear bien!
L1(L1<0) = NaN;

% coordenadas del triangulo
x = [0 1.0 0.5];
y = [0 0 sqrt(0.75)];
X = L1*x(1) + L2*x(2) + L3*x(3);
Y = L1*y(1) + L2*y(2) + L3*y(3);
X = X(:); 
Y = Y(:); 
isnanX = isnan(X);
X(isnanX) = [];
Y(isnanX) = [];

TRI = delaunay(X,Y);
numtriang = size(TRI, 1);   % numero de triangulos
% Se eliminan de la lista de triangulos aquellos cuya area sea negativa o 
% igual a cero
Area = zeros(numtriang,1);
for i = 1:numtriang
   Area(i) = 0.5*det([ 1 X(TRI(i,1)) Y(TRI(i,1))      %Area del EF e
                       1 X(TRI(i,2)) Y(TRI(i,2))
                       1 X(TRI(i,3)) Y(TRI(i,3))]);
end
TRI(Area < 1e-3,:) = [];

[xsp,ysp,zsp] = sphere;
xsp = 0.025*xsp;
ysp = 0.025*ysp;
zsp = 0.025*zsp;

for i=1:nno
   % creo una funcion que pueda evaluar numericamente
   NN = matlabFunction(N{i}, 'Vars', [xi, eta]);

   Z = NN(L1,L2);
   Z = Z(:);
   Z = Z(~isnanX);
   figure
   trimesh(TRI,X,Y,Z,'LineWidth',2);
   hold on
   trisurf(TRI,X,Y,Z);   
   alpha 0.3
   shading interp
   colormap winter
%   axis([-0.1 1.1 -0.1 0.96 -0.4 1.1])
   axis tight
   for j=1:nno
      cxsp = coord(j,:)*x';
      cysp = coord(j,:)*y';
      surf(xsp+cxsp, ysp+cysp, zsp+0);  % sphere centered at (x(1),y(1),0)
   end
   daspect([1 1 1]);
   title(sprintf('N_{%d} = %s',i,char(N{i})),'FontSize',20);
%   print('-dpdf',sprintf('%d.pdf',i));
end

%% Se verifica la condicion de cuerpo rigido: sum(N) == 1
syms a b
suma = subs(sum([N{:}]), {'L1', 'L2', 'L3'}, {1-a-b, a, b});
fprintf('\nSe verifica la condicion de cuerpo rigido: sum(N) == ');
disp(simplify(suma));

%% Bye, bye!
return;

%% Calcular correctamente los polinomios de las funciones de forma 1D
function N = calc_N(xp, yp, var)
    % se ve verifican los tamanios de los vectores xp y yp
    nx = length(xp);
    ny = length(yp);
    assert(nx == ny, 'Los vectores xp y yp deben tener el mismo tamanio');

    % se calculan los coeficientes de los polinomios
    c = polyfit(xp, yp, nx-1);
    
    % se eliminan los errores en la aproximacion numerica, haciendo los
    % coeficientes demasiado pequenios igual a cero
    c(abs(c) < 1e-10) = 0;
    
    % con los coeficientes corregidos se calculan las funciones de forma
    N = poly2sym(c, var);
end