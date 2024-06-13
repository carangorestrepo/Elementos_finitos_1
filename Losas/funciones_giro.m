X = 1; Y = 2;
syms xi eta
A=1*0.2;
I=1*0.2^3/12;
E=1;%4700*sqrt(28)*1000;
G=1;%0.4*E;
EI=E*I;
Ac=A*5/6*G;



L1_xi = cell(2,1); % contenedor para las funciones de forma (en dir XI)
%v=1 xi=-1
L1_xi{1}=(3*Ac*x^2)/(Ac*L^2 + 12*EI) - (4*(Ac*L^2 + 3*EI)*x)/(L*(Ac*L^2 + 12*EI)) + 1;
%v=1 xi=1
L1_xi{2}=(3*Ac*x^2)/(Ac*L^2 + 12*EI) + (2*(6*EI - Ac*L^2)*x)/(L*(Ac*L^2 + 12*EI));

L1_eta = cell(2,1); % contenedor para las funciones de forma (en dir XI)
%v=1 xi=-1
L1_eta{1}=(Ac*eta^3)/(4*Ac + 12*EI) - eta^2/4 + (- (3*EI)/(2*Ac + 6*EI) - (Ac - 6*EI)/(4*Ac + 12*EI))*eta + 1/4;
%v=1 xi=1
L1_eta{2}=(Ac*eta^3)/(4*Ac + 12*EI) + eta^2/4 + (- (3*EI)/(2*Ac + 6*EI) - (Ac - 6*EI)/(4*Ac + 12*EI))*eta - 1/4;

nod = [ ...
%  xi   eta     % nodo   
   -1   -1      %  1
    1   -1      %  2
    1    1      %  3
   -1    1    ];%  9
nno = size(nod, 1);
pos = nod;
pos(nod==-1) = 1;
pos(nod== 1) = 2;
%% Se calculan las funciones de forma bidimensionales
N = cell(nno,1);
for i = 1:nno
   N{i} = simplify(L1_xi{pos(i,1)}*L1_eta{pos(i,2)});
end

%% Se imprimen las funciones de forma
fprintf('Funciones de forma lagrangianas del elemento rectangular de 9 nodos:\n')
for i = 1:nno
   fprintf('N%d = ', i); disp(N{i})
end
%% Se verifica la condicion de cuerpo rigido: sum(N) == 1
suma = 0;
for i = 1:nno
   suma = suma + N{i};
end
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
   %for j=1:nno
   %   surf(xsp+nod(j,X), ysp+nod(j,Y), zsp+(i==j), 'facecolor', 'k')
   %end 
   
   axis tight             % ejes apretados
   daspect([1 1 0.3])       % similar a axis equal pero en 3D
   view(3)                % vista tridimensional
end

