%% Funciones de forma del elemento rectangular lagrangiano de 9 nodos 
clear, clc

X = 1; Y = 2;
syms C1 C2 C3 C4  xi L eta EI Ac L x
L=4;
A=0.4^2;
I=0.4^4/12;
E=4700*sqrt(28)*1000;
G=0.4*E;
EI=E*I;
Ac=A*5/6*G;
q=0;
V=int(q,x)+C1;
M=int(V,x)+C2;
t=int(M/EI,x)+C3;
v=int(t-V/Ac,x)+C4;

[c1,c2,c3,c4]=solve(subs(v,x,0)==0,...
                    subs(v,x,L)==0,...
                    subs(t,x,0)==1,...
                    subs(t,x,L)==0,...
                    [C1,C2,C3,C4]); 
                

v1=subs(t,{C1,C2,C3,C4},{c1,c2,c3,c4});
v2=subs(v,{C1,C2,C3,C4,x},{c1,c2,c3,c4,eta});  

x1=linspace(0,L,200);
y=subs(v1,x,x1);
plot(x1,y)
 %t=1 xi=-1               
t1=(3*xi^2)/16 - xi/2 + 5/16;
 %t=1 eta=-1  
t2=(3*eta^2)/16 - eta/2 + 5/16;

%t=1 xi=1 
t3=(3*xi^2)/16 + xi/2 + 5/16;
%t=1 xi=1 
t4=(3*eta^2)/16 + eta/2 + 5/16;


syms xi eta
L1_xi = cell(2,1); % contenedor para las funciones de forma (en dir XI)
%v=1 xi=-1
L1_xi{1}=xi^3/16 - (9*xi)/16 + 1/2;
%v=1 xi=1
L1_xi{2}=- xi^3/16 + (9*xi)/16 + 1/2;

L1_eta = cell(2,1); % contenedor para las funciones de forma (en dir XI)
%v=1 xi=-1
L1_eta{1}=eta^3/16 - (9*eta)/16 + 1/2;
%v=1 xi=1
L1_eta{2}=- eta^3/16 + (9*eta)/16 + 1/2;

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
   for j=1:nno
      surf(xsp+nod(j,X), ysp+nod(j,Y), zsp+(i==j), 'facecolor', 'k')
   end 
   
   axis tight             % ejes apretados
   daspect([1 1 1])       % similar a axis equal pero en 3D
   view(3)                % vista tridimensional
end

