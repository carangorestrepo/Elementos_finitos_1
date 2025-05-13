clear, clc, close all
syms xi eta
nno=6;
%% Se grafican las funciones de forma
XXI  = linspace(-1, 1, 50);
EETA = linspace(-1, 1, 50);
[XI,ETA] = meshgrid(XXI,EETA);

Na =  [((eta - 1)*(xi - 1))/4;   % N1
                   -((eta - 1)*(xi + 1))/4;    % N2
                    ((eta + 1)*(xi + 1))/4;    % N3
                   -((eta + 1)*(xi - 1))/4;    % N4
                                  1 - xi^2;    % N5
                                 1 - eta^2 ];  % N6
N = cell(nno,1);             
for i = 1:nno
   N{i} = simplify(Na(i,1));
end

%% Se verifica la condicion de cuerpo rigido: sum(N) == 1
suma = 0;
for i = 1:nno
   suma = suma + N{i};
end
fprintf('\nSe verifica la condicion de cuerpo rigido: sum(N) == ');
disp(simplify(suma));

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
  
   axis tight             % ejes apretados
   daspect([1 1 1])       % similar a axis equal pero en 3D
   view(3)                % vista tridimensional
end