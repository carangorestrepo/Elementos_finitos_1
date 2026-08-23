function postproceso_zapata(xnod,LaG,in_col,geo,suelo,resultado,opciones)
%POSTPROCESO_ZAPATA
% Graficas compatibles con versiones antiguas de MATLAB.
% Se usa subplot en lugar de tiledlayout.

fprintf('\nDESPLAZAMIENTOS Y PRESIONES\n');
fprintf('w minimo = %12.6e m\n',min(resultado.w));
fprintf('w maximo = %12.6e m\n',max(resultado.w));
fprintf('p nodal minima = %12.6f kPa\n',min(resultado.presion_nodal));
fprintf('p nodal maxima = %12.6f kPa\n',max(resultado.presion_nodal));

if any(resultado.presion_nodal < -1e-8)
    warning(['Se obtiene traccion en algunos resortes Winkler. ' ...
             'El modelo actual es lineal y permite traccion.']);
end

if ~isfield(opciones,'graficar') || ~opciones.graficar
    return
end

escala=opciones.escala_deformada;
nno=size(xnod,1);

%% Malla
figure
hold on
axis equal
grid on

patch('Faces',LaG(~in_col,:),'Vertices',xnod, ...
      'FaceColor',[0.82 1.00 0.82],'EdgeColor',[0.4 0.4 0.4]);

patch('Faces',LaG(in_col,:),'Vertices',xnod, ...
      'FaceColor',[1.00 0.82 0.82],'EdgeColor',[0.4 0.4 0.4]);

plot([0 geo.Lx geo.Lx 0 0],[0 0 geo.Ly geo.Ly 0], ...
    'k-','LineWidth',2);

plot([geo.cx0 geo.cx0+geo.lcx geo.cx0+geo.lcx geo.cx0 geo.cx0], ...
     [geo.cy0 geo.cy0 geo.cy0+geo.lcy geo.cy0+geo.lcy geo.cy0], ...
     'b--','LineWidth',2);

xlabel('x [m]')
ylabel('y [m]')
title('Malla de la zapata')

%% Deformada
figure
trisurf(LaG,xnod(:,1),xnod(:,2),escala*resultado.w,resultado.w, ...
    'EdgeColor',[0.35 0.35 0.35]);
axis equal
grid on
view(3)
colorbar
xlabel('x [m]')
ylabel('y [m]')
zlabel('w escalado')
title(sprintf('Deformada transversal, escala = %g',escala))

%% Presion Winkler
figure
trisurf(LaG,xnod(:,1),xnod(:,2),zeros(nno,1), ...
    resultado.presion_nodal,'EdgeColor','none');
axis equal
view(2)
colorbar
xlabel('x [m]')
ylabel('y [m]')
title(sprintf('Presion Winkler nodal aproximada, k_s = %g kN/m^3', ...
    suelo.ks))

%% Momentos
figure

subplot(1,3,1)
trisurf(LaG,xnod(:,1),xnod(:,2),zeros(nno,1), ...
    resultado.Mx,'EdgeColor','none');
axis equal
view(2)
colorbar
xlabel('x')
ylabel('y')
title('M_x [kN m/m]')

subplot(1,3,2)
trisurf(LaG,xnod(:,1),xnod(:,2),zeros(nno,1), ...
    resultado.My,'EdgeColor','none');
axis equal
view(2)
colorbar
xlabel('x')
ylabel('y')
title('M_y [kN m/m]')

subplot(1,3,3)
trisurf(LaG,xnod(:,1),xnod(:,2),zeros(nno,1), ...
    resultado.Mxy,'EdgeColor','none');
axis equal
view(2)
colorbar
xlabel('x')
ylabel('y')
title('M_{xy} [kN m/m]')

%% Cortantes
figure

subplot(1,2,1)
trisurf(LaG,xnod(:,1),xnod(:,2),zeros(nno,1), ...
    resultado.Qx,'EdgeColor','none');
axis equal
view(2)
colorbar
xlabel('x')
ylabel('y')
title('Q_x [kN/m]')

subplot(1,2,2)
trisurf(LaG,xnod(:,1),xnod(:,2),zeros(nno,1), ...
    resultado.Qy,'EdgeColor','none');
axis equal
view(2)
colorbar
xlabel('x')
ylabel('y')
title('Q_y [kN/m]')

%% Patron de K
figure
spy(resultado.K)
title('Patron de la matriz global K')

end
