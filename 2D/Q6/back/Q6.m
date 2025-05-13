

%% se leen algunas variables
%T        = readcell(archivo_xlsx, 'Sheet','varios','Range','B1:B9');
rhoe = 7850;        % densidad (kg/m^3)
g        = 9.81;    % [m/s²]   aceleración de la gravedad
U_LONG   = 'm';     % unidades de longitud
U_FUERZA = 'kN';    % unidades de fuerza
U_ESFUER = 'kN/m2'; % unidades de esfuerzo
ESC_UV   = 100;    % factor de escala para los desplazamientos
% 
%% constantes que ayudarán en la lectura del código

X = 1;
Y = 2;
NL1=1; NL2=2; NL3=3; NL4=4; NL5=5; NL6 =6;
r_ = [ 1, 2, 3, 4, 5, 6, 7, 8 ];   % GDL a retener en condensación nodal
e_ = [ 9, 10, 11, 12 ];            % GDL a eliminar en condensación nodal

 %% seleccione la malla a emplear:
%nombre_archivo = 'malla1'    # EJEMPLO CLASE
%nombre='nombre_archivo';
%nombre_archivo = 'malla1_no_estructurada';%    %EJEMPLO CLASE
%xnodxy = xlsread(nombre_archivo,1);% xnod;
%LaG_mat = xlsread(nombre_archivo,2);% LaG;
%restric = xlsread(nombre_archivo,3);% restric;
%carga_punt=xlsread(nombre_archivo,5);% carga_punt;
%prop_mat=xlsread(nombre_archivo,6);% prop_mat;
%carga_distr = xlsread(nombre_archivo,4);% carga_distr;
% %% posición de los nodos:
% xnod: fila=número del nodo, columna=coordenada X=0 o Y=1
%xnod=xnodxy(:,2:3);
%nno  = size(xnod,1);%   %# número de nodos (número de filas de la matriz xnod)

%% definición de los grados de libertad

ngdl = 2*nno;            % número de grados de libertad por nodo = [X, Y]
gdl  = [(1:2:(nno*2))',(1:2:(nno*2))'+1]; % nodos vs grados de libertad

%% definición de elementos finitos con respecto a nodos
% LaG: fila=número del elemento, columna=número del nodo local

%LaG = LaG_mat(:,2:5);
%nef = size(LaG,1);     % número de EFs (número de filas de la matriz LaG)

%% definición de los materiales
%mat  = LaG_mat(:,6);
%Ee   = prop_mat(:,2);     % [Pa]     módulo de elasticidad del sólido
%nue  = prop_mat(:,3);     % [-]      coeficiente de Poisson
%rhoe = prop_mat(:,4);     % [kg/m³]  densidad
%te   = prop_mat(:,5);     % [m]      espesor
%nmat =size(Ee,1);         % número de materiales

%% relación de cargas puntuales
%cp  = carga_punt(:,
%ncp = size(carga_punt,1);    % número de cargas puntuales
%f   = zeros(ngdl,1);           % vector de fuerzas nodales equivalentes global
%for i=1: ncp
%   f(carga_punt(i,1)*2-carga_punt(i,2)+1,1)=carga_punt(i,3);
%end

%% Se dibuja la malla de elementos finitos
cg = zeros(nef,2);%  % almacena el centro de gravedad de los EF
figure()
hold on
for e =1:nef
    % se dibujan las aristas
    nod_ef = LaG(e, [NL1, NL2, NL3, NL4, NL1]);
    plot(xnod(nod_ef, X), xnod(nod_ef, Y), 'b')
    % se calcula la posición del centro de gravedad
    cg(e,:) = [mean(xnod(nod_ef, X)),mean(xnod(nod_ef, Y))];
    % y se reporta el número del elemento actual
    text(cg(e,X), cg(e,Y),num2str(e), 'Color', 'b','HorizontalAlignment', 'center')
end
% en todos los nodos se dibuja un marcador y se reporta su numeraciÃ³n
plot(xnod(:,X), xnod(:,Y), 'ro');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
axis equal tight
title('Malla de elementos finitos');


%% Funciones de forma y sus derivadas del elemento rectangular de 4 nodos (N1 a
%   N4) y las funciones de forma de los modos incompatibles (N5 y N6)
Nforma = @(xi,eta)[ ((eta - 1)*(xi - 1))/4    % N1
                   -((eta - 1)*(xi + 1))/4    % N2
                    ((eta + 1)*(xi + 1))/4    % N3
                   -((eta + 1)*(xi - 1))/4 ]; % N4
%                           b   1 - xi**2,    % N5
%                              1 - eta**2  ]) % N6

% derivadas de las funciones de forma con respecto a xi
dN_dxi = @(xi,eta) [  eta/4 - 1/4    % dN1_dxi
                      1/4 - eta/4    %dN2_dxi
                      eta/4 + 1/4    %dN3_dxi
                    - eta/4 - 1/4    %dN4_dxi
                           - 2*xi    %dN5_dxi
                                0];  %dN6_dxi


%% Derivadas de N con respecto a eta
dN_deta =  @(xi,eta) [  xi/4 - 1/4   % dN1_deta
                       -xi/4 - 1/4   % dN2_deta
                        xi/4 + 1/4   % dN3_deta
                        1/4 - xi/4   % dN4_deta
                                 0   % dN5_deta
                          - 2*eta];  % dN6_deta

%% Cuadratura de Gauss-Legendre
%NOTA: se asumirá aquí el mismo orden de la cuadratura tanto en la dirección
%       de xi como en la dirección de eta
n_gl         = 2; % orden de la cuadratura de Gauss-Legendre
[x_gl, w_gl] = gausslegendre_quad(n_gl);

 %% Ensamblaje la matriz de rigidez global y el vector de fuerzas másicas
%    nodales equivalentes global

% se inicializan la matriz de rigidez global y los espacios en memoria que
%  almacenarán las matrices de forma y de deformación
K       = zeros(ngdl, ngdl);         % matriz de rigidez global
inv_Kee = cell(nef, 4, 4);
Ker     = cell(nef, 4, 8);

N = cell(nef,n_gl,n_gl,2,2*4); % matriz de forma en cada punto de GL
B = cell(nef,n_gl,n_gl,3,2*6); % matriz de deformaciones en cada punto de GL
idx = cell(nef,1);                 % indices asociados a los gdl del EF e

% matriz constitutiva del elemento para TENSION PLANA
%De =  cell(nmat,1);
%be = cell(nmat,1);

%for i = 1:nmat
    % matriz constitutiva para TENSION PLANA
De = [ Ee/(1-nue^2)     Ee*nue/(1-nue^2)  0
       Ee*nue/(1-nue^2) Ee/(1-nue^2)      0
       0                0                 Ee/(2*(1+nue)) ];
    % vector de fuerzas mÃ¡sicas
    %be{i} = [0; -rhoe(i)*g];
%end
% para cada elemento finito en la malla:
for e = 1:nef
   % se calculan con el siguiente ciclo las matrices de rigidez y el vector de
   % fuerzas nodales equivalentes del elemento usando las cuadraturas de GL
   Ke16 = zeros(12);
   fe = zeros(8,1);
   det_Je = zeros(n_gl,n_gl); % matriz para almacenar los jacobianos

   % se determinan las coordenadas de los nodos el EF e
   xe = xnod(LaG(e,:),X);
   ye = xnod(LaG(e,:),Y);

   for p = 1:n_gl
      for q = 1:n_gl
         xi_gl  = x_gl(p);
         eta_gl = x_gl(q);
         
         % Se evaluan las funciones de forma y sus derivadas 
         % en los puntos de integracion de Gauss-Legendre
         NNforma  = Nforma (xi_gl, eta_gl);
         ddN_dxi  = dN_dxi (xi_gl, eta_gl);%% integran 6 veses
         ddN_deta = dN_deta(xi_gl, eta_gl);%% integran 6 veses
         
         dx_dxi  = sum(ddN_dxi(1:4)  .* xe);   dy_dxi  = sum(ddN_dxi(1:4)  .* ye);
         dx_deta = sum(ddN_deta(1:4) .* xe);   dy_deta = sum(ddN_deta(1:4) .* ye);
         
         % Se ensambla la matriz Jacobiana del elemento
         Je = [ dx_dxi   dy_dxi
                dx_deta  dy_deta ];
            
         % Se calcula el determinante del Jacobiano
         det_Je(p,q) = det(Je);
         
         % las matrices de forma y de deformaciÃ³n se evalÃºan y se ensamblan
         % en el punto de Gauss         
         N{e,p,q} = zeros(2, 2*4);
         B{e,p,q} = zeros(3, 2*6);
         for i = 1:4
            % Se ensambla la matriz de funciones de forma N
            N{e,p,q}(:,[2*i-1 2*i]) = [ NNforma(i)  0         
                                        0           NNforma(i) ];
         end
         for i=1:6
            % Se ensambla la matriz de deformacion del elemento B
            dNi_dx = (+dy_deta*ddN_dxi(i) - dy_dxi*ddN_deta(i))/det_Je(p,q);
            dNi_dy = (-dx_deta*ddN_dxi(i) + dx_dxi*ddN_deta(i))/det_Je(p,q);
            B{e,p,q}(:,[2*i-1 2*i]) = [ dNi_dx       0        
                                             0  dNi_dy   
                                        dNi_dy  dNi_dx ];
         end
         % se ensamblan la matriz de rigidez del EF e y el vector de fuerzas
         % nodales equivalentes del EF e asociado a la fuerza mÃ¡sica         
         Ke16 = Ke16 + B{e,p,q}'*De*B{e,p,q} * det_Je(p,q)*te*w_gl(p)*w_gl(q);
         %fe =   fe + N{e,p,q}'*be{mat(e)}   * det_Je(p,q)*te(mat(e))*w_gl(p)*w_gl(q);
      end
   end
   
   % se determina si hay puntos con jacobiano negativo, en caso tal se termina
   % el programa y se reporta   
   if any(any(det_Je <= 0))
      error('Hay puntos con det_Je negativo en el EF %d.\n', e);
   end
   
    % se condensan los GDL jerárquicos u5, v5, u6, v6
    Krr = Ke16(r_,r_); 
    Ker{e}= Ke16(e_,r_);
    Kre = Ke16(r_,e_);
    inv_Kee{e} = (Ke16(e_,e_))^(-1);
    Ke = Krr - Kre * inv_Kee{e} * Ker{e};

   % y se ensambla la matriz de rigidez del elemento y el vector de fuerzas
   % nodales del elemento en sus correspondientes GDL 
   idx{e}=reshape([LaG(e, [NL1, NL2, NL3, NL4])*2-1;LaG(e, [NL1, NL2, NL3, NL4])*2],1,8);
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke;
   %f(idx{e},:)      = f(idx{e},:)      + fe;
end
%% Muestro la configuracion de la matriz K (K es rala)
figure
spy(K);
title('Los puntos representan los elementos diferentes de cero');


%% Cálculo de las cargas nodales equivalentes de las cargas distribuidas:
%cd   = carga_distr;
%nlcd = size(cd,1);     % número de lados con carga distribuida
%ft   = zeros(ngdl,1);  % fuerzas nodales equivalentes de cargas superficiales

%for i =1:nlcd
%   e     = cd(i,1);
%   lado  = cd(i,2);
%   carga = cd(i,3:6);
%   fte = t2ft_R4(xnod(LaG(e,:),:), lado, carga, te(mat(e)));
%   ft(idx{e},:) = ft(idx{e},:) + fte;
%end
%% se definen los apoyos y sus desplazamientos

idxNODO  = [ap';ap'];%restric(:,1);%% nodos 
sizeap=size(ap,2);
dir_desp =[ones(sizeap,1);ones(sizeap,1)*2] ;%% dirección

% grados de libertad del desplazamiento conocidos  
% desplazamientos conocidos en los apoyos
ac = [zeros(sizeap,1);zeros(sizeap,1)];%% desplazamiento

ngdl_res = size(ac,1); % numero de grados de libertad restringidos
c=zeros(ngdl_res,1);
for i = 1:ngdl_res
   %nodo        direccion
   c(i) = gdl(idxNODO(i), dir_desp(i));
end
%c = gdl(sub2ind([nno 2], idxNODO, dir_desp));
% grados de libertad del desplazamiento desconocidos  
d = setdiff(1:ngdl,c(:,1))';

%% extraigo las submatrices y especifico las cantidades conocidas
% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd |  % recuerde que qc=0 (siempre)
%|    | = |         ||    | - |    |
%| qc |   | Kdc Kdd || ad |   | fc |
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc - Kdc*ac);      % desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % fuerzas de equilibrio desconocidas

% armo los vectores de desplazamientos (a) y fuerzas (q)
a = zeros(ngdl,1);  a(c) = ac;  a(d) = ad; % desplazamientos
q = zeros(ngdl,1);  q(c) = qd;             % fuerzas nodales equivalentes

%% Dibujo la estructura original y la deformada
delta = reshape(a,2,nno)';
xdef = xnod + ESC_UV*delta; % posicion de la deformada
figure
hold on
for e = 1:nef
   line(xnod(LaG(e,[1:4 1]),X), xnod(LaG(e,[1:4 1]),Y), 'Color','b'); % original
   line(xdef(LaG(e,[1:4 1]),X), xdef(LaG(e,[1:4 1]),Y), 'Color','r'); % deformada
end
xlabel(['Eje X [' U_LONG ']']);
ylabel(['Eje Y [' U_LONG ']']);
axis equal tight;
legend('Posicion original','Posicion deformada','Location', 'SouthOutside');
title(sprintf('Deformada escalada %d veces', ESC_UV));

%% Se calcula para cada elemento las deformaciones y los esfuerzos
def = cell(nef,n_gl,n_gl);
esf = cell(nef,n_gl,n_gl);

for e = 1:nef
   % desplazamientos de los gdl del elemento e
   ar = a(idx{e});            % desplazamientos de los gdl del elemento e
   ae =  [ar;
          inv_Kee{e}  * Ker{e} * ar];
   for pp = 1:n_gl
      for qq = 1:n_gl
         def{e,pp,qq} = B{e,pp,qq}*ae;           % calculo las deformaciones
         esf{e,pp,qq} = De*def{e,pp,qq}; % calculo los esfuerzos
      end
   end
end

%% Se extrapolan los esfuerzos y las deformaciones a los nodos y se alisan
% adicionalmente se calcula el error en el alisado
[sx,  error_sx ] = extrapolar_esf_def(xnod, LaG, esf, 'sx');
[sy,  error_sy ] = extrapolar_esf_def(xnod, LaG, esf, 'sy');
[txy, error_txy] = extrapolar_esf_def(xnod, LaG, esf, 'txy');
[ex,  error_ex ] = extrapolar_esf_def(xnod, LaG, def, 'ex');
[ey,  error_ey ] = extrapolar_esf_def(xnod, LaG, def, 'ey');
[gxy, error_gxy] = extrapolar_esf_def(xnod, LaG, def, 'gxy');

%% en tension plana ...
sz   = zeros(nno,1);
txz  = zeros(nno,1);
tyz  = zeros(nno,1);

ez   = -(nue/Ee)*(sx+sy);              % deformaciones ez
tmax = sqrt(((sx-sy)/2).^2+txy.^2);  % esfuerzo cortante maximo
s1   = (sx+sy)/2 + tmax;             % esfuerzo normal maximo
s2   = (sx+sy)/2 - tmax;             % esfuerzo normal minimo
ang  = 0.5*atan2(2*txy, sx-sy);      % angulo de inclinacion de s1

s3   = zeros(nno,1);
sv   = sqrt(((s1-s2).^2 + (s2-s3).^2 + (s1-s3).^2)/2); % von Mises

sxx=reshape(sx,Ny,Nx);
syy=reshape(sy,Ny,Nx);
txyy=reshape(txy,Ny,Nx);

%{
figure
hold on
for e = 1:nef
    x=xnod(LaG(e,:),X); 
    y=xnod(LaG(e,:),Y);
    z= sx(LaG(e,:));
    trisurfc(x,y,z,3)
end
%% Se reportan los resultados en un archivo .xlsx
%tabla_aq = array2table([(1:nno)', reshape(a,2,nno)', reshape(q,2,nno)'], ...
%    'VariableNames', {'nodo', ['u_' U_LONG], 'v', ['qx_' U_FUERZA], 'qy'});

%tabla_def = array2table([(1:nno)', ex, ey, ez, gxy], ...
%    'VariableNames', {'nodo', 'ex', 'ey', 'ez', 'gxy_rad'});

%tabla_esf = array2table([(1:nno)', sx, sy, txy, s1, s2, ang, tmax, sv], ...
%    'VariableNames', {'nodo', ['sx_' U_ESFUER], 'sy', 'txy', 's1', 's2', 'ang_rad', 'tmax', 'sv'});

%nombre_archivo_results = ['resultados_' nombre_archivo{2} '.xlsx'];
%writetable(tabla_aq,  nombre_archivo_results, 'Sheet', 'aq')
%writetable(tabla_def, nombre_archivo_results, 'Sheet', 'deformaciones')
%writetable(tabla_esf, nombre_archivo_results, 'Sheet', 'esfuerzos')

%fprintf('Calculo finalizado. Resultados en "%s".\n', nombre_archivo_results);

%% se grafican las deformaciones

figure
subplot(4,1,1); plot_def_esf(xnod, LaG, ex,  '\epsilon_x')
subplot(4,1,2); plot_def_esf(xnod, LaG, ey,  '\epsilon_y')
subplot(4,1,3); plot_def_esf(xnod, LaG, ez,  '\epsilon_z')
subplot(4,1,4); plot_def_esf(xnod, LaG, gxy, '\gamma_{xy} [rad]')

%% se grafican los esfuerzos
figure
subplot(3,1,1); plot_def_esf(xnod, LaG, sx,  ['\sigma_x [' U_ESFUER ']'])
subplot(3,1,2); plot_def_esf(xnod, LaG, sy,  ['\sigma_y [' U_ESFUER ']'])
subplot(3,1,3); plot_def_esf(xnod, LaG, txy, ['\tau_{xy} [' U_ESFUER ']'])

%% se grafican los errores en los esfuerzos
%figure
%subplot(1,3,1); plot_def_esf(xnod, LaG, error_sx,  ['Error \sigma_x [' U_ESFUER ']'])
%subplot(1,3,2); plot_def_esf(xnod, LaG, error_sy,  ['Error \sigma_y [' U_ESFUER ']'])
%subplot(1,3,3); plot_def_esf(xnod, LaG, error_txy, ['Error \tau_{xy} [' U_ESFUER ']'])

%% se grafican los esfuerzos principales y el esfuerzo cortante maximo
figure
subplot(3,1,1); plot_def_esf(xnod, LaG, s1,   ['(\sigma_1)_{xy} [' U_ESFUER ']'], { ang })
subplot(3,1,2); plot_def_esf(xnod, LaG, s2,   ['(\sigma_2)_{xy} [' U_ESFUER ']'], { ang+pi/2 })
subplot(3,1,3); plot_def_esf(xnod, LaG, tmax, ['\tau_{max} [' U_ESFUER ']'],      { ang+pi/4, ang-pi/4 })

%% se grafican los esfuerzos de von Mises

figure
plot_def_esf(xnod, LaG, sv, ['Esfuerzos de von Mises [' U_ESFUER ']']);

%% se exportan los resultados a GiD/Paraview
% Pasando los esfuerzos ya promediados:
%export_to_GiD('c5_ejemplo_a',xnod,LaG,a,q,[sx sy sz txy txz tyz]);

% Pasando los puntos de Gauss [RECOMENDADO] !!!
% export_to_GiD('c5_ejemplo_b',xnod,LaG,a,q,esf);                    

%%
return; % bye, bye!

%% Lee del archivo "nombre_archivo" de EXCEL la hoja "hoja"
function H = leer_excel(archivo_xlsx, hoja)

    if verLessThan('matlab', '9.9') % R2019b or older
        H = readtable(archivo_xlsx, 'Sheet', hoja);
    else
        H = readtable(archivo_xlsx, 'Sheet', hoja, 'format', 'auto');
    end
end

%% Grafica los esfuerzos y las deformaciones
function plot_def_esf(xnod, LaG, variable, texto, angulos)
    X = 1; Y = 2;

    colorbar;
    
    nef = size(LaG, 1);    
    for e = 1:nef
       data.numEF = e; 
       data.LaG_e = LaG(e,:);
       fill(xnod(LaG(e,:),X), xnod(LaG(e,:),Y), variable(LaG(e,:)), ...
           'UserData', data), hold on;
    end
    axis equal tight
    colormap jet
    title(texto);

    esc = 0.5;
    if nargin == 5
        norma = 1; % = variable % si se quiere proporcional
        for i = 1:length(angulos)
            % se indica la flecha de la direccion principal
            quiver(xnod(:,X),xnod(:,Y),...             
                norma.*cos(angulos{i}), norma.*sin(angulos{i}),... 
                esc, ...                  % con una escala esc
                'k',...                   % de color negro
                'ShowArrowHead','off',... % una flecha sin cabeza
                'LineWidth',2, ...        % con un ancho de linea 2
                'Marker','.'), hold on;            % y en el punto (x,y) poner un punto '.'
            
            % la misma flecha girada 180 grados
            quiver(xnod(:,X),xnod(:,Y),...             
                norma.*cos(angulos{i}+pi), norma.*sin(angulos{i}+pi),... 
                esc,'k', 'ShowArrowHead','off', 'LineWidth',2, 'Marker','.'), hold on;                    
        end      
    end

    % para mostrar el tooltip
    dcm = datacursormode;
    dcm.Enable = 'on';
    dcm.SnapToDataVertex = 'on';
    dcm.UpdateFcn = @mostrar_info_nodo;    
end

function txt = mostrar_info_nodo(~,info)
    if strcmp(info.Target.Type, 'patch')
        x      = info.Position(1);
        y      = info.Position(2);
        numEF  = info.Target.UserData.numEF;
        xnod_e = info.Target.Vertices;
    
        % se busca el punto mÃ¡s cercano
        [~, idx_e] = min(hypot(xnod_e(:,1) - x, xnod_e(:,2) - y));
        numNOD  = info.Target.UserData.LaG_e(idx_e);
        val_esf = info.Target.CData(idx_e);

        txt = ['(x,y) = (' num2str(x) ',' num2str(y) '), ' ...
               'nodo = '     num2str(numNOD)          ', ' ...
               'EF = '       num2str(numEF)           ', ' ...
               'variable = ' num2str(val_esf)];
    else
        txt = 'Haga zoom y evite seleccionar el quiver()';
    end
end
%}
%% Extrapola/alisa esfuerzos y deformaciones de puntos de Gauss a los nodos
function [esf, error_esf] = extrapolar_esf_def(xnod, LaG, esfuerzo, tipo_esf)
    nno = size(xnod, 1);
    nef = size(LaG, 1);

    num_elem_ady = zeros(nno,1);  % numero de elementos adyacentes
    esf.sum      = zeros(nno,1);
    esf.max      =  -inf(nno,1);
    esf.min      =   inf(nno,1);

    % matriz de extrapolación
    A = [[  3^(1/2)/2 + 1,            -1/2,            -1/2,   1 - 3^(1/2)/2];
        [            -1/2,   1 - 3^(1/2)/2,   3^(1/2)/2 + 1,            -1/2];
        [   1 - 3^(1/2)/2,            -1/2,            -1/2,   3^(1/2)/2 + 1];
        [            -1/2,   3^(1/2)/2 + 1,   1 - 3^(1/2)/2,            -1/2]];

    switch tipo_esf
        case {'sx',  'ex'},  num_esf = 1;
        case {'sy',  'ey'},  num_esf = 2;
        case {'txy', 'gxy'}, num_esf = 3;
        otherwise,           error('Opcion no soportada');
    end
    
    % se hace la extrapolaciÃ³n de los esfuerzos en cada EF a partir de las 
    % lecturas en los puntos de Gauss
    for e = 1:nef
        esf_EF_e = A * [ esfuerzo{e,1,1}(num_esf)
                         esfuerzo{e,1,2}(num_esf)
                         esfuerzo{e,2,1}(num_esf)
                         esfuerzo{e,2,2}(num_esf) ];        
        
        esf.sum(LaG(e,:),:) = esf.sum(LaG(e,:),:)    + esf_EF_e;
        esf.max(LaG(e,:),:) = max(esf.max(LaG(e,:),:), esf_EF_e);
        esf.min(LaG(e,:),:) = min(esf.max(LaG(e,:),:), esf_EF_e);     
          
        % se lleva un conteo de los elementos adyacentes a un nodo
        num_elem_ady(LaG(e,:),:) = num_elem_ady(LaG(e,:),:) + 1;
    end

    %% alisado (promedio de los esfuerzos en los nodos)
    esf.prom = esf.sum./num_elem_ady;    
    
    %% variables a retornar
    error_esf = (esf.max - esf.min)./esf.prom; % error en el alisado
    error_esf = log10(abs(error_esf));
    error_esf(error_esf < log10(0.1)) = -3;
    esf       = esf.prom;                      % esfuerzo promedio
end


