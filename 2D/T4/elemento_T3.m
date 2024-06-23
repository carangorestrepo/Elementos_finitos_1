%%------------------------------------------------------------------------------
% NOTA: este código SOLO es apropiado para TENSION PLANA usando elementos
%       rectangulares serendípitos de 8 nodos
%-------------------------------------------------------------------------------
%
% Programa para el cálculo de los desplazamientos y las reacciones en los 
% apoyos, las deformaciones y los esfuerzos de la un sólido mediante el método 
% de los elementos finitos

%clear, clc, close all   % borro la memoria, la pantalla y las figuras

%% constantes que ayudarán en la lectura del código
X = 1; Y = 2;

r_ = [ 1, 2, 3, 4, 5, 6 ];   % GDL a retener en condensaci�n nodal
e_ = [ 7, 8];            % GDL a eliminar en condensaci�n nodal

%% se define la estructura a calcular
%nombre_archivo = {'malla_1', 'malla1'};
%nombre_archivo = {'malla_2', 'malla2'};
%nombre_archivo = {'malla_3', 'malla3'};
%nombre_archivo = {'malla_4', 'malla4'};
%nombre_archivo = {'malla_5', 'malla5'};
%archivo_xlsx = fullfile('..', nombre_archivo{1}, [nombre_archivo{2} '.xlsx']);

%% se leen las coordenadas de los nodos
%T = leer_excel(archivo_xlsx, 'xnod');
%idxNODO         = T{:,'nodo'};
%xnod(idxNODO,:) = T{:,{'x','y'}};
nno             = size(xnod,1);   % número de nodos

%% Se definen los grados de libertad
ngdl = 2*nno;        % numero de grados de libertad (dos por nodo)
gdl  = [(1:2:ngdl)' (2:2:ngdl)']; % nodos vs grados de libertad

%% se lee la matriz de conectividad (LaG) y el tipo de material del EF e
%T = leer_excel(archivo_xlsx, 'LaG_mat');
%idxEF        = T{:,'EF'};
%LaG(idxEF,:) = T{:,{'NL1','NL2','NL3','NL4','NL5','NL6','NL7','NL8'}};
%mat          = T{:,'material'};% tipo de material para cada EF
%nef          = size(LaG,1);    % numero de EFs (numero de filas de LaG)

%% se leen los materiales
%T = leer_excel(archivo_xlsx, 'prop_mat');
%E    = T{:,'E'};       % modulo de elasticidad
%nu   = T{:,'nu'};      % coeficiente de Poisson
%rho  = T{:,'rho'};     % densidad
%t    = T{:,'espesor'}; % espesor
nmat = size(E,1);      % numero de materiales

%% se leen las cargas puntuales
%T = leer_excel(archivo_xlsx, 'carga_punt');
%idxNODO = T{:,'nodo'};
%dir_cp  = T{:,'direccion'};
%cp      = T{:,'fuerza_puntual'}; % desplazamientos conocidos en los apoyos

f = zeros(ngdl,1);     % vector de fuerzas nodales equivalentes global
%{
ncp = size(cp,1); % numero de fuerzas puntuales
for i = 1:ncp
%        nodo        direccion      fuerza puntual
   f(gdl(idxNODO(i), dir_cp(i))) = cp(i);
end
%}
f(gdl(sub2ind([nno 2], idxNODOp, dir_cp))) = cp;

%% se leen algunas variables
%T        = readcell(archivo_xlsx, 'Sheet','varios','Range','B1:B9');
g        = 9.81; % aceleracion de la gravedad
U_LONG   = 'm'; % unidades de longitud
U_FUERZA = 'kN'; % unidades de fuerza
U_ESFUER = 'kPa'; % unidades de esfuerzo
ESC_UV   = 100; % factor de escala para los desplazamientos

%% Se dibuja la malla de elementos finitos
figure;
hold on;
cg = zeros(nef,2); % almacena el centro de gravedad de los EFs
for e = 1:nef
   % se dibuja el EF e
   nod_ef = LaG(e,[1:3 1]);
   plot(xnod(nod_ef,X), xnod(nod_ef,Y), 'b');
   
   % se calcula la posición del centro de gravedad del EF e
   cg(e,:) = mean(xnod(LaG(e,:),:));

   % se escribe el número del EF e
 %  text(cg(e,X), cg(e,Y), num2str(e), 'Color', 'b');
end

% en todos los nodos se dibuja un marcador y se reporta su numeración
plot(xnod(:,X), xnod(:,Y), 'ro');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
axis equal tight
title('Malla de elementos finitos');

%% Func. de forma y sus derivadas del EF rectangular serendípito de 8 nodos
% NOTA estas funciones de forma y sus derivadas se encontraron con el
% programa deduccion_funciones_forma/FF_serendipitos_Q4_Q8.m
Nforma = @(xi,eta)[( 1.0 - 9.0 * xi * eta ) * ( 1.0 - xi - eta )
           xi * ( 1.0 - 9.0 * ( 1.0 - xi - eta ) * eta )
           eta * ( 1.0 - 9.0 * ( 1.0 - xi - eta ) * xi )
           27.0 * ( 1.0 - xi - eta ) * xi * eta
           ];
Nforma3 = @(xi,eta)[( 1.0 - 9.0 * xi * eta ) * ( 1.0 - xi - eta )
           xi * ( 1.0 - 9.0 * ( 1.0 - xi - eta ) * eta )
           eta * ( 1.0 - 9.0 * ( 1.0 - xi - eta ) * xi )
           ];      
                
% derivadas de las funciones de forma con respecto a xi
dN_dxi = @(xi,eta) [ ...
     9*eta*(eta + xi - 1) + 9*eta*xi - 1
   9*eta*xi + eta*(9*eta + 9*xi - 9) + 1
                 eta*(9*eta + 18*xi - 9)
 - 27*eta*xi - eta*(27*eta + 27*xi - 27)
                    ];



%% Derivadas de N con respecto a eta
dN_deta = @(xi,eta) [ ... 
     9*xi*(eta + xi - 1) + 9*eta*xi - 1
                 xi*(18*eta + 9*xi - 9)
   9*eta*xi + xi*(9*eta + 9*xi - 9) + 1
 - 27*eta*xi - xi*(27*eta + 27*xi - 27)
                    ];


%% Parametros de la cuadratura de Gauss-Legendre
% se calculan las raices x_gl y los pesos w_gl de polinomios de Legendre
n         = 4; % orden de la cuadratura de Gauss-Legendre
%[x_gl, w_gl] = gausslegendre_quad(n_gl);
xw=TriGaussPoints(n);

x_gl = xw(:,1);
e_gl = xw(:,2);
w_gl =  xw(:,3);
n_gl = size(x_gl,1);  %# N�mero de puntos de Gauss.

%% se ensambla la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
K   = sparse(ngdl,ngdl);   % matriz de rigidez global como RALA (sparse)
inv_Kee = cell(nef, 3, 3);
Ker     = cell(nef, 3, 6);

M   = sparse(ngdl,ngdl);   % matriz de masa global como RALA (sparse)
inv_Mee = cell(nef, 3, 3);
Mer     = cell(nef, 3, 6);

N6   = cell(nef,n_gl,n_gl); % contenedor para las matrices de forma
N3   = cell(nef,n_gl,n_gl); % contenedor para las matrices de forma
B   = cell(nef,n_gl,n_gl); % contenedor para las matrices de deformacion
idx = cell(nef,1);         % indices asociados a los gdl del EF e

% se calcula la matriz constitutiva y el vector de fuerzas másicas para cada material
De = cell(nmat,1);
be = cell(nmat,1);
for i = 1:nmat
    % matriz constitutiva para TENSION PLANA
    De{i} = (E(i)/(1-nu(i)^2)) * [ 1      nu(i)  0
                                   nu(i)  1      0
                                   0      0      (1-nu(i))/2 ];

    % vector de fuerzas másicas
    be{i} = [0; -rho(i)*g];
end

% para cada elemento finito en la malla:
for e = 1:nef
   % se calculan con el siguiente ciclo las matrices de rigidez y el vector de
   % fuerzas nodales equivalentes del elemento usando las cuadraturas de GL
   Ke8 = zeros(4*2);
   fe = zeros(3*2,1);
   Me8=  zeros(4*2,1);
   det_Je = zeros(n_gl,1); % matriz para almacenar los jacobianos

   % se determinan las coordenadas de los nodos el EF e
   xe = [xnod(LaG(e,:),X);cg(e,1)];
   ye = [xnod(LaG(e,:),Y);cg(e,2)];

   for p = 1:n_gl

         xi_gl  = x_gl(p);
         eta_gl = e_gl(p);
         
         % Se evaluan las funciones de forma y sus derivadas 
         % en los puntos de integracion de Gauss-Legendre
         NNforma  = Nforma (xi_gl, eta_gl);
         NNforma3  = Nforma3 (xi_gl, eta_gl);
         ddN_dxi  = dN_dxi (xi_gl, eta_gl);
         ddN_deta = dN_deta(xi_gl, eta_gl);
         
         dx_dxi  = sum(ddN_dxi  .* xe);   dy_dxi  = sum(ddN_dxi  .* ye);
         dx_deta = sum(ddN_deta .* xe);   dy_deta = sum(ddN_deta .* ye);
         
         % Se ensambla la matriz Jacobiana del elemento
         Je = [ dx_dxi   dy_dxi
                dx_deta  dy_deta ];
            
         % Se calcula el determinante del Jacobiano
         det_Je(p) = det(Je);
         
         % las matrices de forma y de deformación se evalúan y se ensamblan
         % en el punto de Gauss         
         N4{e,p} = zeros(2,2*4);
         N3{e,p} = zeros(2,2*3);
         B{e,p}  = zeros(3,2*4);
         
         for i = 1:3
            % Se ensambla la matriz de funciones de forma N
            N3{e,p}(:,[2*i-1 2*i]) = [ NNforma3(i)  0         
                                        0           NNforma3(i) ];
         end
         for i = 1:4
            % Se ensambla la matriz de funciones de forma N
            N4{e,p}(:,[2*i-1 2*i]) = [ NNforma(i)  0         
                                        0         NNforma(i) ];
            % Se ensambla la matriz de deformacion del elemento B
            dNi_dx = (+dy_deta*ddN_dxi(i) - dy_dxi*ddN_deta(i))/det_Je(p);
            dNi_dy = (-dx_deta*ddN_dxi(i) + dx_dxi*ddN_deta(i))/det_Je(p);
            B{e,p}(:,[2*i-1 2*i]) = [ dNi_dx       0        
                                             0  dNi_dy   
                                      dNi_dy  dNi_dx ];
         end
         % se ensamblan la matriz de rigidez del EF e y el vector de fuerzas
         % nodales equivalentes del EF e asociado a la fuerza másica         
         Ke8 = Ke8 + B{e,p}'*De{mat(e)}*B{e,p} * det_Je(p)*t(mat(e))*w_gl(p);
         fe =   fe + N3{e,p}'*be{mat(e)}        * det_Je(p)*t(mat(e))*w_gl(p);
         Me8 = Me8 + N4{e,p}'*rho*N4{e,p}        * det_Je(p)*t(mat(e))*w_gl(p);
    end

   % se determina si hay puntos con jacobiano negativo, en caso tal se termina
   % el programa y se reporta   
   if any(any(det_Je <= 0))
      error('Hay puntos con det_Je negativo en el EF %d.\n', e);
   end
     % se condensan los GDL jer�rquicos u5, v5, u6, v6
    Krr = Ke8(r_,r_); 
    Ker{e}= Ke8(e_,r_);
    Kre = Ke8(r_,e_);
    inv_Kee{e} = (Ke8(e_,e_))^(-1);
    
    Ke = Krr - Kre * inv_Kee{e} * Ker{e};
    
    Mrr = Me8(r_,r_); 
    Mer{e}= Me8(e_,r_);
    MKre = Ke8(r_,e_);
    inv_Mee{e} = (Me8(e_,e_))^(-1);
    
    Me = Mrr - MKre * inv_Mee{e} * Mer{e};
   % y se ensambla la matriz de rigidez del elemento y el vector de fuerzas
   % nodales del elemento en sus correspondientes GDL 
   idx{e}           = reshape(gdl(LaG(e,:),:)', 1, 3*2);
   %idx{e} = [ gdl(LaG(e,1),:)  gdl(LaG(e,2),:) ...
   %          gdl(LaG(e,3),:)  gdl(LaG(e,4),:) ...
   %          gdl(LaG(e,5),:)  gdl(LaG(e,6),:) ...
   %          gdl(LaG(e,7),:)  gdl(LaG(e,8),:) gdl(LaG(e,9),:) ];
   
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke;
   f(idx{e},:)      = f(idx{e},:)      + fe;
   M(idx{e},idx{e}) = M(idx{e},idx{e}) + Me;
end

%% Muestro la configuracion de la matriz K (K es rala)
figure
spy(K);
title('Los puntos representan los elementos diferentes de cero');

%% se leen las cargas distribuidas
%T       = leer_excel(archivo_xlsx, 'carga_distr');
%idxELEM = T{:,'elemento'};
%nodoijk = T{:,{'nodo_i','nodo_j','nodo_k'}};
%carga   = T{:,{'tix','tiy', 'tjx','tjy', 'tkx','tky'}};
nlcd    = size(carga,1); % numero de lados con carga distribuida

%% Relacion de las cargas superficiales (vector ft)
ft = sparse(ngdl,1); % fuerzas nodales equivalentes de cargas superficiales
for i = 1:nlcd
   e     = idxELEM(i);
   LaG_e = LaG(e,:);
   fte   = t2ft_T3(xnod(LaG_e,[X Y]), LaG_e, nodoijk(i,:), carga(i,:), t(mat(e)));
   ft(idx{e},:) = ft(idx{e},:) + fte;
end

% Agrego al vector de fuerzas nodales equivalentes las fuerzas
% superficiales calculadas
f = f + ft;

%% se leen las constantes de balastro k (cimentacion elastica de Winkler)
%T       = leer_excel(archivo_xlsx, 'kWinkler');
%idxELEM = T{:,'elemento'};
%nodoijk = T{:,{'nodo_i','nodo_j','nodo_k'}};
%kwinkl  = T{:,{'kix','kiy', 'kjx','kjy', 'kkx','kky'}};
%nlkW    = size(kwinkl,1); % numero de lados con cimentacion elastica

%% Cálculo de las rigideces asociadas a la cimentación elástica
%for i = 1:nlkW
%   e = idxELEM(i);
%   LaG_e = LaG(e,1:8);
%   He = Hwinkler_8(xnod(LaG_e,[X Y]), LaG_e, nodoijk(i,:), kwinkl(i,:), t(mat(e)));
%   K(idx{e},idx{e}) = K(idx{e},idx{e}) + He;
%end

%% se definen los apoyos y sus desplazamientos
%T = leer_excel(archivo_xlsx, 'restric');
%idxNODOr  = T{:,'nodo'};
%dir_desp = T{:,'direccion'};

% grados de libertad del desplazamiento conocidos  
%{
ngdl_res = size(ac,1); % numero de grados de libertad restringidos
for i = 1:ngdl_res
%             nodo        direccion
   c(i) = gdl(idxNODO(i), dir_desp(i));
end
%}
c = gdl(sub2ind([nno 2], idxNODOr, dir_desp));

% desplazamientos conocidos en los apoyos
%ac = T{:,'desplazamiento'};

% grados de libertad del desplazamiento desconocidos  
d = setdiff(1:ngdl,c)';

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
   line(xnod(LaG(e,[1,2,3,1]),X), xnod(LaG(e,[1,2,3,1]),Y), 'Color','b'); % original
   line(xdef(LaG(e,[1,2,3,1]),X), xdef(LaG(e,[1,2,3,1]),Y), 'Color','r'); % deformada
end
xlabel(['Eje X [' U_LONG ']']);
ylabel(['Eje Y [' U_LONG ']']);
axis equal tight;
legend('Posicion original','Posicion deformada','Location', 'SouthOutside');
title(sprintf('Deformada escalada %d veces', ESC_UV));

%% Se calcula para cada elemento las deformaciones y los esfuerzos
def = cell(nef,n_gl,3);
esf = cell(nef,n_gl,3);

for e = 1:nef
   % desplazamientos de los gdl del elemento e
      % desplazamientos de los gdl del elemento e
   ar = a(idx{e});
   ae = [ar;
         inv_Kee{e}  * Ker{e} * ar];

   for pp = 1:n_gl
         def{e,pp} = B{e,pp}*ae;           % calculo las deformaciones
         esf{e,pp} = De{mat(e)}*def{e,pp}; % calculo los esfuerzos
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

ez   = -(nu/E)*(sx+sy);              % deformaciones ez
tmax = sqrt(((sx-sy)/2).^2+txy.^2);  % esfuerzo cortante maximo
s1   = (sx+sy)/2 + tmax;             % esfuerzo normal maximo
s2   = (sx+sy)/2 - tmax;             % esfuerzo normal minimo
ang  = 0.5*atan2(2*txy, sx-sy);      % angulo de inclinacion de s1

s3   = zeros(nno,1);
sv   = sqrt(((s1-s2).^2 + (s2-s3).^2 + (s1-s3).^2)/2); % von Mises

sxx=reshape(sx,[Ny,Nx]);

syy=reshape(sy,[Ny,Nx]);
txyy=reshape(txy,[Ny,Nx]);
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
%{
%% se grafican las deformaciones
figure
subplot(1,4,1); plot_def_esf(xnod, LaG, ex,  '\epsilon_x')
subplot(1,4,2); plot_def_esf(xnod, LaG, ey,  '\epsilon_y')
subplot(1,4,3); plot_def_esf(xnod, LaG, ez,  '\epsilon_z')
subplot(1,4,4); plot_def_esf(xnod, LaG, gxy, '\gamma_{xy} [rad]')

%% se grafican los esfuerzos
figure
subplot(1,3,1); plot_def_esf(xnod, LaG, sx,  ['\sigma_x [' U_ESFUER ']'])
subplot(1,3,2); plot_def_esf(xnod, LaG, sy,  ['\sigma_y [' U_ESFUER ']'])
subplot(1,3,3); plot_def_esf(xnod, LaG, txy, ['\tau_{xy} [' U_ESFUER ']'])

%% se grafican los errores en los esfuerzos
figure
subplot(1,3,1); plot_def_esf(xnod, LaG, error_sx,  ['Error \sigma_x [' U_ESFUER ']'])
subplot(1,3,2); plot_def_esf(xnod, LaG, error_sy,  ['Error \sigma_y [' U_ESFUER ']'])
subplot(1,3,3); plot_def_esf(xnod, LaG, error_txy, ['Error \tau_{xy} [' U_ESFUER ']'])

%% se grafican los esfuerzos principales y el esfuerzo cortante maximo
figure
subplot(1,3,1); plot_def_esf(xnod, LaG, s1,   ['(\sigma_1)_{xy} [' U_ESFUER ']'], { ang })
subplot(1,3,2); plot_def_esf(xnod, LaG, s2,   ['(\sigma_2)_{xy} [' U_ESFUER ']'], { ang+pi/2 })
subplot(1,3,3); plot_def_esf(xnod, LaG, tmax, ['\tau_{max} [' U_ESFUER ']'],      { ang+pi/4, ang-pi/4 })

%% se grafican los esfuerzos de von Mises
figure
plot_def_esf(xnod, LaG, sv, ['Esfuerzos de von Mises [' U_ESFUER ']']);
%}
%% se exportan los resultados a GiD/Paraview
% Pasando los esfuerzos ya promediados:
%export_to_GiD('c5_ejemplo_a',xnod,LaG,a,q,[sx sy sz txy txz tyz]);

% Pasando los puntos de Gauss [RECOMENDADO] !!!
% export_to_GiD('c5_ejemplo_b',xnod,LaG,a,q,esf);                    

%%
return; % bye, bye!
%{
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
    hold on; 
    colorbar;
    
    nef = size(LaG, 1);    
    for e = 1:nef
       data.numEF = e; 
       data.LaG_e = LaG(e,:);
       fill(xnod(LaG(e,:),X), xnod(LaG(e,:),Y), variable(LaG(e,:)), ...
           'UserData', data);
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
                'Marker','.');            % y en el punto (x,y) poner un punto '.'
            
            % la misma flecha girada 180 grados
            quiver(xnod(:,X),xnod(:,Y),...             
                norma.*cos(angulos{i}+pi), norma.*sin(angulos{i}+pi),... 
                esc,'k', 'ShowArrowHead','off', 'LineWidth',2, 'Marker','.');                    
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
    
        % se busca el punto más cercano
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
    A =[0.126340726488392,-0.638559587411924,-0.638559587411923,1.87365927351158,0.138559587411936,0.138559587411936;
        -0.638559587411907,-0.638559587411945,0.126340726488378,0.138559587411937,0.138559587411936,1.87365927351160;
        -0.638559587411908,0.126340726488377,-0.638559587411944,0.138559587411937,1.87365927351160,0.138559587411936];
    switch tipo_esf
        case {'sx',  'ex'},  num_esf = 1;
        case {'sy',  'ey'},  num_esf = 2;
        case {'txy', 'gxy'}, num_esf = 3;
        otherwise,           error('Opcion no soportada');
    end
    
    % se hace la extrapolación de los esfuerzos en cada EF a partir de las 
    % lecturas en los puntos de Gauss
    for e = 1:nef
        esfa=[esfuerzo{e,:,1}];
        esf_EF_e = A * esfa(num_esf,:)';    
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
