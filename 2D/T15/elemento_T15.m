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
   nod_ef = LaG(e,[1,4,5,6,2,7,8,9,3,10,11,12,1]);
   plot(xnod(nod_ef,X), xnod(nod_ef,Y), 'b');
   
   % se calcula la posición del centro de gravedad del EF e
   cg(e,:) = mean(xnod(LaG(e,:),:));

   % se escribe el número del EF e
   text(cg(e,X), cg(e,Y), num2str(e), 'Color', 'b');
end

% en todos los nodos se dibuja un marcador y se reporta su numeración
plot(xnod(:,X), xnod(:,Y), 'ro');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
axis equal tight
title('Malla de elementos finitos');

%% Func. de forma y sus derivadas del EF rectangular serendípito de 8 nodos
% NOTA estas funciones de forma y sus derivadas se encontraron con el
% programa deduccion_funciones_forma/FF_serendipitos_T10.m
 Nforma =  @(xi,eta)[ ((eta + xi - 1)*(22*eta + 22*xi + 48*(eta + xi - 1)^2 + 32*(eta + xi - 1)^3 - 19))/3
                                                                   (xi*(32*xi^3 - 48*xi^2 + 22*xi - 3))/3
                                                               (eta*(32*eta^3 - 48*eta^2 + 22*eta - 3))/3
                                    -(xi*(16*eta + 16*xi - 16)*(6*eta + 6*xi + 8*(eta + xi - 1)^2 - 5))/3
                                                      xi*(4*xi - 1)*(4*eta + 4*xi - 3)*(4*eta + 4*xi - 4)
                                                        -(xi*(8*xi^2 - 6*xi + 1)*(16*eta + 16*xi - 16))/3
                                                                        (16*eta*xi*(8*xi^2 - 6*xi + 1))/3
                                                                          4*eta*xi*(4*eta - 1)*(4*xi - 1)
                                                                      (16*eta*xi*(8*eta^2 - 6*eta + 1))/3
                                                     -(eta*(8*eta^2 - 6*eta + 1)*(16*eta + 16*xi - 16))/3
                                                    eta*(4*eta - 1)*(4*eta + 4*xi - 3)*(4*eta + 4*xi - 4)
                                   -(eta*(16*eta + 16*xi - 16)*(6*eta + 6*xi + 8*(eta + xi - 1)^2 - 5))/3
                                                          eta*xi*(4*eta + 4*xi - 3)*(32*eta + 32*xi - 32)
                                                                 -eta*xi*(4*xi - 1)*(32*eta + 32*xi - 32)
                                                                -eta*xi*(4*eta - 1)*(32*eta + 32*xi - 32)];

%% Derivadas de N con respecto a xi
dN_dxi =@(xi,eta) [
(128*eta^3)/3 + 128*eta^2*xi - 80*eta^2 + 128*eta*xi^2 - 160*eta*xi + (140*eta)/3 + (128*xi^3)/3 - 80*xi^2 + (140*xi)/3 - 25/3
                                                                                      (128*xi^3)/3 - 48*xi^2 + (44*xi)/3 - 1
                                                                                                                           0
- (128*eta^3)/3 - 256*eta^2*xi + 96*eta^2 - 384*eta*xi^2 + 384*eta*xi - (208*eta)/3 - (512*xi^3)/3 + 288*xi^2 - (416*xi)/3 + 16
                            128*eta^2*xi - 16*eta^2 + 384*eta*xi^2 - 288*eta*xi + 28*eta + 256*xi^3 - 384*xi^2 + 152*xi - 12
                                         64*eta*xi - (224*xi)/3 - (16*eta)/3 - 128*eta*xi^2 + 224*xi^2 - (512*xi^3)/3 + 16/3
                                                                                            (16*eta*(24*xi^2 - 12*xi + 1))/3
                                                                                                4*eta*(4*eta - 1)*(8*xi - 1)
                                                                                            (16*eta*(8*eta^2 - 6*eta + 1))/3
                                                                                           -(16*eta*(8*eta^2 - 6*eta + 1))/3
                                                                                        4*eta*(4*eta - 1)*(8*eta + 8*xi - 7)
                                                          -(16*eta*(24*eta^2 + 48*eta*xi - 36*eta + 24*xi^2 - 36*xi + 13))/3
                                                                  32*eta*(4*eta^2 + 16*eta*xi - 7*eta + 12*xi^2 - 14*xi + 3)
                                                                              -32*eta*(8*eta*xi - 10*xi - eta + 12*xi^2 + 1)
                                                                                        -32*eta*(4*eta - 1)*(eta + 2*xi - 1)];

%% Derivadas de N con respecto a eta
dN_deta = @(xi,eta) [
(128*eta^3)/3 + 128*eta^2*xi - 80*eta^2 + 128*eta*xi^2 - 160*eta*xi + (140*eta)/3 + (128*xi^3)/3 - 80*xi^2 + (140*xi)/3 - 25/3
                                                                                                                           0
                                                                                   (128*eta^3)/3 - 48*eta^2 + (44*eta)/3 - 1
                                                           -(16*xi*(24*eta^2 + 48*eta*xi - 36*eta + 24*xi^2 - 36*xi + 13))/3
                                                                                          4*xi*(4*xi - 1)*(8*eta + 8*xi - 7)
                                                                                              -(16*xi*(8*xi^2 - 6*xi + 1))/3
                                                                                               (16*xi*(8*xi^2 - 6*xi + 1))/3
                                                                                                 4*xi*(8*eta - 1)*(4*xi - 1)
                                                                                           (16*xi*(24*eta^2 - 12*eta + 1))/3
                                       64*eta*xi - (16*xi)/3 - (224*eta)/3 - 128*eta^2*xi + 224*eta^2 - (512*eta^3)/3 + 16/3
                           256*eta^3 + 384*eta^2*xi - 384*eta^2 + 128*eta*xi^2 - 288*eta*xi + 152*eta - 16*xi^2 + 28*xi - 12
- (512*eta^3)/3 - 384*eta^2*xi + 288*eta^2 - 256*eta*xi^2 + 384*eta*xi - (416*eta)/3 - (128*xi^3)/3 + 96*xi^2 - (208*xi)/3 + 16
                                                                   32*xi*(12*eta^2 + 16*eta*xi - 14*eta + 4*xi^2 - 7*xi + 3)
                                                                                          -32*xi*(4*xi - 1)*(2*eta + xi - 1)
                                                                              -32*xi*(8*eta*xi - xi - 10*eta + 12*eta^2 + 1)];

%% Parametros de la cuadratura de Gauss-Legendre
% se calculan las raices x_gl y los pesos w_gl de polinomios de Legendre
n         = 6; % orden de la cuadratura de Gauss-Legendre
%[x_gl, w_gl] = gausslegendre_quad(n_gl);
xw=TriGaussPoints(n);

x_gl = xw(:,1);
e_gl = xw(:,2);
w_gl =  xw(:,3);
n_gl = size(x_gl,1);  %# N�mero de puntos de Gauss.

%% se ensambla la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
K   = zeros(ngdl,ngdl);   % matriz de rigidez global como RALA (sparse)
M   = sparse(ngdl,ngdl);   % matriz de masa global como RALA (sparse)
N   = cell(nef,n_gl,3,2*15); % contenedor para las matrices de forma
B   = cell(nef,n_gl,3,2*15); % contenedor para las matrices de deformacion
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
   Ke = zeros(15*2);
   fe = zeros(15*2,1);
   Me=  zeros(15*2,1);
   det_Je = zeros(n_gl,1); % matriz para almacenar los jacobianos

   % se determinan las coordenadas de los nodos el EF e
   xe = xnod(LaG(e,:),X);
   ye = xnod(LaG(e,:),Y);

   for p = 1:n_gl

         xi_gl  = x_gl(p);
         eta_gl = e_gl(p);
         
         % Se evaluan las funciones de forma y sus derivadas 
         % en los puntos de integracion de Gauss-Legendre
         NNforma  = Nforma (xi_gl, eta_gl);
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
         N{e,p} = zeros(2, 2*15);
         B{e,p} = zeros(3, 2*15);
         for i = 1:15
            % Se ensambla la matriz de funciones de forma N
            N{e,p}(:,[2*i-1 2*i]) = [ NNforma(i)  0         
                                        0           NNforma(i) ];
         
            % Se ensambla la matriz de deformacion del elemento B
            dNi_dx = (+dy_deta*ddN_dxi(i) - dy_dxi*ddN_deta(i))/det_Je(p);
            dNi_dy = (-dx_deta*ddN_dxi(i) + dx_dxi*ddN_deta(i))/det_Je(p);
            B{e,p}(:,[2*i-1 2*i]) = [ dNi_dx       0        
                                             0  dNi_dy   
                                      dNi_dy  dNi_dx ];
         end

         % se ensamblan la matriz de rigidez del EF e y el vector de fuerzas
         % nodales equivalentes del EF e asociado a la fuerza másica         
         Ke = Ke + B{e,p}'*De{mat(e)}*B{e,p} * det_Je(p)*t(mat(e))*w_gl(p);
         fe = fe + N{e,p}'*be{mat(e)}        * det_Je(p)*t(mat(e))*w_gl(p);
         Me = Me + N{e,p}'*rho *N{e,p}       * det_Je(p)*t(mat(e))*w_gl(p);
    end

   % se determina si hay puntos con jacobiano negativo, en caso tal se termina
   % el programa y se reporta   
   if any(any(det_Je <= 0))
      error('Hay puntos con det_Je negativo en el EF %d.\n', e);
   end

   % y se ensambla la matriz de rigidez del elemento y el vector de fuerzas
   % nodales del elemento en sus correspondientes GDL 
   idx{e}           = reshape(gdl(LaG(e,:),:)', 1, 15*2);
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
   fte   = t2ft_T15(xnod(LaG_e,[X Y]), LaG_e, nodoijk(i,:), carga(i,:), t(mat(e)));
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
   line(xnod(LaG(e,[1,4,5,6,2,7,8,9,3,10,11,12,1]),X), xnod(LaG(e,[1,4,5,6,2,7,8,9,3,10,11,12,1]),Y), 'Color','b'); % original
   line(xdef(LaG(e,[1,4,5,6,2,7,8,9,3,10,11,12,1]),X), xdef(LaG(e,[1,4,5,6,2,7,8,9,3,10,11,12,1]),Y), 'Color','r'); % deformada
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
   ae = a(idx{e});
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
    A = [0.960092769143901,-0.407042414326935,-0.407042414326933,2.30338198189438,-0.189111457117860,-0.189111457117858,0.134383168337343,0.657323464607180,-1.32729013701887,0.134383168337342,-1.32729013701887,0.657323464607185;
        -16.1795095672505,7.87322104293654,8.45229646480399,-2.54884946308189,1.04970883586833,3.42429969487222,-2.19188615649463,-5.43448404227668,7.09078669469698,-3.80354518702247,7.56029117503839,-4.29232949209030;
        -16.1795095672505,8.45229646480403,7.87322104293652,-2.54884946308189,3.42429969487224,1.04970883586832,-3.80354518702249,-4.29232949209027,7.56029117503844,-2.19188615649462,7.09078669469696,-5.43448404227672;
        -1.31311182574608,0.654757241422352,0.234460647711091,0.139442524572790,0.0159739229725569,0.0923874673002777,-0.193371860326651,-0.365236573760065,0.188785533893397,-0.0667419834484363,1.73109419099705,-0.118439285588279;
        -2.18368756934197,1.39356491361819,0.176261426743583,-0.403153281271127,0.141917251736234,-0.0500784398904940,-0.428910764405423,0.528777107921811,0.862721505686405,-0.0941195896631753,1.64563880868672,-0.588931369820761;
        -5.76727372779353,3.30213478615734,2.04124500502355,-0.805539316776360,0.412051299319412,0.641291932302574,-0.964370033320805,-0.0358879087745948,2.54617086405574,-0.991727905717610,2.29561428153792,-1.67370927601364;
        5.49951623463161,-2.54140849191087,-3.38200167933338,0.740804212055950,-0.322913692932887,-0.170086604277439,0.619954489712391,1.35612436427419,-2.34590175417990,2.38170256259102,-2.22303117164770,1.38724153101703;
        -5.90361509013581,2.64487693057782,2.64487693057780,-0.949822729350195,0.319254129962406,0.319254129962402,-0.144713272911792,-1.58981798423268,2.69711910634727,-0.144713272911786,2.69711910634726,-1.58981798423269;
        5.49951623463162,-3.38200167933340,-2.54140849191087,0.740804212055952,-0.170086604277444,-0.322913692932884,2.38170256259102,1.38724153101702,-2.22303117164772,0.619954489712388,-2.34590175417990,1.35612436427420;
        -5.76727372779351,2.04124500502355,3.30213478615732,-0.805539316776355,0.641291932302575,0.412051299319407,-0.991727905717612,-1.67370927601363,2.29561428153792,-0.964370033320798,2.54617086405572,-0.0358879087745966;
        -2.18368756934195,0.176261426743578,1.39356491361818,-0.403153281271123,-0.0500784398904946,0.141917251736231,-0.0941195896631738,-0.588931369820752,1.64563880868672,-0.428910764405419,0.862721505686393,0.528777107921814;
        -1.31311182574607,0.234460647711088,0.654757241422346,0.139442524572793,0.0923874673002775,0.0159739229725551,-0.0667419834484353,-0.118439285588274,1.73109419099704,-0.193371860326649,0.188785533893389,-0.365236573760064;
        0.994381139870656,0.00485147204015903,0.00485147204015901,-0.000987426287318188,0.000135900547958122,0.000135900547958114,-0.00119563407257569,-0.000869041525523306,0.000380446218312987,-0.00119563407257568,0.000380446218312971,-0.000869041525523327;
        -0.0202313993022238,0.0143762403663800,1.00993924288682,-0.00361677835962526,0.00156091272896122,0.00134024043926201,-0.00367591265089046,-0.00762837457427984,0.00962005784538494,-0.00528174104069501,0.0103557515515545,-0.00675823989064652;
        -0.0202313993022328,1.00993924288682,0.0143762403663846,-0.00361677835962667,0.00134024043926277,0.00156091272896182,-0.00528174104069638,-0.00675823989064886,0.0103557515515587,-0.00367591265089202,0.00962005784538874,-0.00762837457428235];

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
