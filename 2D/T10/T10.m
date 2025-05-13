clc
clear
close all
%% se leen algunas variables
%T        = readcell(archivo_xlsx, 'Sheet','varios','Range','B1:B9');
g        = 9.81;    % [m/s²]   aceleración de la gravedad
U_LONG   = 'm';     % unidades de longitud
U_FUERZA = 'kN';    % unidades de fuerza
U_ESFUER = 'kN/m2'; % unidades de esfuerzo
ESC_UV   = 50;   % factor de escala para los desplazamientos
% 
%% constantes que ayudarán en la lectura del código

X = 1;
Y = 2;
NL1=1; NL2=2; NL3=3; NL4=4; NL5=5; NL6 =6; NL7 =7;NL8 =8; NL9 =9;NL10 =10;

 %% seleccione la malla a emplear:
%nombre_archivo = 'malla1'    # EJEMPLO CLASE
%nombre='nombre_archivo';
nombre_archivo = 'malla';%    %EJEMPLO CLASE
xnodxy = xlsread(nombre_archivo,1);% xnod;
LaG_mat = xlsread(nombre_archivo,2);% LaG;
restric = xlsread(nombre_archivo,3);% restric;
carga_distr = xlsread(nombre_archivo,4);% carga_distr;
carga_punt=xlsread(nombre_archivo,5);% carga_punt;
prop_mat=xlsread(nombre_archivo,6);% prop_mat;
% %% posición de los nodos:
% xnod: fila=número del nodo, columna=coordenada X=0 o Y=1
xnod=xnodxy(:,2:3);
nno  = size(xnod,1);%   %# número de nodos (número de filas de la matriz xnod)

%% definición de los grados de libertad

ngdl = 2*nno;            % número de grados de libertad por nodo = [X, Y]
gdl  = [(1:2:(nno*2))',(1:2:(nno*2))'+1]; % nodos vs grados de libertad

%% definición de elementos finitos con respecto a nodos
% LaG: fila=número del elemento, columna=número del nodo local

LaG = LaG_mat(:,2:11);
nef = size(LaG,1);     % número de EFs (número de filas de la matriz LaG)

%% definición de los materiales
mat  = LaG_mat(:,12);
Ee   = prop_mat(:,2);     % [KPa]     módulo de elasticidad del sólido
nue  = prop_mat(:,3);     % [-]      coeficiente de Poisson
rhoe = prop_mat(:,4);     % [T/m³]  densidad
te   = prop_mat(:,5);     % [m]      espesor
nmat =size(Ee,1);         % número de materiales

%% relación de cargas puntuales
%cp  = carga_punt(:,
ncp = size(carga_punt,1);    % número de cargas puntuales
f   = zeros(ngdl,1);           % vector de fuerzas nodales equivalentes global
for i=1: ncp
   f(carga_punt(i,1)*2-carga_punt(i,2)+1,1)=carga_punt(i,3);
end

%% Se dibuja la malla de elementos finitos
cg = zeros(nef,2);%  % almacena el centro de gravedad de los EF
figure()
hold on
for e =1:nef
    % se dibujan las aristas
    nod_ef = LaG(e, [NL1, NL4, NL5, NL2, NL6, NL7, NL3, NL8, NL9, NL1]);
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

%% Funciones de forma y sus derivadas del elemento triangular de 10 nodos:

    %''' Evalúa las funciones de forma en un punto determinado de coordenadas
    %    (xi, eta).
    %    Funciones de forma obtenidas del programa FF_triangulo.py.
    %'''
    %L1 = (1 - xi - eta);
    %L2 = xi;
    %L3 = eta;
 Nforma =  @(xi,eta)[(1 - xi - eta)*(4.5*(1 - xi - eta)^2 - 4.5*(1 - xi - eta) + 1.0)
                 xi *(4.5*xi ^2 - 4.5*xi  + 1.0)
                 eta*(4.5*eta^2 - 4.5*eta + 1.0)
                 (1 - xi - eta)*xi *(13.5*(1 - xi - eta) - 4.5)
                 (1 - xi - eta)*xi *(13.5*xi  - 4.5)
                 xi *eta*(13.5*xi  - 4.5)
                 xi *eta*(13.5*eta - 4.5)
                 (1 - xi - eta)*eta*(13.5*eta - 4.5)
                 (1 - xi - eta)*eta*(13.5*(1 - xi - eta) - 4.5)
                 27.0*(1 - xi - eta)*xi *eta];
% derivadas de las funciones de forma con respecto a xi
    %L2 = xi;
    %L3 = eta;
    dN_dxi =@(xi,eta) [
            -13.5*xi ^2 - 27.0*xi *eta + 18.0*xi  - 13.5*eta^2 + 18.0*eta - 5.5
            13.5*xi ^2 - 9.0*xi  + 1.0
            0
            40.5*xi ^2 + 54.0*xi *eta - 45.0*xi  + 13.5*eta^2 - 22.5*eta + 9.0
            -40.5*xi ^2 - 27.0*xi *eta + 36.0*xi  + 4.5*eta - 4.5
            eta*(27.0*xi  - 4.5)
            eta*(13.5*eta - 4.5)
            eta*(-13.5*eta + 4.5)
            eta*(27.0*xi  + 27.0*eta - 22.5)
            27.0*eta*(-2*xi  - eta + 1)
            ];
% derivadas de N con respecto a eta
    %xi  = xi;
    %L3 = eta;
    dN_deta = @(xi,eta) [
            -13.5*xi ^2 - 27.0*xi *eta + 18.0*xi  - 13.5*eta^2 + 18.0*eta - 5.5
            0
            13.5*eta^2 - 9.0*eta + 1.0
            xi *(27.0*xi  + 27.0*eta - 22.5)
            xi *(-13.5*xi  + 4.5)
            xi *(13.5*xi  - 4.5)
            xi *(27.0*eta - 4.5)
            -27.0*xi *eta + 4.5*xi  - 40.5*eta^2 + 36.0*eta - 4.5
            13.5*xi ^2 + 54.0*xi *eta - 22.5*xi  + 40.5*eta^2 - 45.0*eta + 9.0
            27.0*xi *(-xi  - 2*eta + 1)
            ];
n = 4;  %# orden de la cuadratura de Gauss-Legendre (para triángulos)
xw = TriGaussPoints(n);

x_gl = xw(:,1); 
eta_gl=xw(:,2); 
w_gl=xw(:,3);
n_gl = size(x_gl,1);  %# Número de puntos de Gauss.

 %% Ensamblaje la matriz de rigidez global y el vector de fuerzas másicas
%    nodales equivalentes global

% se inicializan la matriz de rigidez global y los espacios en memoria que
%  almacenarán las matrices de forma y de deformación

K       = zeros(ngdl, ngdl);         % matriz de rigidez global
N = cell(nef,n_gl,n_gl,2,2*10); % matriz de forma en cada punto de GL
B = cell(nef,n_gl,n_gl,3,2*10); % matriz de deformaciones en cada punto de GL
idx = cell(nef,1);                 % indices asociados a los gdl del EF e
recorrido=1:2:19;
% matriz constitutiva del elemento para TENSION PLANA
De =  cell(nmat,1);
be = cell(nmat,1);

for i = 1:nmat
    % matriz constitutiva para TENSION PLANA
    De{i} = (Ee(i)/(1-nue(i)^2)) * [ 1     nue(i)  0
                                   nue(i)  1       0
                                   0       0      (1-nue(i))/2 ];
    % vector de fuerzas mÃ¡sicas
    be{i} = [0; -rhoe(i)*g];
end
% para cada elemento finito en la malla:
for e = 1:nef
   % se calculan con el siguiente ciclo las matrices de rigidez y el vector de
   % fuerzas nodales equivalentes del elemento usando las cuadraturas de GL
   Ke = zeros(2*10);
   fe = zeros(2*10,1);
   det_Je = zeros(n_gl,1); % matriz para almacenar los jacobianos
   
   % se determinan las coordenadas de los nodos el EF e
   xe = xnod(LaG(e,:),X);
   ye = xnod(LaG(e,:),Y);
   for i = 1:n_gl
         xi  = x_gl(i);
         eta = eta_gl(i);
         
         % Se evaluan las funciones de forma y sus derivadas 
         % en los puntos de integracion de Gauss-Legendre
         NNforma  = Nforma (xi, eta);
         ddN_dxi  = dN_dxi (xi, eta);
         ddN_deta = dN_deta(xi, eta);
         
         dx_dxi  = sum(ddN_dxi .* xe);   dy_dxi  = sum(ddN_dxi  .* ye);
         dx_deta = sum(ddN_deta .* xe);   dy_deta = sum(ddN_deta .* ye);
         
         % Se ensambla la matriz Jacobiana del elemento
         Je = [ dx_dxi   dy_dxi
                dx_deta  dy_deta ];
            
         % Se calcula el determinante del Jacobiano
         det_Je(i) = det(Je);
         
        % las matrices de forma y de deformación se evalúan y se ensamblan
        % en el punto de Gauss     
         Npq = zeros(2, 2*10);
         Bpq = zeros(3, 2*10);
         % Se ensambla la matriz de funciones de forma N
         Npq(1,recorrido)=NNforma;
         Npq(2,recorrido+1)=NNforma;
         
         dNi_dx=(+dy_deta*ddN_dxi - dy_dxi*ddN_deta)/det_Je(i);
         dNi_dy=(-dx_deta*ddN_dxi + dx_dxi*ddN_deta)/det_Je(i);
         
         Bpq(1,recorrido) = dNi_dx;
         Bpq(2,recorrido+1) = dNi_dy;
         Bpq(3,recorrido) = dNi_dy;
         Bpq(3,recorrido+1) = dNi_dx;
    
         N{e,i}=Npq;
         B{e,i}=Bpq;
         % se ensamblan la matriz de rigidez del EF e y el vector de fuerzas
         % nodales equivalentes del EF e asociado a la fuerza mÃ¡sica         
         Ke = Ke + B{e,i}'*De{mat(e)}*B{e,i} * det_Je(i)*te(mat(e))*w_gl(i);
         fe =   fe + N{e,i}'*be{mat(e)}      * det_Je(i)*te(mat(e))*w_gl(i);
   end
   
   % se determina si hay puntos con jacobiano negativo, en caso tal se termina
   % el programa y se reporta   
   if any(any(det_Je <= 0))
      error('Hay puntos con det_Je negativo en el EF %d.\n', e);
   end
   % y se ensambla la matriz de rigidez del elemento y el vector de fuerzas
   % nodales del elemento en sus correspondientes GDL 
   idx{e}=reshape([LaG(e, [NL1, NL2, NL3, NL4, NL5, NL6, NL7,NL8, NL9, NL10])*2-1;LaG(e, [NL1, NL2, NL3, NL4, NL5, NL6, NL7,NL8, NL9, NL10])*2],1,20);
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke;
   f(idx{e},:)      = f(idx{e},:)      + fe;
end
%% Muestro la configuracion de la matriz K (K es rala)
figure
spy(K);
title('Los puntos representan los elementos diferentes de cero');

%% Cálculo de las cargas nodales equivalentes de las cargas distribuidas:
cd   = carga_distr;
nlcd = size(cd,1);     % número de lados con carga distribuida
ft   = zeros(ngdl,1);  % fuerzas nodales equivalentes de cargas superficiales

for i =1:nlcd
   e     = cd(i,1);
   lado  = cd(i,2);
   carga = cd(i,3:10);
   fte = t2ft_T10(xnod(LaG(e,:),:), lado, carga, te(mat(e)));
   ft(idx{e},:) = ft(idx{e},:) + fte;
end
% Agrego al vector de fuerzas nodales equivalentes las fuerzas
% superficiales calculadas
f = f + ft;
%% se definen los apoyos y sus desplazamientos

idxNODO  = restric(:,1);%% nodos 
dir_desp = restric(:,2);%% dirección

% grados de libertad del desplazamiento conocidos  
% desplazamientos conocidos en los apoyos
ac = restric(:,3);%% desplazamiento

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
   nod_ef = LaG(e, [NL1, NL4, NL5, NL2, NL6, NL7, NL3, NL8, NL9, NL1]);
   plot(xnod(nod_ef, X), xnod(nod_ef, Y), 'b')
   line(xdef(nod_ef, X), xdef(nod_ef, Y), 'Color','r'); % deformada
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
   for i = 1:n_gl
         def{e,i} = B{e,i}*ae;           % calculo las deformaciones
         esf{e,i} = De{mat(e)}*def{e,i}; % calculo los esfuerzos
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


function [esf, error_esf] = extrapolar_esf_def(xnod, LaG, esfuerzo, tipo_esf)
    nno = size(xnod, 1);
    nef = size(LaG, 1);

    num_elem_ady = zeros(nno,1);  % numero de elementos adyacentes
    esf.sum      = zeros(nno,1);
    esf.max      =  -inf(nno,1);
    esf.min      =   inf(nno,1);

    % matriz de extrapolación
    A = [[ 0.12634073, -0.63855959, -0.63855959,  1.87365927,  0.13855959,  0.13855959];
         [-0.63855959, -0.63855959,  0.12634073,  0.13855959,  0.13855959,  1.87365927];
         [-0.63855959,  0.12634073, -0.63855959,  0.13855959,  1.87365927,  0.13855959];
         [-0.20780502,  1.13839679, -0.4627718 ,  0.46755195,  0.17544269, -0.11081461];
         [-0.4627718 ,  1.13839679, -0.20780502, -0.11081461,  0.17544269,  0.46755195];
         [ 1.13839679, -0.4627718 , -0.20780502,  0.17544269, -0.11081461,  0.46755195];
         [ 1.13839679, -0.20780502, -0.4627718 ,  0.17544269,  0.46755195, -0.11081461];
         [-0.4627718 , -0.20780502,  1.13839679, -0.11081461,  0.46755195,  0.17544269];
         [-0.20780502, -0.4627718 ,  1.13839679,  0.46755195, -0.11081461,  0.17544269];
         [ 0.42570639,  0.42570639,  0.42570639, -0.09237306, -0.09237306, -0.09237306]];

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

