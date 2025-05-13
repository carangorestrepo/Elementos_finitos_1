%% Calculo de los desplazamientos verticales y angulos de giro, las 
% reacciones, los momentos flectores y las fuerzas cortantes en una losa de
% Mindlin utilizando los elementos finitos de placa "TL6"

%%
%clear, clc, close all   % borro la memoria, la pantalla y las figuras

%% defino las variables/constantes
X = 1; Y = 2;        % un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; % lectura del codigo
%E  = 210000;          % modulo de elasticidad del solido (Pa) = 210GPa
%nu = 0.3;            % coeficiente de Poisson
%t  = 0.05;           % espesor de la losa (m)
%qdistr = -10;     % carga (N/m^2)

% Definimos la geometria de la losa (creada con "generar_malla_losa.m")
%load malla_losa
% ya tenemos en la memoria las variables
% xnod - posicion (x,y) de los nodos
% LaG  - definicion de elementos finitos con respecto a nodos

nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl = 3*nno;        % numero de grados de libertad (tres por nodo)
gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad
nef  = size(LaG,1);  % numero de EFs (numero de filas de LaG)
nnoef = size(LaG,2);  % numero de nodos por EF

%% Se dibuja la malla de elementos finitos
figure;
hold on;
cg = zeros(nef,2); % almacena el centro de gravedad de los EFs
for e = 1:nef
   % se dibuja el EF e
   nod_ef = LaG(e,[1,4,2,5,3,6,1]);
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
%% Se cargan las funciones de forma junto con sus derivadas
% Se cargan las funciones de forma del elemento lagrangiano de 9 nodos 
% junto con sus derivadas con respecto a xi y a eta
% NOTA estas funciones de forma y sus derivadas se encontraron con el
% programa "codigo/2D/deduccion_funciones_forma/FF_lagrangianos_Q9.m"

%% N son las funciones de forma del elemento lagrangiano de 9 nodos
 Nforma =  @(xi, eta) [...
 (eta + xi - 1)*(2*eta + 2*xi - 1)
                     xi*(2*xi - 1)
                   eta*(2*eta - 1)
            -xi*(4*eta + 4*xi - 4)
                          4*eta*xi
           -eta*(4*eta + 4*xi - 4) 
                             ];
 
%% Derivadas de N con respecto a xi
dN_dxi = @(xi,eta) [ ...
 4*eta + 4*xi - 3
         4*xi - 1
                0
 4 - 8*xi - 4*eta
            4*eta
           -4*eta                     
];  % dN6_dxi   

 
%% Derivadas de N con respecto a eta
%% Derivadas de N con respecto a eta
dN_deta = @(xi,eta) [ ...
 4*eta + 4*xi - 3
                0
        4*eta - 1
            -4*xi
             4*xi
 4 - 4*xi - 8*eta                        
];  % dN6_deta

%% parametros de la cuadratura de Gauss-Legendre (INTEGRACION SELECTIVA)
% se asumira aqui el mismo orden de la cuadratura tanto en la direccion de
% xi como en la direccion de eta

% se utilizara integracion COMPLETA
%{
n_gl_b = 3; % orden de la cuadratura de GL para la integracion de Kb
n_gl_s = 3; % orden de la cuadratura de GL para la integracion de Ks
%}

% se utilizara integracion SELECTIVA
n_gl_b = 3; % orden de la cuadratura de GL para la integracion de Kb
n_gl_s = 2; % orden de la cuadratura de GL para la integracion de Ks

% se utilizara integracion REDUCIDA
%{
n_gl_b = 2; % orden de la cuadratura de GL para la integracion de Kb
n_gl_s = 2; % orden de la cuadratura de GL para la integracion de Ks
%}

% calcula las raices (x_gl) y los pesos (w_gl) de polinomios de Legendre

xw=TriGaussPoints(n_gl_b);

x_gl_b = xw(:,1);
e_gl_b = xw(:,2);
w_gl_b =  xw(:,3);
n_gl_b = size(x_gl_b,1);  %# N�mero de puntos de Gauss.

xw=TriGaussPoints(n_gl_s);

x_gl_s = xw(:,1);
e_gl_s = xw(:,2);
w_gl_s =  xw(:,3);
n_gl_s = size(x_gl_s,1);  %# N�mero de puntos de Gauss.

%[x_gl_b, w_gl_b]  = gausslegendre_quad(n_gl_b);
%[x_gl_s, w_gl_s]  = gausslegendre_quad(n_gl_s);

%% matrices constitutivas del elemento
Db = (E/(1-nu^2))* [ 1    nu   0
                     nu   1    0
                     0    0    (1-nu)/2 ];              
G = E/(2*(1+nu));  % modulo de rigidez
alpha = 5/6;       % coeficiente de distorsion transversal de la losa de RM
Ds = diag([alpha*G, alpha*G]);
               
Dbg = (t^3/12)*Db; % matriz constitutiva generalizada de flexion
Dsg = t*Ds;        % matriz constitutiva generalizada de cortante

%% se reserva la memoria RAM de diferentes variables
K   = sparse(ngdl,ngdl); % matriz de rigidez global como RALA (sparse)
f   = zeros(ngdl,1);     % vector de fuerzas nodales equivalentes global
idx = cell(nef, 1);      % grados de libertad de cada elemento finito

% en los siguientes contenedores se almacenara la matriz respectiva para 
% cada punto de integracion: 
NN = cell(nef,n_gl_b,n_gl_b); % matrices de funciones de forma calculadas con n_gl_b puntos de integracion
Bb = cell(nef,n_gl_b,n_gl_b); % matrices de deformacion generalizada de flexion
Bs = cell(nef,n_gl_s,n_gl_s); % matrices de deformacion generalizada de cortante

%% se ensambla la matriz de rigidez global K y el vector de fuerzas nodales
%% equivalentes global f
for e = 1:nef      % ciclo sobre todos los elementos finitos
   xe = xnod(LaG(e,:),X);   
   ye = xnod(LaG(e,:),Y);    
    
   %% se calcula la matriz de rigidez de flexion Kb del elemento e 
   Kbe = zeros(3*nnoef);
   det_Je_b = zeros(n_gl_b,1); % Jacobianos con n_gl_b puntos de integracion   
   for p = 1:n_gl_b
      %for q = 1:n_gl_b
         xi_gl  = x_gl_b(p);
         eta_gl = e_gl_b(p);
         
         %% Se evaluan las derivadas de las funciones de forma en los puntos
        %% de integracion de Gauss-Legendre
        ddN_dxi  = dN_dxi (xi_gl, eta_gl);
        ddN_deta = dN_deta(xi_gl, eta_gl);

        %% Se utilizan las funciones de forma de w para el calculo de la 
        %% transformacion isoparametrica
        dx_dxi  = sum(ddN_dxi  .* xe);   dy_dxi  = sum(ddN_dxi  .* ye);
        dx_deta = sum(ddN_deta .* xe);   dy_deta = sum(ddN_deta .* ye);

        %% Se ensambla la matriz Jacobiana del elemento
        Je = [ dx_dxi   dy_dxi
               dx_deta  dy_deta ];

        %% Se calcula el determinante del Jacobiano
        det_Je_b(p) = det(Je);
        if det_Je_b(p) <= 0
           error('El det_Je es negativo');
        end

        %% Se ensambla la matriz de deformacion del elemento Bb
        nno_b = length(xe);
        Bbi = zeros(3,3*nno_b);
        for i = 1:nno_b   
           dNi_dx = (+dy_deta*ddN_dxi(i) - dy_dxi*ddN_deta(i))/det_Je_b(p);
           dNi_dy = (-dx_deta*ddN_dxi(i) + dx_dxi*ddN_deta(i))/det_Je_b(p);

           Bbi(:,(3*i-2):(3*i)) = [ 0 -dNi_dx       0    % se ensambla y
                                   0       0 -dNi_dy    % asigna la matriz
                                   0 -dNi_dy -dNi_dx ]; % Bb_i
        end
         %[Bb{e,p}, det_Je_b(p)] = Bb_RM(xi_gl, eta_gl, xe, ye, dN_dxi, dN_deta);
         % se arma la matriz de rigidez del elemento e
         Bb{e,p}=Bbi;
         Kbe = Kbe + Bb{e,p}'*Dbg*Bb{e,p}*det_Je_b(p)*w_gl_b(p);
      %end
   end
   
   %% se calcula la matrix Ks
   Kse = zeros(3*nnoef);   
   det_Je_s = zeros(n_gl_s); % Jacobianos con n_gl_s puntos de integracion
   for p = 1:n_gl_s
      %for q = 1:n_gl_s
         xi_gl  = x_gl_s(p);        
         eta_gl = e_gl_s(p);
         
         %% Se evaluan las funciones de forma en los puntos de integracion
        %% de Gauss-Legendre
        NNforma = Nforma(xi_gl, eta_gl);

        %% Se evaluan las derivadas de las funciones de forma en los puntos
        %% de integracion de Gauss-Legendre
        ddN_dxi  = dN_dxi (xi_gl, eta_gl);
        ddN_deta = dN_deta(xi_gl, eta_gl);

        %% Se utilizan las funciones de forma de w para el calculo de la 
        %% transformacion isoparametrica
        dx_dxi  = sum(ddN_dxi  .* xe);   dy_dxi  = sum(ddN_dxi  .* ye);
        dx_deta = sum(ddN_deta .* xe);   dy_deta = sum(ddN_deta .* ye);

        %% Se ensambla la matriz Jacobiana del elemento
        Je = [ dx_dxi   dy_dxi
               dx_deta  dy_deta ];

        %% Se calcula el determinante del Jacobiano
        det_Je_s(p) = det(Je);
        if det_Je_s(p) <= 0
           error('El det_Je es negativo');
        end

        %% Se ensambla la matriz de deformacion del elemento Bs
        nno_s = length(xe);
        Bsi  = zeros(2,3*nno_s);
        for i = 1:nno_s
           dNi_dx = (+dy_deta*ddN_dxi(i) - dy_dxi*ddN_deta(i))/det_Je_s(p);
           dNi_dy = (-dx_deta*ddN_dxi(i) + dx_dxi*ddN_deta(i))/det_Je_s(p);
           % se ensambla y asigna la matriz Bs_i
           Bsi(:,(3*i-2):(3*i)) = [ dNi_dx  -NNforma(i)  0     
                                   dNi_dy  0            -NNforma(i) ];                        
        end    
        Bs{e,p}=Bsi;
         %[Bs{e,p}, det_Je_s(p)] = Bs_RM(xi_gl, eta_gl, xe, ye, Nforma, dN_dxi, dN_deta);   
         % se arma la matriz de rigidez del elemento e
         Kse = Kse + Bs{e,p}'*Dsg*Bs{e,p}*det_Je_s(p)*w_gl_s(p);         
      %end
   end 
   
   %% se calcula la matriz NN
   Mbe = zeros(3*nnoef); % matriz que se utiliza en el calculo de fe   
   for p = 1:n_gl_b
      %for q = 1:n_gl_b
         xi_gl  = x_gl_b(p);
         eta_gl = e_gl_b(p);
         % Se evaluan las funciones de forma en los puntos de integracion
         % de Gauss-Legendre
         N = Nforma(xi_gl, eta_gl);
         
         % Se ensambla la matriz de funciones de forma N
         NN{e,p} = zeros(3,3*nnoef);
         for i = 1:nnoef            
            NN{e,p}(:,3*i-2:3*i) = diag([N(i) N(i) N(i)]);
         end
   
         % matriz requerida para calcular el vector de fuerzas nodales 
         % equivalentes (se utiliza la integracion completa)
         Mbe = Mbe + NN{e,p}'*NN{e,p}*det_Je_b(p)*w_gl_b(p);                                              % REVISAR !!!!!!!!!!!!!!!!   
      %end
   end  
   
   %% se calcula el vector de fuerzas nodales equivalentes del elemento e      
   %xa = xnod(LaG(e,1),X);   ya = xnod(LaG(e,1),Y);
   %xb = xnod(LaG(e,5),X);   yb = xnod(LaG(e,5),Y);
   %if (xa >= 0 && xb <= 2) && (ya >= 0 && yb <= 4)
      ffe = zeros(nnoef, 3); ffe(:,ww) = qa;
      ffe = reshape(ffe', 3*nnoef,1);                                                                                             % REVISAR !!!!!!!!!!!!!
   %else
    %  ffe = zeros(3*nnoef,1);
   %end  
   fe = Mbe*ffe;   
   
   %% se asocian los grados de libertad del elemento locales a los globales
   idx{e} = [ gdl(LaG(e,1),:)  gdl(LaG(e,2),:)  gdl(LaG(e,3),:)  ...
              gdl(LaG(e,4),:)  gdl(LaG(e,5),:)  gdl(LaG(e,6),:)];

   %% se procede al ensamblaje
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Kbe + Kse;
   f(idx{e})        = f(idx{e}) + fe;
end

%% Muestro la configuracion de la matriz K (K es rala)
figure
spy(K);
title('Los puntos representan los elementos diferentes de cero')

%% grados de libertad del desplazamiento conocidos y desconocidos

% determino los grados de libertad correspondientes a los bordes
%lado_x0 = find(abs(xnod(:,X) - 0) < 1e-4);     
%lado_y0 = find(abs(xnod(:,Y) - 0) < 1e-4);
%lado_x2 = find(abs(xnod(:,X) - 2) < 1e-4);     
%lado_y4 = find(abs(xnod(:,Y) - 4) < 1e-4);

%c = [ gdl(lado_x0,ww); gdl(lado_x0,ty); 
%      gdl(lado_x2,ww); gdl(lado_x2,ty);
%      gdl(lado_y0,ww); gdl(lado_y0,tx);
%      gdl(lado_y4,ww); gdl(lado_y4,tx) ];

d = setdiff(1:ngdl,c)';

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd   |   | Kcc Kcd || ac |   | fd |  Recuerde que qc = 0
%|      | = |         ||    | - |    |
%| qc=0 |   | Kdc Kdd || ad |   | fc | 

%% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = zeros(length(c),1); % desplazamientos conocidos en contorno

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);        % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a = nan(ngdl,1);     a(c) = ac;   a(d) = ad; % desplazamientos
q = zeros(ngdl,1);   q(c) = qd;              % fuerzas nodales equivalentes

%% imprimo los resultados
vect_mov = reshape(a,3,nno)'; % vector de movimientos
%{
format short g
disp('Desplazamientos nodales                      ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
for i = 1:nno
   fprintf('Nodo %3d: w = %12.4g m, tx = %12.4g rad, ty = %12.4g rad\n', ...
      i, vect_mov(i,ww), vect_mov(i,tx), vect_mov(i,ty));
end

disp(' ');
disp('Fuerzas nodales de equilibrio (solo se imprimen los diferentes de cero)');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
q = reshape(q,3,nno)';
for i = 1:nno   
   if ~isequal(q(i,:),[0 0 0])
      fprintf('Nodo %3d W = %12.4g N, Mx = %12.4g N-m, My = %12.4g N-m\n', ...
         i, q(i,ww), q(i,tx), q(i,ty));
   end
end
%}

%% Se dibuja el plano medio de la malla de elementos finitos y las deformaciones de esta
escala = 1;            % factor de escalamiento de la deformada
xdef   = escala*vect_mov; % posicion de la deformada
figure; 
hold on; 
grid on;
colorbar
for e = 1:nef
   fill3(xnod(LaG(e,[1,4,2,5,3,6,1]),X), ...
         xnod(LaG(e,[1,4,2,5,3,6,1]),Y), ...
         xdef(LaG(e,[1,4,2,5,3,6,1]),ww),...
         xdef(LaG(e,[1,4,2,5,3,6,1]),ww)); %deformada
end
daspect([1 1 1]); % similar a "axis equal", pero en 3D
axis tight
%colorbar('YTick',-0.6:0.05:0)
title(sprintf('Deformada escalada %d veces', escala), 'FontSize', 20);
colormap jet
view(3);


%% se calcula de nuevo Bb y Bs en cada punto de GL
Bb = cell(nef,n_gl_b,n_gl_b); % matrices de deformacion generalizada de flexion
Bs = cell(nef,n_gl_s,n_gl_s); % matrices de deformacion generalizada de cortante
for e = 1:nef      % ciclo sobre todos los elementos finitos
    xe = xnod(LaG(e,:),X);
    ye = xnod(LaG(e,:),Y);
    
    %% se calcula la matrix Bb en los puntos de integracion de GL para el
    % calculo de los momentos flectores y torsores
    for p = 1:n_gl_b
        %for q = 1:n_gl_b
            xi_gl  = x_gl_b(p);
            eta_gl = e_gl_b(p);
            Bb{e,p} = Bb_RM(xi_gl, eta_gl, xe, ye, dN_dxi, dN_deta);
        %end
    end
    
    %% se calcula la matrix Bs en los puntos de integracion de GL para el
    % calculo de las fuerzas cortantes
    for p = 1:n_gl_s
        %for q = 1:n_gl_s
            xi_gl  = x_gl_s(p);
            eta_gl = e_gl_s(p);
            Bs{e,p} = Bs_RM(xi_gl, eta_gl, xe, ye, Nforma, dN_dxi, dN_deta);
        %end
    end
end

%% Se calculan los momentos y las fuerzas en los puntos de GL
sigmag_b = cell(nef, n_gl_b, n_gl_b); % momentos flectores y torsores
sigmag_s = cell(nef, n_gl_s, n_gl_s); % fuerzas cortantes
for e = 1:nef               % ciclo sobre todos los elementos finitos
   for p = 1:n_gl_b
      %for q = 1:n_gl_b
         sigmag_b{e,p} = Dbg*Bb{e,p}*a(idx{e});
      %end
   end
   
   for p = 1:n_gl_s
      %for q = 1:n_gl_s
         sigmag_s{e,p} = Dsg*Bs{e,p}*a(idx{e});
      %end
   end
end

%% Se extrapolan los momentos flectores y fuerzas cortantes a los nodos
%% Se extrapolan los esfuerzos y las deformaciones a los nodos
num_elem_ady = zeros(nno,1);  % numero de elementos adyacentes
Mx  = zeros(nno,1);  Qx = zeros(nno,1);
My  = zeros(nno,1);  Qy = zeros(nno,1);
Mxy = zeros(nno,1);

% matriz de extrapolacion de esfuerzos para un elemento lagrangiano de 9
% nodos
    % matriz de extrapolación
Ab=[[    -9,    5/2,    5/2,      5];
    [    -9,    5/2,      5,    5/2];
    [   9/4,    5/4,   -5/4,   -5/4];
    [ 81/16, -35/16, -15/16, -15/16];
    [ -27/8,   15/8,   15/8,    5/8];
    [ -27/8,   15/8,    5/8,   15/8]];

As=[[ -1/3, -1/3,  5/3];
    [ -1/3,  5/3, -1/3];
    [  5/3, -1/3, -1/3]];
 

     
for e = 1:nef
   sigmag_bb=[sigmag_b{e,:}]; 
   Mx(LaG(e,:),:)  = Mx(LaG(e,:),:)  + Ab * sigmag_bb(1,:)';

   My(LaG(e,:),:)  = My(LaG(e,:),:)  + Ab *  sigmag_bb(2,:)';

   Mxy(LaG(e,:),:) = Mxy(LaG(e,:),:) + Ab *  sigmag_bb(3,:)';

   switch n_gl_s
     case 1
       sigmag_ss=[sigmag_s{e,:,:}];  
       Qx(LaG(e,:),:) = Qx(LaG(e,:),:) +  As * sigmag_ss(1,:)';
       Qy(LaG(e,:),:) = Qy(LaG(e,:),:) +  As * sigmag_ss(2,:)';
     case 2
       sigmag_ss=[sigmag_s{e,:,:}];
       Qx(LaG(e,:),:)  = Qx(LaG(e,:),:)  + As * sigmag_ss(1,:)';

       Qy(LaG(e,:),:)  = Qy(LaG(e,:),:)  + As * sigmag_ss(2,:)';
   end                                         

   num_elem_ady(LaG(e,:),:) = num_elem_ady(LaG(e,:),:) + 1;
end

%% Alisado (promedio de los esfuerzos en los nodos)
Mx  =  Mx./num_elem_ady;
My  =  My./num_elem_ady;
Mxy = Mxy./num_elem_ady;
Qx  =  Qx./num_elem_ady;  
Qy  =  Qy./num_elem_ady;

%% Se grafican los momentos
figure
subplot(1,3,1); plot_M_or_Q(nef, xnod, LaG, Mx,  'Momentos Mx (N-m/m)');
subplot(1,3,2); plot_M_or_Q(nef, xnod, LaG, My,  'Momentos My (N-m/m)');
subplot(1,3,3); plot_M_or_Q(nef, xnod, LaG, Mxy, 'Momentos Mxy (N-m/m)');

%% Se grafican los cortantes
figure
subplot(1,2,1); plot_M_or_Q(nef, xnod, LaG, Qx,  'Cortantes Qx (N/m)');
subplot(1,2,2); plot_M_or_Q(nef, xnod, LaG, Qy,  'Cortantes Qy (N/m)');

%% Se calculan y grafican para cada elemento los momentos principales y
%% sus direcciones
Mt_max = sqrt(((Mx-My)/2).^2 + Mxy.^2); % momento torsion maximo
Mf1_xy = (Mx+My)/2 + Mt_max;            % momento flector maximo
Mf2_xy = (Mx+My)/2 - Mt_max;            % momento flector minimo
ang  = 0.5*atan2(2*Mxy, Mx-My);         % angulo de inclinacion de Mf1_xy

%% Mf1_xy, Mf2_xy, Mt_max
figure
subplot(1,3,1); plot_M_or_Q(nef, xnod, LaG, Mf1_xy, 'Mf1_{xy} (N-m/m)', { ang })
subplot(1,3,2); plot_M_or_Q(nef, xnod, LaG, Mf2_xy, 'Mf2_{xy} (N-m/m)', { ang+pi/2 })
subplot(1,3,3); plot_M_or_Q(nef, xnod, LaG, Mt_max, 'Mt_{max} (N-m/m)', { ang+pi/4, ang-pi/4 })

%% Se calculan y grafican los cortantes maximos, junto con su angulo de inclinacion
Q_max = hypot(Qx, Qy);
ang   = atan2(Qy, Qx);

figure
plot_M_or_Q(nef, xnod, LaG, Q_max, 'Q_{max} (N/m)', { ang })

%%
return; % bye, bye!

%%
function plot_M_or_Q(nef, xnod, LaG, variable, texto, angulos)
    X = 1; Y = 2;
    hold on; 
    colorbar;
    % Por simplicidad no se graficaran los resultados asociados al nodo 9
    for e = 1:nef  
       fill(xnod(LaG(e,[1,4,2,5,3,6,1]),X), xnod(LaG(e,[1,4,2,5,3,6,1]),Y), variable(LaG(e,[1,4,2,5,3,6,1])));
    end
    axis equal tight
    colormap jet
    title(texto, 'FontSize',20);
   
    esc = 0.5;
    if nargin == 6
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
end