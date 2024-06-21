%% Calculo de los desplazamientos verticales y angulos de giro, las 
% reacciones, los momentos flectores y las fuerzas cortantes en una losa de
% Mindlin utilizando los elementos finitos de placa "QL9"

%%
%clear, clc, close all   % borro la memoria, la pantalla y las figuras

%% defino las variables/constantes
X = 1; Y = 2;        % un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; % lectura del codigo
r_ = [ 1, 2, 3, 4, 5, 6, 7, 8,10,11,12];   % GDL a retener en condensación nodal
e_ = [ 13, 14, 15, 16 ,17 ,18,19,20,21,22,23,24];            % GDL a eliminar en condensación nodal


E  = 210000;          % modulo de elasticidad del solido (Pa) = 210GPa
nu = 0.3;            % coeficiente de Poisson
t  = 0.05;           % espesor de la losa (m)
qdistr = -10;     % carga (N/m^2)

% Definimos la geometria de la losa (creada con "generar_malla_losa.m")
%load malla_losa
% ya tenemos en la memoria las variables
% xnod - posicion (x,y) de los nodos
% LaG  - definicion de elementos finitos con respecto a nodos

nef   = size(LaG,1);  % numero de EFs (numero de filas de LaG)
nnoef = size(LaG,2);  % numero de nodos por EF
nno   = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl  = 3*nno;        % numero de grados de libertad (tres por nodo)

gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad

%% Se dibuja la malla de elementos finitos
figure; 
hold on;
for e = 1:nef
   line(xnod(LaG(e,[1:4 1]),X), xnod(LaG(e,[1:4 1]),Y));
   
   % Calculo la posicion del centro de gravedad del elemento finito
   cgx = mean(xnod(LaG(e,:), X));
   cgy = mean(xnod(LaG(e,:), Y));
   text(cgx+0.03, cgy+0.03, num2str(e), 'Color', [1 0 0]);
end
plot(xnod(:,X), xnod(:,Y), 'rx');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
axis([-0.5, 2.5, -0.5, 4.5])
title('Malla de una losa con EFs QL9');

%% Se cargan las funciones de forma junto con sus derivadas
% Se cargan las funciones de forma del elemento lagrangiano de 9 nodos 
% junto con sus derivadas con respecto a xi y a eta
% NOTA estas funciones de forma y sus derivadas se encontraron con el
% programa "codigo/2D/deduccion_funciones_forma/FF_lagrangianos_Q9.m"

%% N son las funciones de forma del elemento lagrangiano de 9 nodos
Nforma = @(xi,eta)[ ((eta - 1)*(xi - 1))/4    % N1
                   -((eta - 1)*(xi + 1))/4    % N2
                    ((eta + 1)*(xi + 1))/4    % N3
                   -((eta + 1)*(xi - 1))/4 ]; % N4
               
Nforma6 = @(xi,eta)[ ((eta - 1)*(xi - 1))/4    % N1
                   -((eta - 1)*(xi + 1))/4     % N2
                    ((eta + 1)*(xi + 1))/4     % N3
                   -((eta + 1)*(xi - 1))/4     % N4
                  ((xi^2 - 1)*(eta - 1))/2     % N2
                 -((eta^2 - 1)*(xi + 1))/2     % N4
                 -((xi^2 - 1)*(eta + 1))/2     % N6
                  ((eta^2 - 1)*(xi - 1))/2 ];  % N8

% derivadas de las funciones de forma con respecto a xi
dN_dxi = @(xi,eta) [  eta/4 - 1/4    % dN1_dxi
                      1/4 - eta/4    %dN2_dxi
                      eta/4 + 1/4    %dN3_dxi
                    - eta/4 - 1/4    %dN4_dxi
                      eta*xi - xi    %dN2_dxi
                    1/2 - eta^2/2    %dN4_dxi
                    -xi*(eta + 1)    % dN6_dxi
                    eta^2/2 - 1/2];  % dN8_dxi


%% Derivadas de N con respecto a eta
dN_deta =  @(xi,eta) [  xi/4 - 1/4   % dN1_deta
                       -xi/4 - 1/4   % dN2_deta
                        xi/4 + 1/4   % dN3_deta
                        1/4 - xi/4   % dN4_deta
                      xi^2/2 - 1/2   % dN2_deta
                        -eta*(xi + 1)% dN4_deta
                        1/2 - xi^2/2 % dN6_deta
                        eta*(xi - 1) ];% dN8_deta


%% parametros de la cuadratura de Gauss-Legendre (INTEGRACION SELECTIVA)
% se asumira aqui el mismo orden de la cuadratura tanto en la direccion de
% xi como en la direccion de eta

% se utilizara integracion COMPLETA
%{
n_gl_b = 3; % orden de la cuadratura de GL para la integracion de Kb
n_gl_s = 3; % orden de la cuadratura de GL para la integracion de Ks
%}

% se utilizara integracion SELECTIVA
n_gl_b = 2; % orden de la cuadratura de GL para la integracion de Kb
n_gl_s = 1; % orden de la cuadratura de GL para la integracion de Ks

% se utilizara integracion REDUCIDA
%{
n_gl_b = 2; % orden de la cuadratura de GL para la integracion de Kb
n_gl_s = 2; % orden de la cuadratura de GL para la integracion de Ks
%}

% calcula las raices (x_gl) y los pesos (w_gl) de polinomios de Legendre
[x_gl_b, w_gl_b]  = gausslegendre_quad(n_gl_b);
[x_gl_s, w_gl_s]  = gausslegendre_quad(n_gl_s);

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
inv_Ksee = cell(nef, 4, 4);
Kbeer    = cell(nef, 4, 24);

inv_Kbee = cell(nef, 4, 4);
Kseer     = cell(nef, 4, 24);

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
   fe  = zeros(3*4,1); 
   %% se calcula la matriz de rigidez de flexion Kb del elemento e 
   Kbe = zeros(3*8);
   det_Je_b = zeros(n_gl_b); % Jacobianos con n_gl_b puntos de integracion   
   for p = 1:n_gl_b
      for q = 1:n_gl_b
         xi_gl  = x_gl_b(p);
         eta_gl = x_gl_b(q);
         [Bb{e,p,q}, det_Je_b(p,q)] = Bb_RM(xi_gl, eta_gl, xe, ye, dN_dxi, dN_deta);
         % se arma la matriz de rigidez del elemento e
         Kbe = Kbe + Bb{e,p,q}'*Dbg*Bb{e,p,q}*det_Je_b(p,q)*w_gl_b(p)*w_gl_b(q);
      end
   end
   
   %% se calcula la matrix Ks
   Kse = zeros(3*8);   
   det_Je_s = zeros(n_gl_s); % Jacobianos con n_gl_s puntos de integracion
   for p = 1:n_gl_s
      for q = 1:n_gl_s
         xi_gl  = x_gl_s(p);        
         eta_gl = x_gl_s(q);
         [Bs{e,p,q}, det_Je_s(p,q)] = Bs_RM(xi_gl, eta_gl, xe, ye, Nforma6, dN_dxi, dN_deta);   

         % se arma la matriz de rigidez del elemento e
         Kse = Kse + Bs{e,p,q}'*Dsg*Bs{e,p,q}*det_Je_s(p,q)*w_gl_s(p)*w_gl_s(q);         
      end
   end 
   
   %% se calcula la matriz NN
   Mbe = zeros(3*nnoef); % matriz que se utiliza en el calculo de fe   
   for p = 1:n_gl_b
      for q = 1:n_gl_b
         xi_gl  = x_gl_b(p);
         eta_gl = x_gl_b(q);
         % Se evaluan las funciones de forma en los puntos de integracion
         % de Gauss-Legendre
         N = Nforma(xi_gl, eta_gl);
         
         % Se ensambla la matriz de funciones de forma N
         NN{e,p,q} = zeros(3,3*nnoef);
         for i = 1:nnoef            
            NN{e,p,q}(:,3*i-2:3*i) = diag([N(i) N(i) N(i)]);
         end
         % matriz requerida para calcular el vector de fuerzas nodales 
         % equivalentes (se utiliza la integracion completa)
         %% vector de fuerzas nodales equivalentes        
         if (xe(1) >= 0 && xe(2) <= 2) && ...
           (ye(2) >= 0 && ye(3) <= 4)
            fe = fe + NN{e,p,q}'*[qdistr 0 0]'*det_Je_b(p,q)*w_gl_b(p)*w_gl_b(q);
         end                                                                                       % REVISAR !!!!!!!!!!!!!!!!   
      end
   end  
   % se condensan los GDL jerárquicos
    Kberr = Kbe(r_,r_); 
    Kbeer{e}= Kbe(e_,r_);
    Kbere = Kbe(r_,e_);
    inv_Kbee{e} = (Kbe(e_,e_))^(-1);
    
    Kebec = Kberr - Kbere * inv_Kbee{e} * Kbeer{e};

    
    % se condensan los GDL jerárquicos
    Kserr = Kse(r_,r_); 
    Kseer{e}= Kse(e_,r_);
    Ksere = Kse(r_,e_);
    inv_Ksee{e} = (Kse(e_,e_))^(-1);
    
    Kesec = Kserr - Ksere * inv_Ksee{e} * Kseer{e};

   
   %% se asocian los grados de libertad del elemento locales a los globales
   idx{e} = [ gdl(LaG(e,1),:)  gdl(LaG(e,2),:)  gdl(LaG(e,3),:)  ...
              gdl(LaG(e,4),:)  ];

   %% se procede al ensamblaje
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Kebec + Kesec;
   f(idx{e})        = f(idx{e}) + fe;
end

%% Muestro la configuracion de la matriz K (K es rala)
figure
spy(K);
title('Los puntos representan los elementos diferentes de cero')

%% grados de libertad del desplazamiento conocidos y desconocidos

% determino los grados de libertad correspondientes a los bordes
lado_x0 = find(abs(xnod(:,X) - 0) < 1e-4);     
lado_y0 = find(abs(xnod(:,Y) - 0) < 1e-4);
lado_x2 = find(abs(xnod(:,X) - 2) < 1e-4);     
lado_y4 = find(abs(xnod(:,Y) - 4) < 1e-4);

c = [ gdl(lado_x0,ww); gdl(lado_x0,ty); 
      gdl(lado_x2,ww); gdl(lado_x2,ty);
      gdl(lado_y0,ww); gdl(lado_y0,tx);
      gdl(lado_y4,ww); gdl(lado_y4,tx) ];

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
escala = 10;            % factor de escalamiento de la deformada
xdef   = escala*vect_mov; % posicion de la deformada
figure; 
hold on; 
grid on;
colorbar
for e = 1:nef
   fill3(xnod(LaG(e,[1:4 1]),X), ...
         xnod(LaG(e,[1:4 1]),Y), ...
         xdef(LaG(e,[1:4 1]),ww),...
         xdef(LaG(e,[1:4 1]),ww)); %deformada
end
daspect([1 1 1]); % similar a "axis equal", pero en 3D
axis tight
%colorbar('YTick',-0.6:0.05:0)
title(sprintf('Deformada escalada %d veces', escala), 'FontSize', 20);
colormap jet
view(3);

%% Se dibuja de la malla de elementos finitos y las deformaciones de esta
%figure; 
%hold on; 
%grid on;
%colorbar
%for e = 1:nef
%   dibujar_EF_Q89_RM(xnod(LaG(e,:),X), xnod(LaG(e,:),Y), ...
%      Nforma, a(idx{e}), t, escala, escala);
%end
%daspect([1 1 1]); % similar a "axis equal", pero en 3D
%axis tight
%title(sprintf('Deformada escalada %d veces', escala), 'FontSize', 20);
%colormap jet
%view(3);

%% En los puntos de integracion de Gauss-Legendre calcular:
%% El vector de momentos flectores y torsores (2x2)
%% El vector de fuerzas cortantes (1x1 o 2x2)
n_gl_b = 2; [x_gl_b, w_gl_b]  = gausslegendre_quad(n_gl_b);

% Observe que n_gl_s = 1; interpola mal la fuerza cortante.
n_gl_s = 2; [x_gl_s, w_gl_s]  = gausslegendre_quad(n_gl_s);

%% se calcula de nuevo Bb y Bs en cada punto de GL
Bb = cell(nef,n_gl_b,n_gl_b); % matrices de deformacion generalizada de flexion
Bs = cell(nef,n_gl_s,n_gl_s); % matrices de deformacion generalizada de cortante
for e = 1:nef      % ciclo sobre todos los elementos finitos
    xe = xnod(LaG(e,:),X);
    ye = xnod(LaG(e,:),Y);
    
    %% se calcula la matrix Bb en los puntos de integracion de GL para el
    % calculo de los momentos flectores y torsores
    for p = 1:n_gl_b
        for q = 1:n_gl_b
            xi_gl  = x_gl_b(p);
            eta_gl = x_gl_b(q);
            Bb{e,p,q} = Bb_RM(xi_gl, eta_gl, xe, ye, dN_dxi, dN_deta);
        end
    end
    
    %% se calcula la matrix Bs en los puntos de integracion de GL para el
    % calculo de las fuerzas cortantes
    for p = 1:n_gl_s
        for q = 1:n_gl_s
            xi_gl  = x_gl_s(p);
            eta_gl = x_gl_s(q);
            Bs{e,p,q} = Bs_RM(xi_gl, eta_gl, xe, ye, Nforma, dN_dxi, dN_deta);
        end
    end
end

%% Se calculan los momentos y las fuerzas en los puntos de GL
sigmag_b = cell(nef, n_gl_b, n_gl_b); % momentos flectores y torsores
sigmag_s = cell(nef, n_gl_s, n_gl_s); % fuerzas cortantes
for e = 1:nef               % ciclo sobre todos los elementos finitos
   for p = 1:n_gl_b
      for q = 1:n_gl_b
         sigmag_b{e,p,q} = Dbg*Bb{e,p,q}*a(idx{e});
      end
   end
   
   for p = 1:n_gl_s
      for q = 1:n_gl_s
         sigmag_s{e,p,q} = Dsg*Bs{e,p,q}*a(idx{e});
      end
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
A = [[  3^(1/2)/2 + 1,            -1/2,            -1/2,   1 - 3^(1/2)/2];
    [            -1/2,   1 - 3^(1/2)/2,   3^(1/2)/2 + 1,            -1/2];
    [   1 - 3^(1/2)/2,            -1/2,            -1/2,   3^(1/2)/2 + 1];
    [            -1/2,   3^(1/2)/2 + 1,   1 - 3^(1/2)/2,            -1/2]];

for e = 1:nef

   Mx(LaG(e,:),:)  = Mx(LaG(e,:),:)  + A * [ sigmag_b{e,1,1}(1)
                                             sigmag_b{e,1,2}(1)
                                             sigmag_b{e,2,1}(1)
                                             sigmag_b{e,2,2}(1) ];

   My(LaG(e,:),:)  = My(LaG(e,:),:)  + A * [ sigmag_b{e,1,1}(2)
                                             sigmag_b{e,1,2}(2)
                                             sigmag_b{e,2,1}(2)
                                             sigmag_b{e,2,2}(2) ];

   Mxy(LaG(e,:),:) = Mxy(LaG(e,:),:) + A * [ sigmag_b{e,1,1}(3)
                                             sigmag_b{e,1,2}(3)
                                             sigmag_b{e,2,1}(3)
                                             sigmag_b{e,2,2}(3) ];

   switch n_gl_s
     case 1
       Qx(LaG(e,:),:) = Qx(LaG(e,:),:) + sigmag_s{e}(1);
       Qy(LaG(e,:),:) = Qy(LaG(e,:),:) + sigmag_s{e}(2);
     case 2
       Qx(LaG(e,:),:)  = Qx(LaG(e,:),:)  + A * [ sigmag_s{e,1,1}(1)
                                                 sigmag_s{e,1,2}(1)
                                                 sigmag_s{e,2,1}(1)
                                                 sigmag_s{e,2,2}(1) ];

       Qy(LaG(e,:),:)  = Qy(LaG(e,:),:)  + A * [ sigmag_s{e,1,1}(2)
                                                 sigmag_s{e,1,2}(2)
                                                 sigmag_s{e,2,1}(2)
                                                 sigmag_s{e,2,2}(2) ];
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
       fill(xnod(LaG(e,1:4),X), xnod(LaG(e,1:4),Y), variable(LaG(e,1:4)));
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