%% 
% Calculo de los desplazamientos en una placa utilizando la teoria de
% Reissner-Mindlin y el elemento finito de placa DKQ/DKMQ 
%
% Algoritmo documentado en:
% Katili, I. (1993), A new discrete Kirchhoff-Mindlin element based on 
% Mindlin-Reissner plate theory and assumed shear strain fields-part II: 
% An extended DKQ element for thick-plate bending analysis. Int. J. Numer. 
% Meth. Engng., 36: 1885-1908. https://doi.org/10.1002/nme.1620361107
%
% Este es el algoritmo de losas usado en MIDAS y AUTODESK ROBOT.
% Se programo intentando seguir la nomenclatura del articulo
%
% Por:
% Diego Andres Alvarez Marin (daalvarez@unal.edu.co)
% Sebastian Jaramillo Moreno

%% borro la memoria, la pantalla y las figuras
%clear, clc, %close all 

%% seleccion del tipo de EF de losa a emplear
DKMQ = 1;
DKQ  = 2;
EFtype = DKMQ;

%% defino las variables/constantes
X = 1; Y = 2; Z = 3; % un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; % lectura del codigo


% Definimos la geometria de la losa
%losa
% ya tenemos en la memoria las variables
% xnod - posicion (x,y) de los nodos
% LaG  - definicion de elementos finitos con respecto a nodos
nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl = 3*nno;        % numero de grados de libertad (tres por nodo)
gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad
nef  = size(LaG,1);  % numero de EFs (numero de filas de LaG)

%% Se dibuja la malla de elementos finitos
figure;
hold on;
cgx = zeros(1,nef); cgy = zeros(1,nef); % almacena el centro de gravedad
for e = 1:nef
   line(xnod(LaG(e,[1 2 3 1]),X), xnod(LaG(e,[1 2 3 1]),Y));
   
   % Calculo la posicion del centro de gravedad del EF
   cgx(e) = mean(xnod(LaG(e,:),X));
   cgy(e) = mean(xnod(LaG(e,:),Y));
   text(cgx(e), cgy(e), num2str(e), 'Color', [1 0 0]);
end
plot(xnod(:,X), xnod(:,Y), 'r*');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
%axis equal tight
title('Malla de elementos finitos');

%% Parametros de la cuadratura de Gauss-Legendre
% se asumira aqui el mismo orden de la cuadratura tanto en la direccion de
% xi como en la direccion de eta
n         = 4; % orden de la cuadratura de Gauss-Legendre
%[x_gl, w_gl] = gausslegendre_quad(n_gl);
xw=TriGaussPoints(n);

x_gl = xw(:,1);
e_gl = xw(:,2);
w_gl =  xw(:,3);
n_gl = size(e_gl,1);  %# Número de puntos de Gauss.
  

%% matrices constitutivas
Db = (E*h^3/(12*(1-nu^2)));   % plate rigidity
Hb = Db * [ 1  nu 0           % matriz constitutiva de flexion generalizada
            nu 1  0           % (Dbe en la nomenclatura del curso) 
            0  0  (1-nu)/2 ]; 
G  = E/(2*(1+nu));     % modulo de cortante
Hs = (5/6)*G*h*eye(2); % matriz constitutiva de cortante generalizada (Dse)

T = rho * diag([h, 0, 0]);
%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
Kplaca = sparse(ngdl,ngdl); % rigidez de placa: flexion + cortante
Ksuelo = sparse(ngdl,ngdl); % rigidez Winkler
K   = sparse(ngdl,ngdl);    % matriz total
M   = sparse(ngdl,ngdl);    % matriz de masa global como RALA (sparse)
f   = zeros(ngdl,1);        % vector de fuerzas nodales equivalentes global
N   = cell(nef, n_gl);
Bb  = cell(nef, n_gl);
Bs  = cell(nef, n_gl);
idx = cell(nef, 1);         % grados de libertad de cada elemento finiton
for e = 1:nef               % ciclo sobre todos los elementos finitos
    %% Longitudes de los lados, cosenos y senos (Figura 4)
    xe = xnod(LaG(e,:),X);       ye = xnod(LaG(e,:),Y);
    x21 = xe(2) - xe(1);         y21 = ye(2) - ye(1); 
    x32 = xe(3) - xe(2);         y32 = ye(3) - ye(2);
    x13 = xe(1) - xe(3);         y13 = ye(1) - ye(3);    

    xji = [ x21 x32 x13  ];   yji = [ y21 y32 y13 ];   
    
    Lk = hypot(xji, yji);      Ck =xji./Lk;      Sk = yji./Lk;
    L4=Lk(1);
    L5=Lk(2);
    L6=Lk(3);
    
    C4=Ck(1);
    C5=Ck(2);
    C6=Ck(3);
    S4=Sk(1); 
    S5=Sk(2);
    S6=Sk(3);
    
    %% Ciclo sobre los puntos de Gauss para calcular Kbe, Kse y fe
    Kbe = zeros(9);
    Me = zeros(9);
    He = zeros(9);
    Kse = zeros(9);
    fe  = zeros(9,1);
    det_Je = zeros(n_gl,1); % almacenara los Jacobianos
    
    for pp = 1:n_gl
       % for qq = 1:n_gl           
            %% Se evaluan las funciones de forma y sus derivadas en los 
            % puntos de Gauss
            % y su determinante (el Jacobiano)
            %det_Je(pp) = det(Je);
            xe = xnod(LaG(e,1:3),X);       ye = xnod(LaG(e,1:3),Y);
            x21 = xe(2) - xe(1);         y21 = ye(2) - ye(1); 
            x32 = xe(3) - xe(2);         y32 = ye(3) - ye(2);
            x13 = xe(1) - xe(3);         y13 = ye(1) - ye(3);    
            x1=xe(1);
            x2=xe(2);
            x3=xe(3);
            y1=ye(1);
            y2=ye(2);
            y3=ye(3);
            xi_gl  = x_gl(pp);            eta_gl = e_gl(pp);

            %det_Je(pp) = det(Je);
            det_Je(pp) =x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2;
            %% Se ensambla la matriz de funciones de forma N
            N{e,pp} =[[ 1 - xi_gl - eta_gl,                  0,                  0, xi_gl,     0,     0, eta_gl,      0,      0];
                       [                  0, 1 - xi_gl - eta_gl,                  0,     0, xi_gl,     0,      0, eta_gl,      0];
                       [                  0,                  0, 1 - xi_gl - eta_gl,     0,     0, xi_gl,      0,      0, eta_gl]];
            %% Se calcula An
            % Ecuacion 22b          
            phi_k = (2/((5/6)*(1 - nu))) .* (h./Lk).^2;
            
            phi4=phi_k(1);
            phi5=phi_k(2);
            phi6=phi_k(3);
            % Ecuacion 38
            A_dbeta = diag((2/3) * Lk .* (1+phi_k));

            % Ecuacion 39
            Aw = [[  1, -x21/2, -y21/2, -1, -x21/2, -y21/2,  0,      0,      0]
                  [  0,      0,      0,  1, -x32/2, -y32/2, -1, -x32/2, -y32/2]
                  [ -1, -x13/2, -y13/2,  0,      0,      0,  1, -x13/2, -y13/2]];

            % Ecuacion 37
            An = A_dbeta\Aw;
           %% Se calcula la matriz de deformacion por flexion Bb (eq. 41)
            Bb_dbeta =[[                                                                                                             (4*C4*(y3 - y1 + eta_gl*y1 - eta_gl*y3 + xi_gl*y1 + xi_gl*y2 - 2*xi_gl*y3))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),                                                                   -(4*C5*(eta_gl*y1 - eta_gl*y3 - xi_gl*y1 + xi_gl*y2))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),                                                                                                              -(4*C6*(y2 - y1 + eta_gl*y1 - 2*eta_gl*y2 + eta_gl*y3 + xi_gl*y1 - xi_gl*y2))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)];
                       [                                                                                                            -(4*S4*(x3 - x1 + eta_gl*x1 - eta_gl*x3 + x1*xi_gl + x2*xi_gl - 2*x3*xi_gl))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),                                                                    (4*S5*(eta_gl*x1 - eta_gl*x3 - x1*xi_gl + x2*xi_gl))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),                                                                                                               (4*S6*(x2 - x1 + eta_gl*x1 - 2*eta_gl*x2 + eta_gl*x3 + x1*xi_gl - x2*xi_gl))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)];
                       [ (4*(C4*x1 - C4*x3 - S4*y1 + S4*y3 - C4*eta_gl*x1 + C4*eta_gl*x3 + S4*eta_gl*y1 - S4*eta_gl*y3 - C4*x1*xi_gl - C4*x2*xi_gl + 2*C4*x3*xi_gl + S4*xi_gl*y1 + S4*xi_gl*y2 - 2*S4*xi_gl*y3))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2), (4*(C5*eta_gl*x1 - C5*eta_gl*x3 - S5*eta_gl*y1 + S5*eta_gl*y3 - C5*x1*xi_gl + C5*x2*xi_gl + S5*xi_gl*y1 - S5*xi_gl*y2))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2), -(4*(C6*x1 - C6*x2 - S6*y1 + S6*y2 - C6*eta_gl*x1 + 2*C6*eta_gl*x2 - C6*eta_gl*x3 + S6*eta_gl*y1 - 2*S6*eta_gl*y2 + S6*eta_gl*y3 - C6*x1*xi_gl + C6*x2*xi_gl + S6*xi_gl*y1 - S6*xi_gl*y2))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)]];

            Bb_beta =[[ 0,  (y2 - y3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),                                                          0, 0, -(y1 - y3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),                                                          0, 0,  (y1 - y2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),                                                          0];
                      [ 0,                                                          0, -(x2 - x3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2), 0,                                                          0,  (x1 - x3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2), 0,                                                          0, -(x1 - x2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)];
                      [ 0, -(x2 - x3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),  (y2 - y3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2), 0,  (x1 - x3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2), -(y1 - y3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2), 0, -(x1 - x2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),  (y1 - y2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)]];
            %Bb=Bb_beta + Bb_dbeta*An;
            Bb{e,pp} = Bb_beta + Bb_dbeta*An;
            
            %% Se calcula An


            A_dbeta = diag((2/3) * Lk .* (1+phi_k));
            
            % Ecuacion 39
            Aw = [[  1, -x21/2, -y21/2, -1, -x21/2, -y21/2,  0,      0,      0]
                  [  0,      0,      0,  1, -x32/2, -y32/2, -1, -x32/2, -y32/2]
                  [ -1, -x13/2, -y13/2,  0,      0,      0,  1, -x13/2, -y13/2]];
                              
            % Ecuacion 37
            An = A_dbeta\Aw;
            
            %% Se calcula la matriz de deformacion por flexion Bb (eq. 41)
            Bb{e,pp} = Bb_beta + Bb_dbeta*An;
            
            %% Se calcula la matriz de deformacion por cortante Bs
            L5_phi5 = Lk(1)*phi_k(1);    
            L6_phi6 = Lk(2)*phi_k(2);
            L7_phi7 = Lk(3)*phi_k(3);    
            %L8_phi8 = Lk(4)*phi_k(4);
                       
            % Ecuacion 27
            Bs_dbeta = [[  (2*phi4*((S6*(eta_gl + xi_gl - 1))/(C4*S6 - C6*S4) - (S5*xi_gl)/(C4*S5 - C5*S4)))/3, -(2*phi5*((S6*eta_gl)/(C5*S6 - C6*S5) - (S4*xi_gl)/(C4*S5 - C5*S4)))/3, -(2*phi6*((S4*(eta_gl + xi_gl - 1))/(C4*S6 - C6*S4) - (S5*eta_gl)/(C5*S6 - C6*S5)))/3];
                        [ -(2*phi4*((C6*(eta_gl + xi_gl - 1))/(C4*S6 - C6*S4) - (C5*xi_gl)/(C4*S5 - C5*S4)))/3,  (2*phi5*((C6*eta_gl)/(C5*S6 - C6*S5) - (C4*xi_gl)/(C4*S5 - C5*S4)))/3,  (2*phi6*((C4*(eta_gl + xi_gl - 1))/(C4*S6 - C6*S4) - (C5*eta_gl)/(C5*S6 - C6*S5)))/3]];

            % Ecuacion 43
            Bs{e,pp} = Bs_dbeta*An;
            
            %% se arma la matriz de rigidez del elemento e por flexion (eq. 45)
            Kbe = Kbe + Bb{e,pp}'*Hb*Bb{e,pp}*det_Je(pp)*w_gl(pp);
            
            %% se arma la matriz de rigidez del elemento e por cortante (eq. 47)
            Kse = Kse + Bs{e,pp}'*Hs*Bs{e,pp}*det_Je(pp)*w_gl(pp);
            %% se arma la matriz de masa del elemento e  (eq. 47)
            
            Me = Me + N{e,pp}'*T*N{e,pp}*det_Je(pp)*w_gl(pp);
            
            %% Fundacion Winkler corregida
            % El balasto actua solamente sobre w y no se multiplica por h.
            Ntri = [1-xi_gl-eta_gl, xi_gl, eta_gl];
            Nw = [Ntri(1),0,0, Ntri(2),0,0, Ntri(3),0,0];
            He = He + Nw'*Balastro*Nw*det_Je(pp)*w_gl(pp);

            %% Presion total evaluada en el punto de Gauss
            xgp = Ntri*xe;
            ygp = Ntri*ye;

            if in_col(e)
                qPM_gp = P/Ac ...
                       + signoMy*My/Iy_col*(xgp-xc) ...
                       + signoMx*Mx/Ix_col*(ygp-yc);
                q_gp = qPM_gp + qzapata + qpedestal;
            else
                q_gp = qzapata + qsuelo;
            end

            %% Vector de fuerzas nodales equivalentes consistente
            fe = fe + Nw'*q_gp*det_Je(pp)*w_gl(pp);
        %end
    end
    %% se verifica que todos los determinantes sean positivos
    if any(det_Je(:) <= 0)
        error('Existen elementos con det(Je(xi,eta)) <= 0 %d.\n', e);
    end
    
    %% ensamblaje matricial
    idx{e} = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:)];    
    Kplaca(idx{e},idx{e}) = Kplaca(idx{e},idx{e}) + Kbe + Kse;
    Ksuelo(idx{e},idx{e}) = Ksuelo(idx{e},idx{e}) + He;
    M(idx{e},idx{e}) = M(idx{e},idx{e}) + Me;
    f(idx{e},:)      = f(idx{e},:)      + fe;
end

%% Matriz total y reaccion de suelo
K = Kplaca + Ksuelo;
K = 0.5*(K+K');

%% Muestro la configuracion de la matriz K (K es rala)
figure
spy(K);
title('Los puntos representan los elementos diferentes de cero')

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
aa = zeros(ngdl,1); aa(c) = ac;  aa(d) = ad; % desplazamientos
qq = zeros(ngdl,1); qq(c) = qd;              % fuerzas nodales equivalentes

%% Reaccion fisica de la fundacion Winkler
Rsuelo = -Ksuelo*aa;
Fw = f(1:3:end);
Rw = Rsuelo(1:3:end);

Mx_ext = sum(Fw.*(xnod(:,Y)-yc));
My_ext = sum(Fw.*(xnod(:,X)-xc));
Mx_suelo = sum(Rw.*(xnod(:,Y)-yc));
My_suelo = sum(Rw.*(xnod(:,X)-xc));

fprintf('\nEQUILIBRIO GLOBAL\n');
fprintf('F externa + R suelo = %.6e kN\n',sum(Fw)+sum(Rw));
fprintf('Mx externa + Mx suelo = %.6e kN*m\n',Mx_ext+Mx_suelo);
fprintf('My externa + My suelo = %.6e kN*m\n',My_ext+My_suelo);

vect_mov = reshape(aa,3,nno)'; % vector de movimientos

%% Dibujo la malla de elementos finitos y las deformaciones de esta
%escala = 10; % factor de escalamiento de la deformada
xdef = escala*vect_mov; % posicion de la deformada
figure;
grid on;

% Una sola superficie grafica para toda la malla.
% Esto evita crear un objeto fill3 por cada elemento.
trisurf(LaG, ...
        xnod(:,X), ...
        xnod(:,Y), ...
        xdef(:,ww), ...
        xdef(:,ww), ...
        'EdgeColor','none', ...
        'FaceColor','interp');

daspect([1 1 1]); % similar a axis equal, pero en 3D
axis tight
colormap jet
colorbar
title(sprintf('Deformada escalada %d veces',escala),'FontSize',20)
view(3)

% Reduce las actualizaciones innecesarias de la ventana grafica.
try
    drawnow limitrate
catch
    drawnow
end

%% Se calcula para cada elemento el vector de momentos en los puntos
%% de Gauss (ecuacion 49)
MxMyMxy = cell(nef,n_gl,n_gl);
for e = 1:nef
    for pp = 1:n_gl
        %for qq = 1:n_gl
            MxMyMxy{e,pp} = Hb*Bb{e,pp}*aa(idx{e});
        %end
    end
end

%% Se calcula para cada elemento el vector de cortantes en los puntos
%% de Gauss (ecuacion 50)
if EFtype == DKMQ
    QxQy = cell(nef,n_gl,n_gl);
    for e = 1:nef
        for pp = 1:n_gl
            %for qq = 1:n_gl
                QxQy{e,pp} = Hs*Bs{e,pp}*aa(idx{e});
            %end
        end
    end
end

%% Se extrapolan los momentos y cortantes a los nodos
num_elem_ady = zeros(nno,1);  % numero de elementos adyacentes
Mx  = zeros(nno,1);
My  = zeros(nno,1);
Mxy = zeros(nno,1);
Qx  = zeros(nno,1);
Qy  = zeros(nno,1);
    A =[0.126340726488392,-0.638559587411924,-0.638559587411923,1.87365927351158,0.138559587411936,0.138559587411936;
        -0.638559587411907,-0.638559587411945,0.126340726488378,0.138559587411937,0.138559587411936,1.87365927351160;
        -0.638559587411908,0.126340726488377,-0.638559587411944,0.138559587411937,1.87365927351160,0.138559587411936];

for e = 1:nef  
    Mx1=[MxMyMxy{e,:,:}];
    Mx(LaG(e,:),:) = Mx(LaG(e,:),:)   + A * Mx1(1,:)';
   
    My1=[MxMyMxy{e,:,:}];
    My(LaG(e,:),:) = My(LaG(e,:),:)   + A *My1(2,:)';
                                                                              
    Mxy1=[MxMyMxy{e,:,:}];                                     
    Mxy(LaG(e,:),:) = Mxy(LaG(e,:),:) + A * Mxy1(3,:)';
                                          
   num_elem_ady(LaG(e,:),:) = num_elem_ady(LaG(e,:),:) + 1;
end 

if EFtype == DKMQ
    for e = 1:nef                             
        Qx1=[QxQy{e,:,:}];
        Qx(LaG(e,:),:) = Qx(LaG(e,:),:)   + A *  Qx1(1,:)';
        Qy1=[QxQy{e,:,:}];
        Qy(LaG(e,:),:) = Qy(LaG(e,:),:)   + A *   Qy1(2,:)';
    end
end 
 
%% Alisado (promedio de los momentos y cortantes en los nodos)
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
if EFtype == DKMQ
    figure
    subplot(1,2,1); plot_M_or_Q(nef, xnod, LaG, Qx,  'Cortantes Qx (N/m)');
    subplot(1,2,2); plot_M_or_Q(nef, xnod, LaG, Qy,  'Cortantes Qy (N/m)');
end
    
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
if EFtype == DKMQ
    Q_max = hypot(Qx, Qy);
    ang   = atan2(Qy, Qx);
    figure
    plot_M_or_Q(nef, xnod, LaG, Q_max, 'Q_{max} (N/m)', { ang })
end

%% Se calculan los momentos de disenio de Wood y Armer
[Mxast_sup, Myast_sup, Mxast_inf, Myast_inf] = arrayfun(@WoodArmer, Mx, My, Mxy);
Mmax = max(abs([Mxast_sup; Myast_sup; Mxast_inf; Myast_inf]));

% se graficaran los momentos de disenio utilizando la misma escala de
% colores en valor absoluto, de este modo Mxast_sup=+100 y Mxast_inf=-100 
% tendran el mismo color
figure
subplot(1,2,1); plot_M_or_Q(nef, xnod, LaG, Mxast_sup,  'Momentos M_x^* sup (N-m/m)');
caxis([0 Mmax]);                                % misma escala de colores
colorbar('ylim', [0 max(Mxast_sup)]);           % rango de colores a mostrar

subplot(1,2,2); plot_M_or_Q(nef, xnod, LaG, Myast_sup,  'Momentos M_y^* sup (N-m/m)');
caxis([0 Mmax]);                                % misma escala de colores
colorbar('ylim', [0 max(Myast_sup)]);           % rango de colores a mostrar


figure
subplot(1,2,1); plot_M_or_Q(nef, xnod, LaG, Mxast_inf,  'Momentos M_x^* inf (N-m/m)');
caxis([-Mmax 0]);                               % misma escala de colores
oldcmap = colormap; colormap(flipud(oldcmap));  % invierto mapa de colores
colorbar('ylim', [min(Mxast_inf) 0]);           % rango de colores a mostrar
subplot(1,2,2); plot_M_or_Q(nef, xnod, LaG, Myast_inf,  'Momentos M_y^* inf (N-m/m)');
caxis([-Mmax 0]);                               % misma escala de colores
oldcmap = colormap; colormap(flipud(oldcmap));  % invierto mapa de colores
colorbar('ylim', [min(Myast_inf) 0]);           % rango de colores a mostrar



%%
return; % bye, bye!

%%
function plot_M_or_Q(nef, xnod, LaG, variable, texto, angulos)
%PLOT_M_OR_Q
% Grafica rapidamente una variable nodal sobre una malla triangular.
%
% Se conserva la firma original para no modificar las llamadas existentes:
%
%   plot_M_or_Q(nef,xnod,LaG,variable,texto)
%   plot_M_or_Q(nef,xnod,LaG,variable,texto,angulos)
%
% La mejora principal consiste en sustituir el ciclo de "fill", que creaba
% un objeto grafico por elemento, por una unica llamada a "patch".
%
% "nef" se conserva por compatibilidad, aunque ya no es necesario para
% construir la superficie.

    %#ok<INUSD>
    X = 1;
    Y = 2;

    variable = variable(:);

    if length(variable) ~= size(xnod,1)
        error(['La variable que se desea graficar debe tener un valor ' ...
               'por cada nodo de la malla.']);
    end

    % Una sola llamada grafica para todos los elementos.
    patch('Faces',LaG, ...
          'Vertices',xnod(:,[X Y]), ...
          'FaceVertexCData',variable, ...
          'FaceColor','interp', ...
          'EdgeColor','none');

    axis equal tight
    box on
    colormap jet
    colorbar
    title(texto,'FontSize',20);

    % Las direcciones principales pueden producir miles de segmentos.
    % Para mantener la grafica fluida se limita automaticamente el numero
    % de nodos en los que se dibujan.
    if nargin == 6 && ~isempty(angulos)

        hold on

        nno = size(xnod,1);

        % Se muestran como maximo unas 400 direcciones.
        max_flechas = 400;
        salto = max(1,ceil(nno/max_flechas));
        nod_graf = 1:salto:nno;

        esc = 0.5;
        norma = 1;

        for i = 1:length(angulos)

            ang_i = angulos{i};

            if length(ang_i) ~= nno
                error(['Cada vector de angulos debe contener un valor ' ...
                       'por cada nodo.']);
            end

            ang_i = ang_i(:);

            % Direccion positiva.
            quiver(xnod(nod_graf,X), ...
                   xnod(nod_graf,Y), ...
                   norma*cos(ang_i(nod_graf)), ...
                   norma*sin(ang_i(nod_graf)), ...
                   esc, ...
                   'k', ...
                   'ShowArrowHead','off', ...
                   'LineWidth',1, ...
                   'Marker','.');

            % Misma direccion girada 180 grados.
            quiver(xnod(nod_graf,X), ...
                   xnod(nod_graf,Y), ...
                   norma*cos(ang_i(nod_graf)+pi), ...
                   norma*sin(ang_i(nod_graf)+pi), ...
                   esc, ...
                   'k', ...
                   'ShowArrowHead','off', ...
                   'LineWidth',1, ...
                   'Marker','.');

        end

        hold off
    end

    try
        drawnow limitrate
    catch
        drawnow
    end
end
