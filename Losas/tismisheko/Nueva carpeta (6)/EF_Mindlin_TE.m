function resultado = EF_Mindlin_TE(xnod,LaG,gdl,c,d,E,nu,h,q,rho, ...
    kWinkler,Nx0,Ny0,Nxy0,xqi,xqf,yqi,yqf,escala)
%==========================================================================
% LOSA MINDLIN Q4 CON FUNCIONES DE FORMA OBTENIDAS A PARTIR DE
% VIGAS DE TIMOSHENKO EN LOS CUATRO BORDES.
%
% Se conserva la estructura convencional de integracion de Gauss:
%
%   Kb = sum(Bb'*Hb*Bb*detJ*w_i*w_j)
%   Ks = sum(Bs'*Hs*Bs*detJ*w_i*w_j)
%   Kw = sum(Nw'*kWinkler*Nw*detJ*w_i*w_j)
%   Kg = sum(BG'*N0*BG*detJ*w_i*w_j)
%
% Convencion:
%   Nx0, Ny0 > 0  -> compresion
%   Ke = Kb + Ks + Kw - Kg
%==========================================================================

%% constantes
X = 1; Y = 2; Z = 3;
ww = 1; tx = 2; ty = 3;

%% ========================================================================
% OPCIONES DE GRAFICA
% ========================================================================
graficar_malla       = true;
graficar_matriz      = true;
graficar_deformada   = true;
graficar_momentos    = true;
graficar_cortantes   = true;
graficar_principales = true;
graficar_qmax        = true;
graficar_woodarmer   = true;
graficar_movimientos = true;

% La numeracion de nodos y elementos es costosa para mallas grandes.
mostrar_num_nodos     = false;
mostrar_num_elementos = false;

% Numero maximo aproximado de vectores para graficas de direcciones.
max_vectores = 300;

nno  = size(xnod,1);
ngdl = 3*nno;
nef  = size(LaG,1);

%% ========================================================================
% Se dibuja la malla de elementos finitos - VERSION RAPIDA
%
% Se utiliza un unico objeto PATCH para toda la malla en lugar de crear
% un objeto LINE por cada elemento. La numeracion se deja como opcion.
% ========================================================================
if graficar_malla
    figure;
    hold on;
    box on;

    patch('Faces',LaG, ...
          'Vertices',xnod(:,1:2), ...
          'FaceColor','none', ...
          'EdgeColor',[0.35 0.35 0.35], ...
          'LineWidth',0.35);

    plot(xnod(:,X),xnod(:,Y),'.r','MarkerSize',5);

    if mostrar_num_elementos
        cgx = mean(reshape(xnod(LaG',X),4,nef),1);
        cgy = mean(reshape(xnod(LaG',Y),4,nef),1);
        text(cgx(:),cgy(:),num2str((1:nef)'), ...
             'Color','b','FontSize',6, ...
             'HorizontalAlignment','center');
    end

    if mostrar_num_nodos
        text(xnod(:,X),xnod(:,Y),num2str((1:nno)'), ...
             'FontSize',6,'Color',[0.5 0 0]);
    end

    axis equal tight
    title('Malla de elementos finitos - Mindlin TE');
    drawnow limitrate
end

%% ========================================================================
% Cuadratura de Gauss-Legendre
% ========================================================================
n_gl = 2;
[x_gl,w_gl] = gausslegendre_quad(n_gl);

%% ========================================================================
% Funciones Q4 para la geometria
% ========================================================================
Nforma = @(xi,eta) [1/4*(1-xi)*(1-eta)
                    1/4*(1+xi)*(1-eta)
                    1/4*(1+xi)*(1+eta)
                    1/4*(1-xi)*(1+eta)];

dN_dxi = @(xi,eta) [-1/4*(1-eta)
                      1/4*(1-eta)
                      1/4*(1+eta)
                     -1/4*(1+eta)];

dN_deta = @(xi,eta) [-1/4*(1-xi)
                      -1/4*(1+xi)
                       1/4*(1+xi)
                       1/4*(1-xi)];

%% ========================================================================
% Matrices constitutivas
% ========================================================================
Db0 = E*h^3/(12*(1-nu^2));
Hb = Db0*[1 nu 0
          nu 1 0
          0 0 (1-nu)/2];

G = E/(2*(1+nu));
Hs = (5/6)*G*h*eye(2);

% Estado membranal inicial para matriz geometrica
N0 = [Nx0 Nxy0
      Nxy0 Ny0];

% Matriz de masa: se conserva tu criterio original
T = rho*diag([h,0,0]);

%% ========================================================================
% Ensamblaje global
% ========================================================================
K = sparse(ngdl,ngdl);
M = sparse(ngdl,ngdl);
f = zeros(ngdl,1);

Kb_global = sparse(ngdl,ngdl);
Ks_global = sparse(ngdl,ngdl);
Kw_global = sparse(ngdl,ngdl);
Kg_global = sparse(ngdl,ngdl);

N  = cell(nef,n_gl,n_gl);
Bb = cell(nef,n_gl,n_gl);
Bs = cell(nef,n_gl,n_gl);
BG = cell(nef,n_gl,n_gl);
idx = cell(nef,1);

%% ========================================================================
% Ciclo sobre elementos
% ========================================================================
for e = 1:nef
    xe = xnod(LaG(e,:),X);
    ye = xnod(LaG(e,:),Y);

    % Esta primera formulacion usa rectangulos alineados con los ejes.
    Lxe = hypot(xe(2)-xe(1),ye(2)-ye(1));
    Lye = hypot(xe(4)-xe(1),ye(4)-ye(1));

    Kbe = zeros(12);
    Kse = zeros(12);
    Kwe = zeros(12);
    Kge = zeros(12);
    Me  = zeros(12);
    fe  = zeros(12,1);
    det_Je = zeros(n_gl,n_gl);

    for pp = 1:n_gl
        for qq = 1:n_gl
            xi_gl  = x_gl(pp);
            eta_gl = x_gl(qq);

            %% funciones Q4 geometricas
            NN = Nforma(xi_gl,eta_gl);
            ddN_dxi  = dN_dxi(xi_gl,eta_gl);
            ddN_deta = dN_deta(xi_gl,eta_gl);

            %% Jacobiano
            dx_dxi  = sum(ddN_dxi.*xe);
            dy_dxi  = sum(ddN_dxi.*ye);
            dx_deta = sum(ddN_deta.*xe);
            dy_deta = sum(ddN_deta.*ye);

            Je = [dx_dxi  dy_dxi
                  dx_deta dy_deta];
            det_Je(pp,qq) = det(Je);

            if det_Je(pp,qq) <= 0
                error('Existen elementos con det(Je)<=0. Elemento %d.',e);
            end

            %% funciones de forma Mindlin-Timoshenko
            [Nw,Ntx,Nty,dNw_dx,dNw_dy,dNtx_dx,dNtx_dy,dNty_dx,dNty_dy] = ...
                funciones_forma_Mindlin_TE(xi_gl,eta_gl,Lxe,Lye,Db0,(5/6)*G*h);

            %% matriz N completa [w tx ty]
            N{e,pp,qq} = [Nw;Ntx;Nty];

            %% matriz de flexion
            Bb{e,pp,qq} = [dNtx_dx
                            dNty_dy
                            dNtx_dy + dNty_dx];

            %% matriz de cortante
            Bs{e,pp,qq} = [Ntx - dNw_dx
                            Nty - dNw_dy];

            %% matriz geometrica
            BG{e,pp,qq} = [dNw_dx
                            dNw_dy];

            peso = det_Je(pp,qq)*w_gl(pp)*w_gl(qq);

            %% rigideces
            Kbe = Kbe + Bb{e,pp,qq}'*Hb*Bb{e,pp,qq}*peso;
            Kse = Kse + Bs{e,pp,qq}'*Hs*Bs{e,pp,qq}*peso;
            Kwe = Kwe + Nw'*kWinkler*Nw*peso;
            Kge = Kge + BG{e,pp,qq}'*N0*BG{e,pp,qq}*peso;

            %% masa
            Me = Me + N{e,pp,qq}'*T*N{e,pp,qq}*peso;

            %% coordenadas fisicas del punto de Gauss
            xpg = NN.'*xe(:);
            ypg = NN.'*ye(:);

            %% vector de fuerzas nodales equivalentes
            if (xpg >= xqi && xpg <= xqf) && (ypg >= yqi && ypg <= yqf)
                fe = fe + Nw'*q*peso;
            end
        end
    end

    %% matriz total del elemento
    Ke = Kbe + Kse + Kwe - Kge;

    %% simetrizacion numerica
    Kbe = 0.5*(Kbe+Kbe.');
    Kse = 0.5*(Kse+Kse.');
    Kwe = 0.5*(Kwe+Kwe.');
    Kge = 0.5*(Kge+Kge.');
    Ke  = 0.5*(Ke+Ke.');

    %% ensamblaje
    idx{e} = [gdl(LaG(e,1),:) gdl(LaG(e,2),:) ...
              gdl(LaG(e,3),:) gdl(LaG(e,4),:)];

    K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke;
    M(idx{e},idx{e}) = M(idx{e},idx{e}) + Me;
    f(idx{e}) = f(idx{e}) + fe;

    Kb_global(idx{e},idx{e}) = Kb_global(idx{e},idx{e}) + Kbe;
    Ks_global(idx{e},idx{e}) = Ks_global(idx{e},idx{e}) + Kse;
    Kw_global(idx{e},idx{e}) = Kw_global(idx{e},idx{e}) + Kwe;
    Kg_global(idx{e},idx{e}) = Kg_global(idx{e},idx{e}) + Kge;
end

%% ========================================================================
% Configuracion de la matriz global
% ========================================================================
if graficar_matriz
    figure
    spy(K)
    title('Los puntos representan los elementos diferentes de cero')
    drawnow limitrate
end

%% ========================================================================
% Resolucion del sistema
% ========================================================================
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

ac = zeros(length(c),1);
ad = Kdd\(fc-Kdc*ac);
qd = Kcc*ac + Kcd*ad - fd;

aa = zeros(ngdl,1); aa(c)=ac; aa(d)=ad;
qq = zeros(ngdl,1); qq(c)=qd;

vect_mov = reshape(aa,3,nno)';

%% ========================================================================
% Deformada - VERSION RAPIDA
%
% Se dibuja toda la superficie deformada mediante un unico PATCH.
% ========================================================================
xdef = escala*vect_mov;
if graficar_deformada
    figure
    hold on
    grid on
    box on

    zdef = xdef(:,ww);

    patch('Faces',LaG, ...
          'Vertices',[xnod(:,X),xnod(:,Y),zdef], ...
          'FaceVertexCData',zdef, ...
          'FaceColor','interp', ...
          'EdgeColor',[0.35 0.35 0.35], ...
          'LineWidth',0.20);

    daspect([1 1 1])
    axis tight
    colorbar
    title(sprintf('Deformada escalada %g veces',escala),'FontSize',14)
    xlabel('x')
    ylabel('y')
    zlabel('w')
    view(3)
    drawnow limitrate
end

%% ========================================================================
% Momentos en puntos de Gauss
% ========================================================================
MxMyMxy = cell(nef,n_gl,n_gl);
for e = 1:nef
    for pp = 1:n_gl
        for qq = 1:n_gl
            MxMyMxy{e,pp,qq} = Hb*Bb{e,pp,qq}*aa(idx{e});
        end
    end
end

%% ========================================================================
% Cortantes en puntos de Gauss
% ========================================================================
QxQy = cell(nef,n_gl,n_gl);
for e = 1:nef
    for pp = 1:n_gl
        for qq = 1:n_gl
            QxQy{e,pp,qq} = Hs*Bs{e,pp,qq}*aa(idx{e});
        end
    end
end

%% ========================================================================
% Extrapolacion a nodos
% ========================================================================
num_elem_ady = zeros(nno,1);
Mx  = zeros(nno,1);
My  = zeros(nno,1);
Mxy = zeros(nno,1);
Qx  = zeros(nno,1);
Qy  = zeros(nno,1);

Aext = [ ...
   sqrt(3)/2 + 1,          -1/2,          -1/2, 1 - sqrt(3)/2
          -1/2, 1 - sqrt(3)/2, sqrt(3)/2 + 1,          -1/2
   1 - sqrt(3)/2,          -1/2,          -1/2, sqrt(3)/2 + 1
          -1/2, sqrt(3)/2 + 1, 1 - sqrt(3)/2,          -1/2];

for e = 1:nef
    Mx(LaG(e,:)) = Mx(LaG(e,:)) + Aext*[MxMyMxy{e,1,1}(1);MxMyMxy{e,1,2}(1);MxMyMxy{e,2,1}(1);MxMyMxy{e,2,2}(1)];
    My(LaG(e,:)) = My(LaG(e,:)) + Aext*[MxMyMxy{e,1,1}(2);MxMyMxy{e,1,2}(2);MxMyMxy{e,2,1}(2);MxMyMxy{e,2,2}(2)];
    Mxy(LaG(e,:)) = Mxy(LaG(e,:)) + Aext*[MxMyMxy{e,1,1}(3);MxMyMxy{e,1,2}(3);MxMyMxy{e,2,1}(3);MxMyMxy{e,2,2}(3)];

    Qx(LaG(e,:)) = Qx(LaG(e,:)) + Aext*[QxQy{e,1,1}(1);QxQy{e,1,2}(1);QxQy{e,2,1}(1);QxQy{e,2,2}(1)];
    Qy(LaG(e,:)) = Qy(LaG(e,:)) + Aext*[QxQy{e,1,1}(2);QxQy{e,1,2}(2);QxQy{e,2,1}(2);QxQy{e,2,2}(2)];

    num_elem_ady(LaG(e,:)) = num_elem_ady(LaG(e,:)) + 1;
end

Mx  = Mx./num_elem_ady;
My  = My./num_elem_ady;
Mxy = Mxy./num_elem_ady;
Qx  = Qx./num_elem_ady;
Qy  = Qy./num_elem_ady;

%% ========================================================================
% Graficas de momentos
% ========================================================================
if graficar_momentos
    figure
    subplot(1,3,1); plot_M_or_Q(nef,xnod,LaG,Mx,'Momentos Mx (kN-m/m)',[],max_vectores);
    subplot(1,3,2); plot_M_or_Q(nef,xnod,LaG,My,'Momentos My (kN-m/m)',[],max_vectores);
    subplot(1,3,3); plot_M_or_Q(nef,xnod,LaG,Mxy,'Momentos Mxy (kN-m/m)',[],max_vectores);
    drawnow limitrate
end

%% ========================================================================
% Graficas de cortantes
% ========================================================================
if graficar_cortantes
    figure
    subplot(1,2,1); plot_M_or_Q(nef,xnod,LaG,Qx,'Cortantes Qx (kN/m)',[],max_vectores);
    subplot(1,2,2); plot_M_or_Q(nef,xnod,LaG,Qy,'Cortantes Qy (kN/m)',[],max_vectores);
    drawnow limitrate
end

%% ========================================================================
% Momentos principales y direcciones
% ========================================================================
Mt_max = sqrt(((Mx-My)/2).^2 + Mxy.^2);
Mf1_xy = (Mx+My)/2 + Mt_max;
Mf2_xy = (Mx+My)/2 - Mt_max;
angM = 0.5*atan2(2*Mxy,Mx-My);

if graficar_principales
    figure
    subplot(1,3,1); plot_M_or_Q(nef,xnod,LaG,Mf1_xy,'Mf1_{xy} (kN-m/m)',{angM},max_vectores)
    subplot(1,3,2); plot_M_or_Q(nef,xnod,LaG,Mf2_xy,'Mf2_{xy} (kN-m/m)',{angM+pi/2},max_vectores)
    subplot(1,3,3); plot_M_or_Q(nef,xnod,LaG,Mt_max,'Mt_{max} (kN-m/m)',{angM+pi/4,angM-pi/4},max_vectores)
    drawnow limitrate
end

%% ========================================================================
% Cortante maximo y direccion
% ========================================================================
Q_max = hypot(Qx,Qy);
angQ = atan2(Qy,Qx);
if graficar_qmax
    figure
    plot_M_or_Q(nef,xnod,LaG,Q_max,'Q_{max} (kN/m)',{angQ},max_vectores)
    drawnow limitrate
end

%% ========================================================================
% Wood-Armer
% ========================================================================
[Mxast_sup,Myast_sup,Mxast_inf,Myast_inf] = arrayfun(@WoodArmer,Mx,My,Mxy);
Mmax = max(abs([Mxast_sup;Myast_sup;Mxast_inf;Myast_inf]));

if graficar_woodarmer
    figure
    subplot(1,2,1); plot_M_or_Q(nef,xnod,LaG,Mxast_sup,'Momentos M_x^* sup (kN-m/m)',[],max_vectores);
    if Mmax > 0, caxis([0 Mmax]); end
    subplot(1,2,2); plot_M_or_Q(nef,xnod,LaG,Myast_sup,'Momentos M_y^* sup (kN-m/m)',[],max_vectores);
    if Mmax > 0, caxis([0 Mmax]); end

    figure
    subplot(1,2,1); plot_M_or_Q(nef,xnod,LaG,Mxast_inf,'Momentos M_x^* inf (kN-m/m)',[],max_vectores);
    if Mmax > 0, caxis([-Mmax 0]); end
    subplot(1,2,2); plot_M_or_Q(nef,xnod,LaG,Myast_inf,'Momentos M_y^* inf (kN-m/m)',[],max_vectores);
    if Mmax > 0, caxis([-Mmax 0]); end
    drawnow limitrate
end

%% ========================================================================
% Graficas adicionales de desplazamiento y giros nodales
% ========================================================================
if graficar_movimientos
    figure
    subplot(1,3,1); plot_M_or_Q(nef,xnod,LaG,vect_mov(:,ww),'Desplazamiento w (m)',[],max_vectores);
    subplot(1,3,2); plot_M_or_Q(nef,xnod,LaG,vect_mov(:,tx),'Giro theta_x (rad)',[],max_vectores);
    subplot(1,3,3); plot_M_or_Q(nef,xnod,LaG,vect_mov(:,ty),'Giro theta_y (rad)',[],max_vectores);
    drawnow limitrate
end

%% ========================================================================
% Resultados principales
% ========================================================================
[wmin,iwmin] = min(vect_mov(:,ww));
[wmax,iwmax] = max(vect_mov(:,ww));

fprintf('\n============================================================\n')
fprintf(' RESULTADOS PRINCIPALES\n')
fprintf('============================================================\n')
fprintf('Numero de nodos      = %d\n',nno)
fprintf('Numero de elementos  = %d\n',nef)
fprintf('Numero de GDL        = %d\n',ngdl)
fprintf('w minimo             = %.12e   nodo %d\n',wmin,iwmin)
fprintf('w maximo             = %.12e   nodo %d\n',wmax,iwmax)
fprintf('||Kb||_F             = %.12e\n',norm(Kb_global,'fro'))
fprintf('||Ks||_F             = %.12e\n',norm(Ks_global,'fro'))
fprintf('||Kw||_F             = %.12e\n',norm(Kw_global,'fro'))
fprintf('||Kg||_F             = %.12e\n',norm(Kg_global,'fro'))
fprintf('Error simetria K     = %.12e\n',norm(K-K.','fro'))

%% ========================================================================
% Salidas
% ========================================================================
resultado.K = K;
resultado.Kb = Kb_global;
resultado.Ks = Ks_global;
resultado.Kw = Kw_global;
resultado.Kg = Kg_global;
resultado.M = M;
resultado.f = f;
resultado.a = aa;
resultado.q = qq;
resultado.vect_mov = vect_mov;
resultado.Bb = Bb;
resultado.Bs = Bs;
resultado.BG = BG;
resultado.MxMyMxy = MxMyMxy;
resultado.QxQy = QxQy;
resultado.Mx = Mx;
resultado.My = My;
resultado.Mxy = Mxy;
resultado.Qx = Qx;
resultado.Qy = Qy;
resultado.Mf1 = Mf1_xy;
resultado.Mf2 = Mf2_xy;
resultado.Mtmax = Mt_max;
resultado.Qmax = Q_max;
resultado.Mxast_sup = Mxast_sup;
resultado.Myast_sup = Myast_sup;
resultado.Mxast_inf = Mxast_inf;
resultado.Myast_inf = Myast_inf;
resultado.idx = idx;

end

%% ========================================================================
% Funcion grafica - VERSION RAPIDA
% ========================================================================
function plot_M_or_Q(nef,xnod,LaG,variable,texto,angulos,max_vectores)
%=========================================================================
% PLOT_M_OR_Q
%
% Version optimizada para mallas Q4.
%
% La version anterior creaba un objeto FILL por cada elemento. En mallas
% grandes eso genera cientos o miles de objetos graficos y hace muy lenta
% la visualizacion. Aqui toda la malla se representa mediante un unico
% objeto PATCH.
%
% Si se suministran angulos, los vectores principales se submuestrean para
% evitar crear miles de flechas QUIVER.
%=========================================================================

    X = 1;
    Y = 2;

    if nargin < 6
        angulos = [];
    end

    if nargin < 7 || isempty(max_vectores)
        max_vectores = 300;
    end

    variable = variable(:);

    hold on
    box on

    %% ====================================================================
    % UN SOLO PATCH PARA TODA LA MALLA
    % ====================================================================
    patch('Faces',LaG, ...
          'Vertices',xnod(:,1:2), ...
          'FaceVertexCData',variable, ...
          'FaceColor','interp', ...
          'EdgeColor','none');

    axis equal tight
    title(texto,'FontSize',13)
    colorbar

    %% ====================================================================
    % DIRECCIONES PRINCIPALES / CORTANTES
    % ====================================================================
    if ~isempty(angulos)

        nno = size(xnod,1);
        salto = max(1,ceil(nno/max_vectores));
        nodos_plot = 1:salto:nno;

        Lx = max(xnod(:,X))-min(xnod(:,X));
        Ly = max(xnod(:,Y))-min(xnod(:,Y));
        Lref = 0.025*max(Lx,Ly);

        for ii = 1:length(angulos)

            ang = angulos{ii};
            ang = ang(:);

            quiver(xnod(nodos_plot,X), ...
                   xnod(nodos_plot,Y), ...
                   Lref*cos(ang(nodos_plot)), ...
                   Lref*sin(ang(nodos_plot)), ...
                   0,'k', ...
                   'ShowArrowHead','off', ...
                   'LineWidth',0.65);

            quiver(xnod(nodos_plot,X), ...
                   xnod(nodos_plot,Y), ...
                  -Lref*cos(ang(nodos_plot)), ...
                  -Lref*sin(ang(nodos_plot)), ...
                   0,'k', ...
                   'ShowArrowHead','off', ...
                   'LineWidth',0.65);
        end
    end
end
