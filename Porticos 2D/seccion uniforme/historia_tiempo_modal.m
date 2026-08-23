clc
clear
close all
format long g

%% ========================================================================
% MATRICES DEL SISTEMA
% =========================================================================
% Ejemplo Luis Enrique Garcia - pagina 288

k = [216760,-306770,105490,-19561,4282.20000000000,-510.880000000000;
    -306770,668240,-475140,137940,-29375,5385.70000000000;
     105490,-475140,731370,-493230,159600,-29327;
     -19561,137940,-493230,749020,-494470,145710;
     4282.20000000000,-29375,159600,-494470,738110,-511900;
     -510.880000000000,5385.70000000000,-29327,145710,-511900,889940];

m = 256*eye(6);

%% ========================================================================
% DATOS DEL ANALISIS
% =========================================================================

n = 6;                  % numero de modos
iter = 200;             % iteraciones de Stodola
xi = 0.05;              % amortiguamiento modal
dt = 0.02;              % paso temporal del acelerograma [s]

%% ========================================================================
% ANALISIS MODAL - SE CONSERVA TU FUNCION ESTODOLA_2
% =========================================================================

[Phi,Omega,Te] = estodola_2(k,m,iter,n);

% Periodos, suponiendo que Omega esta en rad/s
Te = 2*pi./Omega;

% Frecuencias en Hz
frecuencia = Omega/(2*pi);

%% ========================================================================
% NORMALIZACION MODAL RESPECTO A M
% =========================================================================

for i = 1:n
    Phi(:,i) = Phi(:,i)/sqrt(Phi(:,i)'*m*Phi(:,i));
end

% Matrices modales de control
Mmodal = Phi'*m*Phi;
Kmodal = Phi'*k*Phi;

%% ========================================================================
% VECTOR DE INFLUENCIA
% =========================================================================
% Para este ejemplo todos los GDL son traslaciones en la direccion del sismo.

r = ones(size(m,1),1);

%% ========================================================================
% PROPIEDADES MODALES
% =========================================================================

Mi = diag(Phi'*m*Phi);       % masa modal
Li = Phi'*m*r;               % carga modal generalizada
Gamma = Li./Mi;              % factor de participacion modal
Meff = Li.^2./Mi;            % masa modal efectiva
Mtotal = r'*m*r;             % masa total participante
participacion = Meff/Mtotal;
participacion_acumulada = cumsum(participacion);

%% ========================================================================
% ACELEROGRAMA
% =========================================================================

ag = load('ELCENTRO.TXT');
ag = ag(:);
size_ag = length(ag);
t = (0:size_ag-1)*dt;

%% ========================================================================
% RESPUESTA DE CADA OSCILADOR MODAL MEDIANTE TU RK4_45
% =========================================================================
% D     = desplazamiento relativo del oscilador unitario
% VD    = velocidad relativa
% ADD   = aceleracion relativa
% ADABS = aceleracion absoluta del oscilador unitario

D     = zeros(n,size_ag);
VD    = zeros(n,size_ag);
ADD   = zeros(n,size_ag);
ADABS = zeros(n,size_ag);

for i = 1:n
    [D(i,:),VD(i,:),ADD(i,:),ADABS(i,:)] = RK4_45_1(xi,Te(i),ag,dt);
end

%% ========================================================================
% COORDENADAS MODALES
% =========================================================================

GammaM = repmat(Gamma,1,size_ag);

q   = GammaM.*D;
qd  = GammaM.*VD;
qdd = GammaM.*ADD;

%% ========================================================================
% RECONSTRUCCION DE LA RESPUESTA FISICA
% =========================================================================

u     = Phi*q;             % desplazamiento relativo
v     = Phi*qd;            % velocidad relativa
a_rel = Phi*qdd;           % aceleracion relativa

% Aceleracion absoluta de la estructura
% a_abs = a_rel + r*ag(t)
a_abs = a_rel + r*ag.';

%% ========================================================================
% FUERZAS ELASTICAS
% =========================================================================

Fs = k*u;

%% ========================================================================
% CONTRIBUCIONES MODALES INDIVIDUALES
% =========================================================================

umodal = zeros(size(m,1),size_ag,n);
Fsmodal = zeros(size(m,1),size_ag,n);

for i = 1:n
    umodal(:,:,i) = Phi(:,i)*q(i,:);
    Fsmodal(:,:,i) = k*umodal(:,:,i);
end

%% ========================================================================
% COMPROBACION DE SUPERPOSICION MODAL
% =========================================================================

u_comprobacion = sum(umodal,3);
error_superposicion = norm(u-u_comprobacion,'fro')/(norm(u,'fro')+eps);

%% ========================================================================
% MAXIMOS Y SIGNOS
% =========================================================================

[umax,ind_umax] = max(abs(u),[],2);
[vmax,ind_vmax] = max(abs(v),[],2);
[arelmax,ind_arelmax] = max(abs(a_rel),[],2);
[aabsmax,ind_aabsmax] = max(abs(a_abs),[],2);
[Fsmax,ind_Fsmax] = max(abs(Fs),[],2);

umax_signo = zeros(size(m,1),1);
vmax_signo = zeros(size(m,1),1);
arelmax_signo = zeros(size(m,1),1);
aabsmax_signo = zeros(size(m,1),1);
Fsmax_signo = zeros(size(m,1),1);

for i = 1:size(m,1)
    umax_signo(i) = u(i,ind_umax(i));
    vmax_signo(i) = v(i,ind_vmax(i));
    arelmax_signo(i) = a_rel(i,ind_arelmax(i));
    aabsmax_signo(i) = a_abs(i,ind_aabsmax(i));
    Fsmax_signo(i) = Fs(i,ind_Fsmax(i));
end

%% ========================================================================
% RESULTADOS MODALES
% =========================================================================

fprintf('\n')
fprintf('================================================================================\n')
fprintf('                         RESULTADOS MODALES\n')
fprintf('================================================================================\n\n')
fprintf('Modo      Omega        f[Hz]        T[s]        Gamma        Meff       Part.Acum\n')
fprintf('--------------------------------------------------------------------------------\n')

for i = 1:n
    fprintf('%3d   %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  %11.5f\n',...
        i,Omega(i),frecuencia(i),Te(i),Gamma(i),Meff(i),participacion_acumulada(i));
end

fprintf('\nPhi''*M*Phi =\n\n')
disp(Mmodal)

fprintf('\nPhi''*K*Phi =\n\n')
disp(Kmodal)

fprintf('\nError de reconstruccion modal = %.6e\n',error_superposicion)

%% ========================================================================
% RESULTADOS MAXIMOS
% =========================================================================

fprintf('\n')
fprintf('================================================================================\n')
fprintf('                    MAXIMOS DE LA RESPUESTA HISTORIA-TIEMPO\n')
fprintf('================================================================================\n\n')
fprintf('GDL      umax        vmax       a_rel_max    a_abs_max      Fsmax\n')
fprintf('--------------------------------------------------------------------------\n')

for i = 1:size(m,1)
    fprintf('%3d  %11.4e  %11.4e  %11.4e  %11.4e  %11.4e\n',...
        i,umax(i),vmax(i),arelmax(i),aabsmax(i),Fsmax(i));
end

%% ========================================================================
% GRAFICAS
% =========================================================================

figure
plot(t,ag,'k')
grid on
xlabel('Tiempo [s]')
ylabel('a_g')
title('Acelerograma')

figure
plot(t,u')
grid on
xlabel('Tiempo [s]')
ylabel('Desplazamiento relativo')
title('Historia de desplazamientos')

figure
plot(t,v')
grid on
xlabel('Tiempo [s]')
ylabel('Velocidad relativa')
title('Historia de velocidades')

figure
plot(t,a_abs')
grid on
xlabel('Tiempo [s]')
ylabel('Aceleracion absoluta')
title('Historia de aceleraciones absolutas')

figure
plot(t,Fs')
grid on
xlabel('Tiempo [s]')
ylabel('Fuerza elastica')
title('Historia de fuerzas elasticas')
