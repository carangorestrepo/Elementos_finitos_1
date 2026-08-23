%% MAIN_ZAPATA_DKMT
clc
clear
close all
format long g

%% 1. GEOMETRIA
geo.Lx = 1.50;
geo.Ly = 1.70;
geo.hz = 0.40;

geo.lcx = 0.35;
geo.lcy = 0.40;

geo.cx0 = geo.Lx/2-geo.lcx/2;
geo.cy0 = geo.Ly/2-geo.lcy/2;

geo.pd = 2.50;
geo.hped = geo.pd-geo.hz;

%% 2. MATERIAL Y SUELO
mat.E = 24870062.3;
mat.nu = 0.25;
mat.kappa = 5/6;

suelo.gamma = 18.0;
suelo.ks = 6000.0;

concreto.gamma = 24.0;

%% 3. ACCIONES
acciones.P = -120.0;
acciones.Mx = 15.0;
acciones.My = 20.0;

acciones.P_incluye_pedestal = false;
acciones.P_incluye_zapata = false;

% Convencion:
% Mx = integral q*(y-yc)dA
% My = integral q*(x-xc)dA
acciones.signo_Mx = 1;
acciones.signo_My = 1;

%% 4. MALLA
malla.h_fuera = 0.06;
malla.h_huella = 0.03;
malla.diagonal_alternada = true;

[xnod,LaG,in_col,info_malla] = generar_malla_zapata(geo,malla);

fprintf('\nMALLA\n');
fprintf('Numero de nodos       = %d\n',size(xnod,1));
fprintf('Numero de elementos   = %d\n',size(LaG,1));
fprintf('Elementos en huella   = %d\n',sum(in_col));
fprintf('Area mallada          = %.8f m^2\n',info_malla.area_total);
fprintf('Area teorica zapata   = %.8f m^2\n',geo.Lx*geo.Ly);
fprintf('Area mallada huella   = %.8f m^2\n',info_malla.area_huella);
fprintf('Area teorica huella   = %.8f m^2\n',geo.lcx*geo.lcy);

%% 5. CARGAS
cargas = definir_cargas_zapata(geo,acciones,suelo,concreto);

fprintf('\nPRESIONES UNIFORMES\n');
fprintf('q_zapata   = %12.6f kN/m^2\n',cargas.q_zapata);
fprintf('q_suelo    = %12.6f kN/m^2\n',cargas.q_suelo);
fprintf('q_pedestal = %12.6f kN/m^2\n',cargas.q_pedestal);
fprintf('q_fuera    = %12.6f kN/m^2\n',cargas.q_fuera);

%% 6. OPCIONES
opciones.orden_gauss = 4;
opciones.gdl_restringidos = [];
opciones.escala_deformada = 100;
opciones.graficar = true;
opciones.tolerancia_equilibrio = 1e-6;

%% 7. SOLUCION
resultado = resolver_zapata_DKMT_modular( ...
    xnod,LaG,in_col,geo,mat,suelo,cargas,opciones);

%% 8. POSTPROCESO
postproceso_zapata(xnod,LaG,in_col,geo,suelo,resultado,opciones);
