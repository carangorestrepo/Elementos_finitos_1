function cargas = definir_cargas_zapata(geo,acciones,suelo,concreto)
% Convierte P, Mx y My en presion lineal sobre la huella y define pesos.

cargas.Ac = geo.lcx*geo.lcy;
cargas.xc = geo.cx0+geo.lcx/2;
cargas.yc = geo.cy0+geo.lcy/2;

cargas.Ix = geo.lcx*geo.lcy^3/12;
cargas.Iy = geo.lcy*geo.lcx^3/12;

cargas.P  = acciones.P;
cargas.Mx = acciones.Mx;
cargas.My = acciones.My;

cargas.signo_Mx = acciones.signo_Mx;
cargas.signo_My = acciones.signo_My;

cargas.q_zapata = -concreto.gamma*geo.hz;
if acciones.P_incluye_zapata
    cargas.q_zapata=0;
end

h_suelo=geo.pd-geo.hz;
cargas.q_suelo=-suelo.gamma*h_suelo;

cargas.q_pedestal=-concreto.gamma*geo.hped;
if acciones.P_incluye_pedestal
    cargas.q_pedestal=0;
end

cargas.q_fuera=cargas.q_zapata+cargas.q_suelo;

cargas.W_zapata=cargas.q_zapata*geo.Lx*geo.Ly;
cargas.W_suelo=cargas.q_suelo*(geo.Lx*geo.Ly-cargas.Ac);
cargas.W_pedestal=cargas.q_pedestal*cargas.Ac;

cargas.Fz_teorica=acciones.P+cargas.W_zapata+ ...
                   cargas.W_suelo+cargas.W_pedestal;
end
