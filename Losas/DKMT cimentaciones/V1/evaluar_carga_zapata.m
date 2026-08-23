function [qtotal,componentes] = evaluar_carga_zapata(cargas,x,y,en_huella)
% Evalua la presion total y sus componentes.

componentes.qPM=0;
componentes.q_zapata=cargas.q_zapata;
componentes.q_suelo=0;
componentes.q_pedestal=0;

if en_huella
    componentes.qPM = cargas.P/cargas.Ac ...
        + cargas.signo_My*cargas.My/cargas.Iy*(x-cargas.xc) ...
        + cargas.signo_Mx*cargas.Mx/cargas.Ix*(y-cargas.yc);

    componentes.q_pedestal=cargas.q_pedestal;
else
    componentes.q_suelo=cargas.q_suelo;
end

qtotal=componentes.qPM+componentes.q_zapata+ ...
       componentes.q_suelo+componentes.q_pedestal;
end
