function [x,V,M,w]=deformada(a,q,EIz,Acz,L,puntos_graficas)
%txi =a(1);
tyi =a(2);
wzi =a(3);
%txf =a(4);
tyf =a(5);
wzf =a(6);
q1=q;
q2=q;
Ac=Acz;
EI=EIz;

v1=a(3);
v2=a(6);
t1=a(2);
t2=a(5);

x=linspace(0,L,puntos_graficas);
%V =qZ*x - (24*Acz*EIz*wzf - 24*Acz*EIz*wzi + Acz*L^4*qZ + 12*EIz*L^2*qZ - 12*Acz*EIz*L*tyf - 12*Acz*EIz*L*tyi)/(2*L*(Acz*L^2 + 12*EIz));
%M=(qZ*x.^2)/2 - ((24*Acz*EIz*wzf - 24*Acz*EIz*wzi + Acz*L^4*qZ + 12*EIz*L^2*qZ - 12*Acz*EIz*L*tyf - 12*Acz*EIz*L*tyi)*x)/(2*L*(Acz*L^2 + 12*EIz)) + (144*EIz^2*tyf - 144*EIz^2*tyi + Acz*L^5*qZ + 12*EIz*L^3*qZ - 24*Acz*EIz*L^2*tyf - 48*Acz*EIz*L^2*tyi + 72*Acz*EIz*L*wzf - 72*Acz*EIz*L*wzi)/(12*L*(Acz*L^2 + 12*EIz));
%w=(qZ.*x.^4)/(24*EIz) - ((24.*Acz.*EIz.*wzf - 24.*Acz.*EIz.*wzi + Acz.*L^4.*qZ + 12.*EIz.*L^2.*qZ - 12.*Acz.*EIz.*L.*tyf - 12.*Acz.*EIz.*L.*tyi).*x.^3)/(12.*EIz.*L.*(Acz.*L^2 + 12.*EIz)) - ((144.*Acz.*EIz^2.*tyi - 144.*Acz.*EIz^2.*tyf - Acz^2.*L^5.*qZ + 144.*EIz^2.*L.*qZ - 72.*Acz^2.*EIz.*L.*wzf + 72.*Acz^2.*EIz.*L.*wzi + 24.*Acz^2.*EIz.*L^2.*tyf + 48.*Acz^2.*EIz.*L^2.*tyi).*x.^2)/(24.*Acz.*EIz.*L.*(Acz.*L^2 + 12.*EIz)) + ((2.*Acz^2.*L^3.*tyi + 24.*Acz.*EIz.*wzf - 24.*Acz.*EIz.*wzi + Acz.*L^4.*qZ + 12.*EIz.*L^2.*qZ - 12.*Acz.*EIz.*L.*tyf + 12.*Acz.*EIz.*L.*tyi).*x)/(2.*Acz.*L.*(Acz.*L^2 + 12.*EIz)) + wzi;

V =q1.*x - (240.*Ac.*EI.*v2 - 240.*Ac.*EI.*v1 + 7.*Ac.*L.^4.*q1 + 3.*Ac.*L.^4.*q2 + 80.*EI.*L.^2.*q1 + 40.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t1 - 120.*Ac.*EI.*L.*t2)/(20.*L.*(Ac.*L.^2 + 12.*EI)) - (x.^2.*(q1 - q2))/(2.*L);
M =(720.*EI.^2.*t2 - 720.*EI.^2.*t1 + 3.*Ac.*L.^5.*q1 + 2.*Ac.*L.^5.*q2 + 30.*EI.*L.^3.*q1 + 30.*EI.*L.^3.*q2 - 240.*Ac.*EI.*L.^2.*t1 - 120.*Ac.*EI.*L.^2.*t2 - 360.*Ac.*EI.*L.*v1 + 360.*Ac.*EI.*L.*v2)/(60.*L.*(Ac.*L.^2 + 12.*EI)) - (q1.*x.^2)/2 - x.*(((q1 - q2).*x.^2)/(2.*L) - q1.*x + (240.*Ac.*EI.*v2 - 240.*Ac.*EI.*v1 + 7.*Ac.*L.^4.*q1 + 3.*Ac.*L.^4.*q2 + 80.*EI.*L.^2.*q1 + 40.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t1 - 120.*Ac.*EI.*L.*t2)/(20.*L.*(Ac.*L.^2 + 12.*EI))) + (x.^3.*(q1 - q2))/(3.*L); 
w =(((q1 - q2).*x.^5)/(30.*L) - (q1.*x.^4)/24 + (EI.*(720.*EI.^2.*t2 - 720.*EI.^2.*t1 + 60.*Ac.^2.*L.^3.*v1 + 3.*Ac.*L.^5.*q1 + 2.*Ac.*L.^5.*q2 + 30.*EI.*L.^3.*q1 + 30.*EI.*L.^3.*q2 - 240.*Ac.*EI.*L.^2.*t1 - 120.*Ac.*EI.*L.^2.*t2 + 360.*Ac.*EI.*L.*v1 + 360.*Ac.*EI.*L.*v2))/(60.*Ac.*L.*(Ac.*L.^2 + 12.*EI)))/EI - (1/Ac - x.^2/(2.*EI)).*(((q1 - q2).*x.^3)/(3.*L) - (q1.*x.^2)/2 + (720.*EI.^2.*t2 - 720.*EI.^2.*t1 + 3.*Ac.*L.^5.*q1 + 2.*Ac.*L.^5.*q2 + 30.*EI.*L.^3.*q1 + 30.*EI.*L.^3.*q2 - 240.*Ac.*EI.*L.^2.*t1 - 120.*Ac.*EI.*L.^2.*t2 - 360.*Ac.*EI.*L.*v1 + 360.*Ac.*EI.*L.*v2)/(60.*L.*(Ac.*L.^2 + 12.*EI))) + (- x.^3/(6.*EI) + x/Ac).*(((q1 - q2).*x.^2)/(2.*L) - q1.*x + (240.*Ac.*EI.*v2 - 240.*Ac.*EI.*v1 + 7.*Ac.*L.^4.*q1 + 3.*Ac.*L.^4.*q2 + 80.*EI.*L.^2.*q1 + 40.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t1 - 120.*Ac.*EI.*L.*t2)/(20.*L.*(Ac.*L.^2 + 12.*EI))) + (x.*(- ((q1 - q2).*x.^4)/(8.*L) + (q1.*x.^3)/6 + EI.*t1))/EI;
