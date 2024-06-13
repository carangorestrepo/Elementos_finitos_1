clc
clear
close all
%% 5.4 Lévy’s Solution for Rectangular Plates
a=2;
b=4;
E=210*1000;
t=0.05;
nu=0.3;
q=-10;
%%RRxRRy 
%%RRxEEy
%%EExEEy
ap='RRxEEy';

D = (E*t^3)/(12*(1 - nu^2));

deltax=0.05*2;
deltay=0.05*2;
n=2;
rho=2.3;

Nx=round(a/deltax,0);
dv=round(Nx/n,0);
Nx=dv*n+1;%% numero de nudos Y

Ny=round(b/deltay,0);
dv=round(Ny/n,0);
Ny=dv*n+1;%% numero de nudos Y

x=linspace(0,a,Nx);
y=linspace(-b/2,b/2,Ny);
[Xe,Ye]=meshgrid(x,y);%% grilla coordenadas nudos elemento
w=zeros(Nx,Ny);

m=1:2:100;
alfam=m*pi*b/(2*a);
t=tanh(alfam);
c=cosh(alfam);

Am=-2*(alfam.*tanh(alfam)+2)./(pi^5*m.^5.*cosh(alfam));
Bm=2./(pi^5*m.^5.*cosh(alfam));
s=0;
smx=0;
smy=0;
smxy=0;
sqx=0;
sqy=0;
for i=1:50
    s1=sinh(2*alfam(i)*Ye/b);
    c1=cosh(2*alfam(i)*Ye/b);
    %% apoyado en las 4 esquinas
    if strcmp(ap,'RRxRRy')==1
        s=s+1./m(i).^5.*(1-(m(i)*pi*b/(2*a).*tanh(alfam(i))+2)./(2*cosh(alfam(i))).*cosh(2*alfam(i)*Ye/b)...
        +1./(2*cosh(alfam(i))).*m(i)*pi.*Ye/a.*sinh(2*alfam(i)*Ye/b)).*sin(m(i).*pi.*Xe/a);
        smx=smx+D.*((pi.^2.*sin((pi.*Xe.*m(i))/a).*((Ye.*m(i).*pi.*sinh((2.*Ye.*alfam(i))/b))/(2.*a.*cosh(alfam(i))) - (cosh((2.*Ye.*alfam(i))/b).*(4.*a + pi.*b.*m(i).*tanh(alfam(i))))/(4.*a.*cosh(alfam(i))) + 1))/(a.^2.*m(i).^3) + (alfam(i).*nu.*sin((pi.*Xe.*m(i))/a).*(4.*a.*alfam(i).*cosh((2.*Ye.*alfam(i))/b) - 2.*b.*m(i).*pi.*cosh((2.*Ye.*alfam(i))/b) - 2.*Ye.*alfam(i).*m(i).*pi.*sinh((2.*Ye.*alfam(i))/b) + alfam(i).*b.*m(i).*pi.*cosh((2.*Ye.*alfam(i))/b).*tanh(alfam(i))))/(a.*b.^2.*m(i).^5.*cosh(alfam(i))));
        smy=smy+D.*((nu.*pi.^2.*sin((pi.*Xe.*m(i))/a).*((Ye.*m(i).*pi.*sinh((2.*Ye.*alfam(i))/b))/(2.*a.*cosh(alfam(i))) - (cosh((2.*Ye.*alfam(i))/b).*(4.*a + pi.*b.*m(i).*tanh(alfam(i))))/(4.*a.*cosh(alfam(i))) + 1))/(a.^2.*m(i).^3) + (alfam(i).*sin((pi.*Xe.*m(i))/a).*(4.*a.*alfam(i).*cosh((2.*Ye.*alfam(i))/b) - 2.*b.*m(i).*pi.*cosh((2.*Ye.*alfam(i))/b) - 2.*Ye.*alfam(i).*m(i).*pi.*sinh((2.*Ye.*alfam(i))/b) + alfam(i).*b.*m(i).*pi.*cosh((2.*Ye.*alfam(i))/b).*tanh(alfam(i))))/(a.*b.^2.*m(i).^5.*cosh(alfam(i))));
        smxy=smxy+(D.*pi.*cos((pi.*Xe.*m(i))/a).*(nu - 1).*(4.*a.*alfam(i).*sinh((2.*Ye.*alfam(i))/b) - b.*m(i).*pi.*sinh((2.*Ye.*alfam(i))/b) - 2.*Ye.*alfam(i).*m(i).*pi.*cosh((2.*Ye.*alfam(i))/b) + alfam(i).*b.*m(i).*pi.*sinh((2.*Ye.*alfam(i))/b).*tanh(alfam(i))))/(2.*a.^2.*b.*m(i).^4.*cosh(alfam(i)));
        sqx=sqx+D.*((pi.^3.*cos((pi.*Xe.*m(i))/a).*((Ye.*m(i).*pi.*sinh((2.*Ye.*alfam(i))/b))/(2.*a.*cosh(alfam(i))) - (cosh((2.*Ye.*alfam(i))/b).*((pi.*b.*m(i).*tanh(alfam(i)))/(2.*a) + 2))/(2.*cosh(alfam(i))) + 1))/(a.^3.*m(i).^2) - (pi.*cos((pi.*Xe.*m(i))/a).*((2.*alfam(i).*m(i).*pi.*cosh((2.*Ye.*alfam(i))/b))/(a.*b.*cosh(alfam(i))) - (2.*alfam(i).^2.*cosh((2.*Ye.*alfam(i))/b).*((pi.*b.*m(i).*tanh(alfam(i)))/(2.*a) + 2))/(b.^2.*cosh(alfam(i))) + (2.*Ye.*alfam(i).^2.*m(i).*pi.*sinh((2.*Ye.*alfam(i))/b))/(a.*b.^2.*cosh(alfam(i)))))/(a.*m(i).^4));
        sqy=sqy-D.*((sin((pi.*Xe.*m(i))/a).*((6.*alfam(i).^2.*m(i).*pi.*sinh((2.*Ye.*alfam(i))/b))/(a.*b.^2.*cosh(alfam(i))) - (4.*alfam(i).^3.*sinh((2.*Ye.*alfam(i))/b).*((pi.*b.*m(i).*tanh(alfam(i)))/(2.*a) + 2))/(b.^3.*cosh(alfam(i))) + (4.*Ye.*alfam(i).^3.*m(i).*pi.*cosh((2.*Ye.*alfam(i))/b))/(a.*b.^3.*cosh(alfam(i)))))/m(i).^5 - (pi.^2.*sin((pi.*Xe.*m(i))/a).*((m(i).*pi.*sinh((2.*Ye.*alfam(i))/b))/(2.*a.*cosh(alfam(i))) - (alfam(i).*sinh((2.*Ye.*alfam(i))/b).*((pi.*b.*m(i).*tanh(alfam(i)))/(2.*a) + 2))/(b.*cosh(alfam(i))) + (Ye.*alfam(i).*m(i).*pi.*cosh((2.*Ye.*alfam(i))/b))/(a.*b.*cosh(alfam(i)))))/(a.^2.*m(i).^3));
    elseif strcmp(ap,'RRxEEy')==1
   %% apoyado inf y sup empotrado izq y empotrado dere
        s=s+1./m(i).^5.*sin(m(i).*pi.*Xe./a).*(1-(alfam(i).*tanh(alfam(i))+2)./(2.*cosh(alfam(i))).*cosh(m(i).*pi.*Ye./a)+1./(2.*cosh(alfam(i))).*m(i).*pi.*Ye./a.*sinh(m(i).*pi.*Ye./a));
        smx=smx+(-(D.*pi.^2.*sin((pi.*Xe.*m(i))./a).*(2.*a.*cosh((pi.*Ye.*m(i))./a) - 2.*a.*cosh(alfam(i)) - Ye.*m(i).*pi.*sinh((pi.*Ye.*m(i))./a) + a.*alfam(i).*cosh((pi.*Ye.*m(i))./a).*tanh(alfam(i)) - a.*alfam(i).*nu.*cosh((pi.*Ye.*m(i))./a).*tanh(alfam(i)) + Ye.*m(i).*nu.*pi.*sinh((pi.*Ye.*m(i))./a)))./(2.*a.^3.*m(i).^3.*cosh(alfam(i))));
        smy=smy+(-(D.*pi.^2.*sin((pi.*Xe.*m(i))./a).*(2.*a.*nu.*cosh((pi.*Ye.*m(i))./a) - 2.*a.*nu.*cosh(alfam(i)) + Ye.*m(i).*pi.*sinh((pi.*Ye.*m(i))./a) - a.*alfam(i).*cosh((pi.*Ye.*m(i))./a).*tanh(alfam(i)) + a.*alfam(i).*nu.*cosh((pi.*Ye.*m(i))./a).*tanh(alfam(i)) - Ye.*m(i).*nu.*pi.*sinh((pi.*Ye.*m(i))./a)))./(2.*a.^3.*m(i).^3.*cosh(alfam(i))));
        smxy=smxy+((D.*pi.^2.*cos((pi.*Xe.*m(i))./a).*(nu - 1).*(a.*sinh((pi.*Ye.*m(i))./a).*cosh(alfam(i)) + a.*alfam(i).*sinh((pi.*Ye.*m(i))./a).*sinh(alfam(i)) - Ye.*m(i).*pi.*cosh((pi.*Ye.*m(i))./a).*cosh(alfam(i))))./(2.*a.^3.*m(i).^3.*cosh(alfam(i)).^2));
        sqx=sqx+(-D.*((pi.*cos((pi.*Xe.*m(i))./a).*((m(i).^2.*pi.^2.*cosh((pi.*Ye.*m(i))./a))./(a.^2.*cosh(alfam(i))) - (m(i).^2.*pi.^2.*cosh((pi.*Ye.*m(i))./a).*(alfam(i).*tanh(alfam(i)) + 2))./(2.*a.^2.*cosh(alfam(i))) + (Ye.*m(i).^3.*pi.^3.*sinh((pi.*Ye.*m(i))./a))./(2.*a.^3.*cosh(alfam(i)))))./(a.*m(i).^4) - (pi.^3.*cos((pi.*Xe.*m(i))./a).*((Ye.*m(i).*pi.*sinh((pi.*Ye.*m(i))./a))./(2.*a.*cosh(alfam(i))) - (cosh((pi.*Ye.*m(i))./a).*(alfam(i).*tanh(alfam(i)) + 2))./(2.*cosh(alfam(i))) + 1))./(a.^3.*m(i).^2)));
        sqy=sqy+(-D.*((sin((pi.*Xe.*m(i))./a).*((3.*m(i).^3.*pi.^3.*sinh((pi.*Ye.*m(i))./a))./(2.*a.^3.*cosh(alfam(i))) - (m(i).^3.*pi.^3.*sinh((pi.*Ye.*m(i))./a).*(alfam(i).*tanh(alfam(i)) + 2))./(2.*a.^3.*cosh(alfam(i))) + (Ye.*m(i).^4.*pi.^4.*cosh((pi.*Ye.*m(i))./a))./(2.*a.^4.*cosh(alfam(i)))))./m(i).^5 - (pi.^2.*sin((pi.*Xe.*m(i))./a).*((m(i).*pi.*sinh((pi.*Ye.*m(i))./a))./(2.*a.*cosh(alfam(i))) - (m(i).*pi.*sinh((pi.*Ye.*m(i))./a).*(alfam(i).*tanh(alfam(i)) + 2))./(2.*a.*cosh(alfam(i))) + (Ye.*m(i).^2.*pi.^2.*cosh((pi.*Ye.*m(i))./a))./(2.*a.^2.*cosh(alfam(i)))))./(a.^2.*m(i).^3)));    
    elseif strcmp(ap,'EExEEy')==1
        
    end
        
        
end 
figure
wa=q*4*a^4/(pi^5*D)*s;
Mx=q*4*a^4/(pi^5*D)*smx;
My=q*4*a^4/(pi^5*D)*smy;
Mxy=q*4*a^4/(pi^5*D)*smxy;
Qx=q*4*a^4/(pi^5*D)*sqx;
Qy=q*4*a^4/(pi^5*D)*sqy;

surf(Xe,Ye,wa)
axis equal tight
colormap jet
title('Deformada','FontSize',20)
figure
subplot(1,3,1)
contourf(Xe,Ye,Mx,12)
colorbar
axis equal tight
colormap jet
title('Momentos M_x (KN-m/m)', 'FontSize',20);
subplot(1,3,2)
contourf(Xe,Ye,My,12)
colorbar
axis equal tight
colormap jet
title('Momentos M_y (KN-m/m)', 'FontSize',20);
subplot(1,3,3)
contourf(Xe,Ye,Mxy,12)
colorbar
axis equal tight
colormap jet
title('Momentos M_xy (KN-m/m)', 'FontSize',20);

figure
subplot(1,3,1)
contourf(Xe,Ye,Qx,12)
colorbar
axis equal tight
colormap jet
title('Cortantes Q_x (KN/m)', 'FontSize',20);

subplot(1,3,2)
contourf(Xe,Ye,Qy,12)
colorbar
axis equal tight
colormap jet
title('Cortantes Q_y (KN/m)', 'FontSize',20);

subplot(1,3,3)
contourf(Xe,Ye,max(abs(Qy),abs(Qx)),12)
colorbar
axis equal tight
colormap jet
title('Cortantes Q_max (KN/m)', 'FontSize',20);




