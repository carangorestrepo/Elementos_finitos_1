function deformada(w,d,L,esc_def,esc_faxial,esc_V,esc_M3,esc_M2,x1,y1,x2,y2,z1,z2,~,T,EA,GJ,Acy,Acz,EIz,EIy)
%% deformaciones en nudos 
daA1 = d(1);
daV2 = d(2);
daV3 = d(3);
%daT4 = d(4);
gaM2 = d(5);
gaM3 = d(6);

dbA1 = d(7);
dbV2 = d(8);
dbV3 = d(9);
%dbT4 = d(10);
gbM2 = d(11);
gbM3 = d(12);
X=1;
Y=2;
Z=3;
datos=50;
x=linspace(0,L,datos);
%[~,~,v2,u2,V2,M3,axialu]=valores(L,Ava,Avb,Avc,Ixa,Iyze,Ixc,la,lb,w,E,G,ks,daA1,daV2,gaM3,dbA1,dbV2,gbM3,datos);
%[x,s,v3,~,V3,M2,~]=valores(L,Ava,Avb,Avc,Ixa,Ixze,Ixc,la,lb,0,E,G,ks,daA1,daV3,gaM2,dbA1,dbV3,gbM2,datos); 

[V2,M3,v2,u2,axialu]=cor_mome_def(w(1),w(2),Acz,EIz,L,gaM3,gbM3,daV2,dbV2,EA,0,0,daA1,dbA1,x);

[V3,M2,v3,~,~]=cor_mome_def(0,0,Acy,EIy,L,gaM2,gbM2,daV3,dbV3,EA,0,0,daA1,dbA1,x);

%% Dibujar de deformada
figure(1)
s=[0,x,L];
u2=[daA1,u2,dbA1];
v2=[daV2,v2,dbV2];
v3=[daV3,v3,dbV3];

pos = T*[ s + esc_def*u2; esc_def*v2; esc_def*v3 ];
%pos = T*[ s + esc_def*u2; esc_def*v3; esc_def*v2 ];
xx = pos(X,:)+x1;
yy = pos(Z,:)+y1;
zz = pos(Y,:)+z1;
plot3(xx,yy,zz,'r-','LineWidth',2),hold on;
hold on; title('Deformada exagerada');    xlabel('x, m'); ylabel('y, m'); xlabel('x, m'); zlabel('z, m'); axis equal
%figure(3); hold on; title('Fuerza axial [kN]');      xlabel('x, m'); ylabel('y, m'); axis equal
%figure(4); hold on; title('Fuerza cortante [kN]');   xlabel('x, m'); ylabel('y, m'); axis equal
grid on
axis equal
grid minor

%% graficas diagrmas de cortante V2
figure(2)
%s=[0,s,L];
%s=[x1,s,x2];
%V2=[x1,V2,x2];
%v3=[daV3,v3,dbV3];
ceros=zeros(1,datos);
pos = T*[ x ; esc_V*V2; ceros ];
%pos = T*[ s + esc_def*u2; esc_def*v3; esc_def*v2 ];
xx = [x1,pos(X,:)+x1,x2];
yy = [y1,pos(Z,:)+y1,y2];
zz = [z1,pos(Y,:)+z1,z2];
xe=[x1,x2];
ye=[y1,y2];
ze=[z1,z2];
plot3(xx,yy,zz,'r-','LineWidth',2),hold on;
plot3(xe,ye,ze,'--k','LineWidth', 2),hold on;
hold on; title('Fuerza cortante VZ [kN]');   xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); axis equal
grid on
axis equal
grid minor

%% graficas diagrmas de cortante V3
figure(3)
%s=[0,s,L];
%s=[x1,s,x2];
%V2=[x1,V2,x2];
%v3=[daV3,v3,dbV3];
ceros=zeros(1,datos);
pos = T*[ x ; ceros;esc_V*V3 ];
%pos = T*[ s + esc_def*u2; esc_def*v3; esc_def*v2 ];
xx = [x1,pos(X,:)+x1,x2];
yy = [y1,pos(Z,:)+y1,y2];
zz = [z1,pos(Y,:)+z1,z2];
xe=[x1,x2];
ye=[y1,y2];
ze=[z1,z2];
plot3(xx,yy,zz,'r-','LineWidth',2),hold on;
plot3(xe,ye,ze,'--k','LineWidth', 2),hold on;
hold on; title('Fuerza cortante VZ [kN]');   xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); axis equal
grid on
axis equal
grid minor

%% graficas diagrmas de cortante M3
figure(4)
%s=[0,s,L];
%s=[x1,s,x2];
%V2=[x1,V2,x2];
%v3=[daV3,v3,dbV3];
ceros=zeros(1,datos);
pos = T*[ x ; -esc_M3*M3; ceros ];
%pos = T*[ s + esc_def*u2; esc_def*v3; esc_def*v2 ];
xx = [x1,pos(X,:)+x1,x2];
yy = [y1,pos(Z,:)+y1,y2];
zz = [z1,pos(Y,:)+z1,z2];
xe=[x1,x2];
ye=[y1,y2];
ze=[z1,z2];
plot3(xx,yy,zz,'r-','LineWidth',2),hold on;
plot3(xe,ye,ze,'--k','LineWidth', 2),hold on;
hold on; title('Momento flector Mz[kN-m]'); xlabel('x, m'); ylabel('y, m');zlabel('z, m'); axis equal
grid on
axis equal
grid minor

%% graficas diagrmas de cortante M2
figure(5)
%s=[0,s,L];
%s=[x1,s,x2];
%V2=[x1,V2,x2];
%v3=[daV3,v3,dbV3];
ceros=zeros(1,datos);
pos = T*[ x ;ceros;-esc_M2*M2 ];
%pos = T*[ s + esc_def*u2; esc_def*v3; esc_def*v2 ];
xx = [x1,pos(X,:)+x1,x2];
yy = [y1,pos(Z,:)+y1,y2];
zz = [z1,pos(Y,:)+z1,z2];
xe=[x1,x2];
ye=[y1,y2];
ze=[z1,z2];
plot3(xx,yy,zz,'r-','LineWidth',2),hold on;
plot3(xe,ye,ze,'--k','LineWidth', 2),hold on;
hold on; title('Momento flector Mxy [kN-m]'); xlabel('x, m'); ylabel('y, m');zlabel('z, m'); axis equal
grid on
axis equal
grid minor

%% graficas diagrmas de cortante V2
figure(6)
%s=[0,s,L];
%s=[x1,s,x2];
%V2=[x1,V2,x2];
%v3=[daV3,v3,dbV3];
ceros=zeros(1,datos);
pos = T*[ x ; esc_faxial*axialu; ceros ];
%pos = T*[ s + esc_def*u2; esc_def*v3; esc_def*v2 ];
xx = [x1,pos(X,:)+x1,x2];
yy = [y1,pos(Z,:)+y1,y2];
zz = [z1,pos(Y,:)+z1,z2];
xe=[x1,x2];
ye=[y1,y2];
ze=[z1,z2];
plot3(xx,yy,zz,'r-','LineWidth',2),hold on;
plot3(xe,ye,ze,'--k','LineWidth', 2),hold on;
hold on; title('Fuerza axial [kN]');   xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); axis equal
grid on
axis equal
grid minor


