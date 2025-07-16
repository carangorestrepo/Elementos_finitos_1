function deformada(EA,EI,Ac,w,d,L,ang,esc_def,esc_faxial,esc_V,esc_M,x1,y1,x2,y2,puntos_graficas)
%syms x       
         %deformada(Ava,Avb,Avc,Ixa,Ixb,Ixc,la,lb,w,E,G,ks,d,L,ang,esc_def,esc_faxial,esc_V,esc_M,x1,y1,x2,y2,qe)

X1 = 1; Y1 = 2; M1 = 3; X2 = 4; Y2 = 5; M2   = 6;
daA = d(1);
da = d(2);
%ga = d(3);
dbA = d(4);
db = d(5);
%gb = d(6);

t1 = d(1);
u1 = d(2);
v1 = d(3);
t2 = d(4);
u2 = d(5);
v2 = d(6);

b1=0;
q1=-w(3);
b2=0;
q2=-w(6);
x=linspace(0,L,puntos_graficas);
s=linspace(0,L,puntos_graficas);

% momentos, cortantes,axiales y deformadas
V =q1.*x - (240.*Ac.*EI.*v2 - 240.*Ac.*EI.*v1 + 7.*Ac.*L.^4.*q1 + 3.*Ac.*L.^4.*q2 + 80.*EI.*L.^2.*q1 + 40.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t1 - 120.*Ac.*EI.*L.*t2)/(20.*L.*(Ac.*L.^2 + 12.*EI)) - (x.^2.*(q1 - q2))/(2.*L);
M =(720.*EI.^2.*t2 - 720.*EI.^2.*t1 + 3.*Ac.*L.^5.*q1 + 2.*Ac.*L.^5.*q2 + 30.*EI.*L.^3.*q1 + 30.*EI.*L.^3.*q2 - 240.*Ac.*EI.*L.^2.*t1 - 120.*Ac.*EI.*L.^2.*t2 - 360.*Ac.*EI.*L.*v1 + 360.*Ac.*EI.*L.*v2)/(60.*L.*(Ac.*L.^2 + 12.*EI)) - (q1.*x.^2)/2 - x.*(((q1 - q2).*x.^2)/(2.*L) - q1.*x + (240.*Ac.*EI.*v2 - 240.*Ac.*EI.*v1 + 7.*Ac.*L.^4.*q1 + 3.*Ac.*L.^4.*q2 + 80.*EI.*L.^2.*q1 + 40.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t1 - 120.*Ac.*EI.*L.*t2)/(20.*L.*(Ac.*L.^2 + 12.*EI))) + (x.^3.*(q1 - q2))/(3.*L);      
v =(((q1 - q2).*x.^5)/(30.*L) - (q1.*x.^4)/24 + (EI.*(720.*EI.^2.*t2 - 720.*EI.^2.*t1 + 60.*Ac.^2.*L.^3.*v1 + 3.*Ac.*L.^5.*q1 + 2.*Ac.*L.^5.*q2 + 30.*EI.*L.^3.*q1 + 30.*EI.*L.^3.*q2 - 240.*Ac.*EI.*L.^2.*t1 - 120.*Ac.*EI.*L.^2.*t2 + 360.*Ac.*EI.*L.*v1 + 360.*Ac.*EI.*L.*v2))/(60.*Ac.*L.*(Ac.*L.^2 + 12.*EI)))/EI - (1/Ac - x.^2/(2.*EI)).*(((q1 - q2).*x.^3)/(3.*L) - (q1.*x.^2)/2 + (720.*EI.^2.*t2 - 720.*EI.^2.*t1 + 3.*Ac.*L.^5.*q1 + 2.*Ac.*L.^5.*q2 + 30.*EI.*L.^3.*q1 + 30.*EI.*L.^3.*q2 - 240.*Ac.*EI.*L.^2.*t1 - 120.*Ac.*EI.*L.^2.*t2 - 360.*Ac.*EI.*L.*v1 + 360.*Ac.*EI.*L.*v2)/(60.*L.*(Ac.*L.^2 + 12.*EI))) + (- x.^3/(6.*EI) + x/Ac).*(((q1 - q2).*x.^2)/(2.*L) - q1.*x + (240.*Ac.*EI.*v2 - 240.*Ac.*EI.*v1 + 7.*Ac.*L.^4.*q1 + 3.*Ac.*L.^4.*q2 + 80.*EI.*L.^2.*q1 + 40.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t1 - 120.*Ac.*EI.*L.*t2)/(20.*L.*(Ac.*L.^2 + 12.*EI))) + (x.*(- ((q1 - q2).*x.^4)/(8.*L) + (q1.*x.^3)/6 + EI.*t1))/EI;
%u =(- ((b1 - b2).*x.^3)/(3.*L) + (b1.*x.^2)/2 + EA.*u1)/EA + (x.*(((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L)))/EA;
%fax =((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L);


figure(2)
xs=[0,x,L];
ys=zeros(1,puntos_graficas+2);
vs=[v(1),v,v(end)];
% rotacion de la solucion antes de dibujar
T   = [ cos(ang)  -sin(ang)    % matriz de rotacion
        sin(ang)   cos(ang) ];
%% Dibujar de deformada
pos = T*[ xs ; ys];
xx = pos(1,:) + x1;
yy = pos(2,:) + y1;
plot3([x1 x2], [y1 y2],[0 0], 'b-', xx, yy,vs*esc_def, 'r-','LineWidth',2),hold on,grid on

%% Dibujar los diagramas de fuerza torsional 
%figure(3)
%pos = T*[ s; esc_faxial*axial ]; % escalamiento del diagrama

%ss = pos(1,:) + x1;
%aa = pos(2,:) + y1;

%plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 aa y2], 'r-','LineWidth',2),hold on;
%text(ss(1),   aa(1),   num2str(axial(1,1)));
%text(ss(end), aa(end), num2str(axial(1,end)));

%% Dibujar los diagramas de fuerza cortante
ys=zeros(1,puntos_graficas);
pos = T*[ s; ys]; % escalamiento del diagrama
ss = pos(1,:) + x1;
vv = pos(2,:) + y1;
figure(3)
V=[V(1),V,V(end)];
plot3([x1 x2], [y1 y2],[0 0], 'b-', [x1 ss x2], [y1 vv y2],-esc_V*V , 'r-','LineWidth',2),hold on;

stem3([x1 ss x2],[y1 vv y2],-esc_V*V,'fill','b','markersize',1), hold on,grid on
text(ss(1),y1,-esc_V*V(1),num2str(V(1,1)));
text(ss(end),y2,-esc_V*V(end),num2str(V(1,end)));
%% Dibujar los diagramas de momento flector
pos = T*[ s;ys]; % escalamiento del diagrama
ss = pos(1,:) + x1;
mm = pos(2,:) + y1;
figure(4)
plot3([x1 x2], [y1 y2],[0 0], 'b-', [x1 ss x2], [y1 mm y2],-esc_M*[M(1),M,M(end)], 'r-','LineWidth',2),hold on;

stem3([x1 ss x2],[y1 mm y2],-esc_M*[M(1),M,M(end)],'fill','b','markersize',1), hold on,grid on
text(ss(1),y1,-esc_M*M(1),num2str(M(1,1)));
text(ss(end),y2,-esc_M*M(end), num2str(M(1,end)));
[minM,idminM] = min(M); 
text(ss(idminM),mm(idminM),-esc_M*M(idminM), num2str(minM));
[maxM,idmaxM] = max(M); 
text(ss(idmaxM),mm(idmaxM),-esc_M*(maxM), num2str(maxM));

  