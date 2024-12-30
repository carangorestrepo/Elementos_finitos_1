function deformada(d,ang,esc_def,esc_faxial,esc_V,esc_M,x1,y1,x2,y2,M,V,fax,L)
%syms x       
         %deformada(Ava,Avb,Avc,Ixa,Ixb,Ixc,la,lb,w,E,G,ks,d,L,ang,esc_def,esc_faxial,esc_V,esc_M,x1,y1,x2,y2,qe)

X1 = 1; Y1 = 2; M1 = 3; X2 = 4; Y2 = 5; M2   = 6;
daA = d(1);
da = d(2);
ga = d(3);
dbA = d(4);
db = d(5);
gb = d(6);

u1 = d(1);
v1 = d(2);
t1 = d(3);
u2 = d(4);
v2 = d(5);
t2 = d(6);

%v=double(subs(v,x,xs));
%u=double(subs(u,x,xs));

%V=double(subs(V,x,xs));
%M=double(subs(M,x,xs));
%% Dibujar de deformada
figure(2)
xyv=[0,da;
     L,db]; %deformaciion vertical

xyu=[0,daA;
     L,dbA];%deformaciion horizontal

xs=[0,L];
s =[0,L];
% rotacion de la solucion antes de dibujar
T   = [ cos(ang)  -sin(ang)    % matriz de rotacion
        sin(ang)   cos(ang) ];

pos = T*[ xs + esc_def*xyu(:,2)';
           esc_def*xyv(:,2)'];
xx = pos(1,:) + x1;
yy = pos(2,:) + y1;

plot([x1 x2], [y1 y2], 'b-', xx, yy, 'r-','LineWidth',2),hold on

%% Dibujar los diagramas de fuerza axial 
figure(3)
axial=fax;
pos = T*[ s; esc_faxial*axial ]; % escalamiento del diagrama

ss = pos(1,:) + x1;
aa = pos(2,:) + y1;

plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 aa y2], 'r-','LineWidth',2),hold on;
text(ss(1),   aa(1),   num2str(axial(1,1)));
text(ss(end), aa(end), num2str(axial(1,end)));

%% Dibujar los diagramas de fuerza cortante

pos = T*[ s; esc_V*V ]; % escalamiento del diagrama
ss = pos(1,:) + x1;
vv = pos(2,:) + y1;
figure(4)
plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 vv y2], 'r-','LineWidth',2),hold on;
text(ss(1),   vv(1),   num2str(V(1,1)));
text(ss(end), vv(end), num2str(V(1,end)));
%% Dibujar los diagramas de momento flector


pos = T*[ s; -esc_M*M ]; % escalamiento del diagrama
ss = pos(1,:) + x1;
mm = pos(2,:) + y1;
figure(5)
plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 mm y2], 'r-','LineWidth',2),hold on;
text(ss(1),   mm(1),   num2str(M(1,1)));
text(ss(end), mm(end), num2str(M(1,end)));
[minM,idminM] = min(M); text(ss(idminM), mm(idminM), num2str(minM));
[maxM,idmaxM] = max(M); text(ss(idmaxM), mm(idmaxM), num2str(maxM));


  