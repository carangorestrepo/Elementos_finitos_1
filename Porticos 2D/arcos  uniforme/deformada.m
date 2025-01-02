function [xx,yy,ssA,aaA,ssV,vvV,ssM,mm]=deformada(d,ang,ang2,esc_def,esc_faxial,esc_V,esc_M,x1,y1,x2,y2,M,V,fax,L,e,elementos)
%syms x       
         %deformada(Ava,Avb,Avc,Ixa,Ixb,Ixc,la,lb,w,E,G,ks,d,L,ang,esc_def,esc_faxial,esc_V,esc_M,x1,y1,x2,y2,qe)

X1 = 1; Y1 = 2; M1 = 3; X2 = 4; Y2 = 5; M2   = 6;
daA = d(1);
da = d(2);
ga = d(3);
dbA = d(4);
db = d(5);

gb = d(6);
% rotacion de la solucion antes de dibujar
T   = [ cos(ang)  -sin(ang)    % matriz de rotacion
        sin(ang)   cos(ang) ];


%v=double(subs(v,x,xs));
%u=double(subs(u,x,xs));

%V=double(subs(V,x,xs));
%M=double(subs(M,x,xs));
%% Dibujar de deformada
%%{
figure(2)
xyv=[0,da;
     L,db]; %deformaciion vertical

xyu=[0,daA;
     L,dbA];%deformaciion horizontal

xs=[0,L];
s =[0,L];


pos = T*[ xs + esc_def*xyu(:,2)';
           esc_def*xyv(:,2)'];
xx = pos(1,:) + x1;
yy = pos(2,:) + y1;

%plot([x1 x2], [y1 y2], 'b-', xx, yy, 'r-','LineWidth',2),hold on
%}
figure(2)
plot([x1 x2], [y1 y2], 'b-'),hold on

%% Dibujar los diagramas de fuerza axial 
figure(3)
axial=fax;
pos = T*[ s; esc_faxial*axial ]; % escalamiento del diagrama

ssA = pos(1,:) + x1;
aaA = pos(2,:) + y1;

plot([x1 x2], [y1 y2], 'b-', [x1 ssA x2], [y1 aaA y2], 'r-','LineWidth',2),hold on;
text(ssA(1),   aaA(1),   num2str(axial(1,1)),'Rotation',ang*180/pi+90);
if e==elementos
    text(ssA(end), aaA(end), num2str(axial(1,end)),'Rotation',ang*180/pi+90);
end

%% Dibujar los diagramas de fuerza cortante

pos = T*[ s; esc_V*V ]; % escalamiento del diagrama
ssV = pos(1,:) + x1;
vvV = pos(2,:) + y1;
figure(4)
plot([x1 x2], [y1 y2], 'b-', [x1 ssV x2], [y1 vvV y2], 'r-','LineWidth',2),hold on;
text(ssV(1),   vvV(1),   num2str(V(1,1)),'Rotation',ang*180/pi+90);
if e==elementos
    text(ssV(end), vvV(end), num2str(V(1,end)),'Rotation',ang*180/pi+90);
end
%% Dibujar los diagramas de momento flector


pos = T*[ s; -esc_M*M ]; % escalamiento del diagrama
ssM = pos(1,:) + x1;
mm = pos(2,:) + y1;
figure(5)
plot([x1 x2], [y1 y2], 'b-', [x1 ssM x2], [y1 mm y2], 'r-','LineWidth',2),hold on;
text(ssM(1),   mm(1),   num2str(M(1,1)),'Rotation',ang*180/pi+90);
if e==elementos
    text(ssM(end), mm(end), num2str(M(1,end)),'Rotation',ang*180/pi+90);
end
%[minM,idminM] = min(M); text(ssM(idminM), mm(idminM), num2str(minM));
%[maxM,idmaxM] = max(M); 
%text(ssM(idmaxM), mm(idmaxM), num2str(maxM),'Rotation',ang*180/pi+90);


  