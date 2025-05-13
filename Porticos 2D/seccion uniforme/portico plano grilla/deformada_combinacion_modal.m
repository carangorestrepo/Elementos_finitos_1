function deformada_combinacion_modal(EA,EI,Ac,w,d,L,ang,kw,P,lalb,esc_def,esc_faxial,esc_V,esc_M,x1,y1,x2,y2,tipo_conti,iteraciones,ite,puntos_graficas,n)
%syms x       
         %deformada(Ava,Avb,Avc,Ixa,Ixb,Ixc,la,lb,w,E,G,ks,d,L,ang,esc_def,esc_faxial,esc_V,esc_M,x1,y1,x2,y2,qe)

X1 = 1; Y1 = 2; M1 = 3; X2 = 4; Y2 = 5; M2   = 6;

b1=w{1,1}(1);
q1=w{1,1}(2);
b2=w{1,1}(4);
q2=w{1,1}(5);
xs=linspace(0,L,puntos_graficas);
s=linspace(0,L,puntos_graficas);
Vc=zeros(n,puntos_graficas);
Mc=zeros(n,puntos_graficas);
Faxc=zeros(n,puntos_graficas);
vc=zeros(n,puntos_graficas);
uc=zeros(n,puntos_graficas);

daA=zeros(1,n);
da=zeros(1,n);
dbA=zeros(1,n);
db=zeros(1,n);
for i=1:n
    daA(i) = d(1,i);
    da(i) = d(2,i);
    %ga = d(3,i);
    dbA(i) = d(4,i);
    db(i) = d(5,i);
    %gb = d(6,i);
    u1 = d(1,i);
    v1 = d(2,i);
    t1 = d(3,i);
    u2 = d(4,i);
    v2 = d(5,i);
    t2 = d(6,i);
    [~,~,~,V,M,fax,v,u]=matriz_rigidez(tipo_conti,EA,Ac,EI,L,b1,b2,q1,q2,1,v1,v2,t1,t2,u1,u2,xs,kw,P,lalb,puntos_graficas);
    Vc(i,:)=V;
    Mc(i,:)=M;
    Faxc(i,:)=fax;
    vc(i,:)=v;
    uc(i,:)=u;
end
M = -rssq(Mc,1);
V = rssq(Vc,1);
fax = rssq(Faxc,1);
v = rssq(vc,1).*sign(vc(1,:));
u = rssq(uc,1).*sign(uc(1,:));
daA = rssq(daA)*sign(d(1,1));
da=rssq(da)*sign(d(2,1));
dbA = rssq(dbA)*sign(d(4,1));
db = rssq(db)*sign(d(5,1));
if ite==iteraciones
    %v=double(subs(v,x,xs));
    %u=double(subs(u,x,xs));
    axial=fax;
    %V=double(subs(V,x,xs));
    %M=double(subs(M,x,xs));
    figure(2)
    xyv=[0,da;xs',v';L,db];
    xyu=[0,daA;xs',u';L,dbA];
    xs=[0,xs,L];
    % rotacion de la solucion antes de dibujar
    T   = [ cos(ang)  -sin(ang)    % matriz de rotacion
            sin(ang)   cos(ang) ];
    %% Dibujar de deformada
    pos = T*[ xs + esc_def*xyu(:,2)'; esc_def*xyv(:,2)'];
    xx = pos(1,:) + x1;
    yy = pos(2,:) + y1;
    plot([x1 x2], [y1 y2], 'b-', xx, yy, 'r-','LineWidth',2),hold on

    %% Dibujar los diagramas de fuerza axial 
    figure(3)
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
end

  