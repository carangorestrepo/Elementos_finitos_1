% se definen algunas constantes que hacen el codigo mas legible
X = 1;
Y = 2;
G = 3;
XY = 12;
XYG=123;
g=9.8066502;      % aceleracion de la gravedad

%% escalas de Dibujo la estructura y su deformada
esc_def    = 7;          % escalamiento de la deformada
esc_faxial = 0.006;      % escalamiento del diagrama de axiales
esc_V      = 1;          % escalamiento del diagrama de cortantes
esc_M      = 0.9;        % escalamiento del diagrama de momentos
%fesct=max(max(xyc))*11/16*3;%% escala textos
%fesc=max(max(xyc))*0.4/16;%% escala apoyos

% Parámetros de entrada
qa = 13;% 13;           % Carga distribuida (N/m)
qb = 10; %10;           % Carga distribuida (N/m)
bv = 25/1000;           % Ancho de la viga (m)
hv = 50/1000;           % Altura de la viga (m)
I = bv * hv^3 / 12;     % Inercia de la viga (m^4)
E = 210000 * 1000;      % Módulo de elasticidad del material (Pa)
v = 0.3;                % Coeficiente de Poisson
Gm = E / (2 * (1 + v)); % Módulo de rigidez en corte (Pa)
EI = E * I;             % Producto de módulo de elasticidad y momento de inercia (N·m²)
Ac = bv*hv*5/6*Gm;
%% tipo de apoyo 
% retriccion,nudo
%%X horizontal, Y vertical, XY horozonal y vertical, XYG horizontal vertica y giro
tipo_apoyo=[XYG,1;
            XYG,45];

%% cordenadas de arco
%xn=[-3 ,3.2,5.6];
%yn=[2.5,5,3.5];
xn=[0 ,4.3,8.6];
yn=[0,3,0];
L=xn(3)-xn(1);
xa=xn(1);
xb=xn(2);
xc=xn(3);
ya=yn(1);
yb=yn(2);
yc=yn(3);

%% divido arco en tramos                    1  2  3  4  6 
datos=89;%% debe ser un numero impar mayor a 3, 5 ,7 ,9 11
xx=linspace(xa,xc,datos)';
%% ecuacion de arco en cordenadas no traladadas
yr =(xx*xa^2*yb - xx*xb^2*ya - xx.^2*xa*yb + xx.^2*xb*ya - xx*xa^2*yc + xx*xc^2*ya + xx.^2*xa*yc - xx.^2*xc*ya + xx*xb^2*yc - xx*xc^2*yb - xx.^2*xb*yc + xx.^2*xc*yb - xa*xb^2*yc + xa*xc^2*yb - xb*xc^2*ya + xa^2*xb*yc - xa^2*xc*yb + xb^2*xc*ya)/((xa - xb)*(xa - xc)*(xb - xc));

pen=((xa^2*yb - xb^2*ya - xa^2*yc + xc^2*ya + xb^2*yc - xc^2*yb - 2*xx*xa*yb + 2*xx*xb*ya + 2*xx*xa*yc - 2*xx*xc*ya - 2*xx*xb*yc + 2*xx*xc*yb)/((xa - xb)*(xa - xc)*(xb - xc)));
 
xyc=[xx,yr];

recorrido=4:2:datos;
sizerecorrido=size(recorrido,2);
elementos=sizerecorrido+1;% numeor de arcos en los que divido arco inicial
xe=zeros(elementos,3);
ye=zeros(elementos,3);
ang1=zeros(elementos,3);
xe(1,1:3)=xx(1:3,1)';
ye(1,1:3)=yr(1:3,1)';
ang1(1,1:3)=pen(1:3,1)';
i=2;
for e=1:sizerecorrido
   xe(i,:)=xx((recorrido(e)-1):(recorrido(e)+1),1)'; % cordenadas de un tramo de arco x  
   ye(i,:)=yr((recorrido(e)-1):(recorrido(e)+1),1)'; % cordenadas de un tramo de arco y
   ang1(i,:)=pen((recorrido(e)-1):(recorrido(e)+1),1)'; % cordenadas de un tramo de arco y
   i=i+1;
end
ninf=[(1:elementos)',(2:(elementos+1))'];
nudos=elementos+1;
%% numero grados de livertad
ngdl=nudos*3;
%% separo memoria
Kloce = cell(elementos,1);
T = cell(elementos,1);
ang = zeros(elementos,1); 
ang2 = zeros(elementos,2); 
Fe = cell(elementos,1);
GLe = zeros(elementos,6);
xxe = zeros(elementos,4);
q=zeros(elementos,2);
K = zeros(nudos*3);
M = zeros(nudos*3);
f = zeros(nudos*3,1);
ca = zeros(nudos*3,1);
napoyo=size(tipo_apoyo,1);
Myy = cell(elementos,1);
Le = zeros(elementos,1);
qe_glob = cell(elementos(1,1),1);
for e=1:elementos
    GLe(e,:) = [(ninf(e,X)*3-2):(ninf(e,X)*3),(ninf(e,Y)*3-2):(ninf(e,Y)*3)];
    xaa=xe(e,1);xbb=xe(e,2);xcc=xe(e,3);
    yaa=ye(e,1);ybb=ye(e,2);ycc=ye(e,3);
    Kloce{e}=k_rigidez(xaa,xbb,xcc,yaa,ybb,ycc,EI,Ac);
    qaa=(qb-qa)/(xc-xa)*(xaa-xa)+qa;
    qbb=(qb-qa)/(xc-xa)*(xcc-xa)+qa;
    q(e,:)=[qaa,qbb];
    Le(e)=hypot(xe(e,3)-xe(e,1),ye(e,3)-ye(e,1));
    ang(e) = atan2(ye(e,3)-ye(e,1), xe(e,3)-xe(e,1));
    ang2(e,:)=atan(ang1(e,[1,3]));
    
    ss = (xe(e,3)-xe(e,1))/Le(e);   
    cc = (ye(e,3)-ye(e,1))/Le(e);
    % matriz de transformacion de coordenadas para la barra e   
    T{e} = [ cc  ss  0  0  0  0
            -ss  cc  0  0  0  0
             0  0  1  0  0  0
             0  0  0  cc  ss  0
             0  0  0 -ss  cc  0
             0  0  0  0  0  1];
    Fe{e}=-fuerzas_nodales_equivalentes(qaa,qbb,xaa,xbb,xcc,yaa,ybb,ycc,EI,Ac);
    Fee=Fe{e};
    %% matriz de rigidez local en coordenadas globales
    %% calculo de  la matriz de masa concentrada
     Myy{e}=   [abs(Fee(2,1))/g,              0,              0,              0,              0,            0;
                            0,abs(Fee(2,1))/g,              0,              0,              0,            0;
                            0,              0,              0,              0,              0,            0;
                            0,              0,              0,abs(Fee(5,1))/g,              0,            0;
                            0,              0,              0,              0,abs(Fee(5,1))/g,            0;
                            0,              0,              0,              0,              0,           0];
    K(GLe(e,:),GLe(e,:)) = K(GLe(e,:),GLe(e,:)) + Kloce{e}; % ensambla Ke{e} en K global
    f(GLe(e,:))          = f(GLe(e,:))          + Fe{e}; % ensambla fe{e} en f global
    M(GLe(e,:),GLe(e,:)) = M(GLe(e,:),GLe(e,:)) + Myy{e}; % ensambla My{e} en M global
end

%% grados de libertad que no tienen ceros
smy=sum(abs(K),1);
[~,fy]=find(smy~=0);
% extraigo las submatrices y especifico las cantidades conocidas
for ap=1:napoyo
    if tipo_apoyo(ap,1)==XYG
        ca(tipo_apoyo(ap,2)*3-2:tipo_apoyo(ap,2)*3,1)=tipo_apoyo(ap,2)*3-2:tipo_apoyo(ap,2)*3;
    elseif tipo_apoyo(ap,1)==XY
        ca(tipo_apoyo(ap,2)*3-2:tipo_apoyo(ap,2)*3-1,1)=tipo_apoyo(ap,2)*3-2:tipo_apoyo(ap,2)*3-1;
    elseif tipo_apoyo(ap,1)==Y
        ca(tipo_apoyo(ap,2)*3-1,1)=tipo_apoyo(ap,2)*3-1;
    elseif tipo_apoyo(ap,1)==X 
        ca(tipo_apoyo(ap,2)*3-2,1)=tipo_apoyo(ap,2)*3-2;
    end

        %figure(1)
        %xr=xyc(tipo_apoyo(ap,2),1);
        %yr=xyc(tipo_apoyo(ap,2),2);
        %[x,y,xn,yn]=apoyos(xr,yr,fesc,tipo_apoyo(ap,1));
        %plot(x,y,'-r',xn,yn,'-r'), hold on;
        %fill(x,y,'r'),hold on
        %alpha(.5)

end
Kel=K(fy,fy);
fel=f(fy,1);
%%
% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
%| qd   |   | Kcc Kcd || ac |   | fd |  Recuerde que qc = 0
%|      | = |         ||    | - |    |
%| qc=0 |   | Kdc Kdd || ad |   | fc | 
c = setdiff(ca,0); % desplazamientos conocidos
apoyos1=size(c,1);
d = setdiff(fy, c);
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);
% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = zeros(apoyos1,1); % desplazamientos conocidos en contorno
%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);        % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
Mdd=M(d,d);
%MMexdd=MMex(d,d);
% armo los vectores de desplazamientos (a) y fuerzas (q)
a  = zeros(ngdl,1);   a(c) = ac;  a(d) = ad; % desplazamientos
q  = zeros(ngdl,1);   q(c) = qd; % fuerzas nodales de equilibrio
xxeA = zeros(elementos,4);
xxeV = zeros(elementos,4);
xxeM = zeros(elementos,4);

for e=1:elementos(1,1)   % para cada barra
    qe_glob{e} = (Kloce{e}*a(GLe(e,:)) - Fe{e});
    x1 = xe(e,1);  x2 = xe(e,3);
    y1 = ye(e,1);  y2 = ye(e,3);
    sol= qe_glob{e};
    M=-[sol(3),-sol(6)];
    Vv=[sol(2),-sol(5)];
    Aa=[sol(1),-sol(4)];
    s=sin(atan(ang1(e,[1,3])));
    c=cos(atan(ang1(e,[1,3])));
    A=(-Aa.*c-Vv.*s);
    V=-(-Vv.*c+Aa.*s);
    %T{e}*a(GLe(e,:))
    [xx,yy,ssA,aaA,ssV,vvV,ssM,mm]=deformada(a(GLe(e,:)),ang(e),ang2(e,:),esc_def,esc_faxial,esc_V,esc_M,x1,y1,x2,y2,M,V,A,Le(e),e,elementos(1,1));
    xxe(e,:)=[xx,yy];
    xxeA(e,:)=[ssA,aaA];
    xxeV(e,:)=[ssV,vvV];
    xxeM(e,:)=[ssM,mm];
    %qe_loc{e}=T{e}*qe_glob{e};
end
xA=[xxeA(:,1);xxeA(end,2)];
Aa=[xxeA(:,3);xxeA(end,4)];

xVv=[xxeV(:,1);xxeV(end,2)];
Vv=[xxeV(:,3);xxeV(end,4)];

xMm=[xxeM(:,1);xxeM(end,2)];
Mm=[xxeM(:,3);xxeM(end,4)];


figure(2)
plot([xxe(:,1);xxe(end,2)],[xxe(:,3);xxe(end,4)], 'r-','LineWidth',2)
%figure(5)
%plot([x1 x2], [y1 y2], 'b-', [x1 ssM x2], [y1 Mm y2], 'r-','LineWidth',2),hold on;
%text(xMm(1),   Mm(1),   num2str(Mm(1,1)));
%text(xMm(end), Mm(end), num2str(Mm(end,1)));
%[minM,idminM] = min(Mm); 
%text(xMm(idminM), Mm(idminM), num2str(minM));
%[maxM,idmaxM] = max(Mm); 
%text(xMm(idmaxM), Mm(idmaxM), num2str(maxM));


figure(2); hold on; title('Deformada exagerada');    xlabel('x, m'); ylabel('y, m'); axis equal
figure(3); hold on; title('Fuerza axial [kN]');      xlabel('x, m'); ylabel('y, m'); axis equal
figure(4); hold on; title('Fuerza cortante [kN]');   xlabel('x, m'); ylabel('y, m'); axis equal
figure(5); hold on; title('Momento flector [kN-m]'); xlabel('x, m'); ylabel('y, m'); axis equal
