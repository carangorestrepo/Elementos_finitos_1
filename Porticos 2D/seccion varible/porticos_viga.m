function [u,v,M,V,P,x]=porticos_viga(EG,mat_elem,xyc,ninf,sec,elemento_carga,tipo_apoyo,tf,tw,bh,nudo_carga,aa)

% se definen algunas constantes que hacen el codigo mas legible
X = 1;
Y = 2;
G = 3;
XY = 12;
XYG=123;
%% porticos_calculos(EG,mat_elem,xyc,ninf,sec,sec_ele,elemento_carga,nudo_carga,tipo_apoyo,tipo_conti,X,Y)
elementos = size(ninf,1);
%% numero nudos
nudos = size(xyc,1);
%Grafico nudos
%figure(1)
%plot(xyc(:,1),xyc(:,2),'o','MarkerSize',6,...
%'MarkerEdgeColor','blue',...
%'MarkerFaceColor','blue'), hold on
%xlabel('x(m)')
%ylabel('y(m)')
%h = text(xyc(:,1)+0.1,xyc(:,2)+0.1, num2str((1:nudos)'),'FontSize',fesct);
%set(h, 'Color','blue');
%% numero grados de livertad
ngdl=nudos*3;
%% separo memoria
xye = zeros(elementos,4);
GLe = zeros(elementos,6);
Le = zeros(elementos,1);
Ee = zeros(elementos,1);
Ge = zeros(elementos,1);
%Ace = zeros(elementos,1);
T = cell(elementos,1);
wex = zeros(elementos,2);
wey = zeros(elementos,2);
Wex = cell(elementos,1); 
Wey = cell(elementos,1); 
Cex = cell(elementos,1); 
Cey = cell(elementos,1); 
Ce = cell(elementos,1); 
ang = zeros(elementos,1); 
Fexx = cell(elementos,1);
Feyy = cell(elementos,1);
Fe = cell(elementos,1);
Kea = cell(elementos,1);
Kloce = cell(elementos,1);
Meex = cell(elementos,1);
Meey = cell(elementos,1);
K = zeros(nudos*3);
napoyo=size(tipo_apoyo,1);
%apoyo=zeros(napoyo,1);
ca = zeros(nudos*3,1);
f = zeros(nudos*3,1);

qe_glob = cell(elementos(1,1),1);
qe_loc  = cell(elementos(1,1),1);

size_elemento_carga=size(elemento_carga,1);

%% vector de cargas distribuidas
for i=1:size_elemento_carga
    if elemento_carga(i,4)==X
        wex(elemento_carga(i,3),1:2)=wex(elemento_carga(i,3),1:2)+elemento_carga(i,1:2);
    end
    if elemento_carga(i,4)==Y
        wey(elemento_carga(i,3),1:2)=wey(elemento_carga(i,3),1:2)+elemento_carga(i,1:2);
    end  
end
%% vector de cargas puntuales
sizenudo_carga=size(nudo_carga,1);
for i=1:sizenudo_carga
    if  nudo_carga(i,3)==X
        f(nudo_carga(i,2)*3-2,1)=nudo_carga(i,1);
    elseif  nudo_carga(i,3)==Y 
        f(nudo_carga(i,2)*3-1,1)=nudo_carga(i,1);
    elseif nudo_carga(i,3)==G 
        f(nudo_carga(i,2)*3,1)=nudo_carga(i,1);
    end
end

%% Asignacion de 

for e=1:elementos % para cada barra
    % saco los 4 cordenadas de la barra e
    xye(e,:) = [xyc(ninf(e,1),1),xyc(ninf(e,1),2),xyc(ninf(e,2),1),xyc(ninf(e,2),2)];
    GLe(e,:) = [(ninf(e,X)*3-2):(ninf(e,X)*3),(ninf(e,Y)*3-2):(ninf(e,Y)*3)];
    L = hypot(xye(e,3)-xye(e,1),xye(e,4)-xye(e,2));
    Le(e)=L;
    ang(e) = atan2(xye(e,4)-xye(e,2), xye(e,3)-xye(e,1));

    %% materiales y secciones elementos
    Ee(e) = EG(mat_elem(e),1); %% E
    Ge(e) = EG(mat_elem(e),2);%% G
    wx = wex(e,:);
    wy = wey(e,:);
    nbe=bh(e,5);
    nhe=bh(e,6);
    %% seco y coseno direcctor
    s = (xye(e,4)-xye(e,2))/L;   
    c = (xye(e,3)-xye(e,1))/L;
    nqe=elemento_carga(e,5);
    naxie=elemento_carga(e,6);
    % matriz de transformacion de coordenadas para la barra e   
    T{e} = [ c  s  0  0  0  0
           -s  c  0  0  0  0
            0  0  1  0  0  0
            0  0  0  c  s  0
            0  0  0 -s  c  0
            0  0  0  0  0  1];
    % matriz de rigidez nudos rigidos local expresada en el sistema de coordenadas locales
    % para la barra e
    %[Kloc,~,Mex]=matriz_rigidez(tipo_conti(e,:),EA,Ac,EI,L,wex(e,1),wex(e,2),wex(e,1),wex(e,2),0,0,0,0,0,0,0,0,kw(e),P(e),lalb(e,:),puntos_graficas);
    %[~,~,Mey]   =matriz_rigidez(tipo_conti(e,:),EA,Ac,EI,L,wey(e,1),wey(e,2),wey(e,1),wey(e,2),0,0,0,0,0,0,0,0,kw(e),P(e),lalb(e,:),puntos_graficas);
    
    %Kloc=Ke(nqe,sec,bh(e,1),bh(e,2),nbe,bh(e,3),bh(e,4),nhe,Ge(e),Ee(e),tf,tw,L);
    Kloc=Ke_1(bh(e,1),bh(e,2),nbe,bh(e,3),bh(e,4),nhe,Ee(e),Ge(e),L,tw,tf,sec);
     
    Kloce{e}=Kloc;
    
    [MEEya]=Ax_Ix_As2(bh(e,1),bh(e,2),nbe,bh(e,3),bh(e,4),nhe,wey(e,1),wey(e,2),nqe,wey(e,1),wey(e,2),naxie,0,L,Ee(e),Ge(e),L,tw,tf,sec,0,0);
    [MEExa]=Ax_Ix_As2(bh(e,1),bh(e,2),nbe,bh(e,3),bh(e,4),nhe,wex(e,1),wex(e,2),nqe,wex(e,1),wex(e,2),naxie,0,L,Ee(e),Ge(e),L,tw,tf,sec,0,0);
    
    Fey=MEEya;
    Fex=MEExa;
    %[MEEy]=funcion(wey(e,1),wey(e,2),nqe,sec,bh(e,1),bh(e,2),nbe,bh(e,3),bh(e,4),nhe,Ge(e),Ee(e),tf,tw,L,0,0,0,L);
    %[MEEx]=funcion(wex(e,1),wex(e,2),nqe,sec,bh(e,1),bh(e,2),nbe,bh(e,3),bh(e,4),nhe,Ge(e),Ee(e),tf,tw,L,0,0,0,L);
    %[Ry]=funcion1(tf,tw,bh(e,1),bh(e,2),nbe,bh(e,3),bh(e,4),nhe,wey(e,1),wey(e,2),naxie,L,Ee(e),Ge(e),0,0,L,sec);
    %[Rx]=funcion1(tf,tw,bh(e,1),bh(e,2),nbe,bh(e,3),bh(e,4),nhe,wex(e,1),wex(e,2),naxie,L,Ee(e),Ge(e),0,0,L,sec);
    
    %Fey=[Ry(1);MEEy(1:2);Ry(2);MEEy(3:4)];
    
    %Fex=[Rx(1);MEEx(1:2);Rx(2);MEEx(3:4)];
    %[Kloc,FE,Me,V,M,fax,v,u]=fuerzas_empotramiento_bwinckler(EA,Ac,EI,L,we(e,1),we(e,2),0,0,0,0,0,0,0)
    %% direccion de carga distribuida
    %if we(e,4)==X   %DIRECCÓN HORIZONTAL
    direcargax=[1;0;0;1;0;0];
    Wex{e}=T{e}*[wex(e,1);0;0;wex(e,2);0;0];
    Cex{e}=[wx(1)*c;-wx(1)*s;-wx(1)*s;wx(2)*c;-wx(2)*s;-wx(2)*s];
    %if we(e,4)==Y     %DIRECCÓN VERTICAL
    direcargay=[0;1;0;0;1;0]; 
    Wey{e}=T{e}*[0;wey(e,1);0;0;wey(e,2);0];
    Cey{e}=[wy(1)*s;wy(1)*c;wy(1)*c;wy(2)*s;wy(2)*c;wy(2)*c];
    Ce{e}=Cey{e}+Cex{e};
    %end
    %elseif we(e,4)==0
    %    direcarga=[0;0;0;0;0;0]; 
    %    We{e}=[0;0;0;0;0;0]; 
    %    Ce{e}=[0;0;0;0;0;0]; 
    %end
    %% matriz de rigidez local en coordenadas globales
        % matriz de mometos de empotramiento
    Mey=[Fey(1)  ,0     ,0,0     ,0     ,0;
             0   ,Fey(2),0,0     ,0     ,0;
             0   ,Fey(3),0,0     ,0     ,0;
             0   ,0     ,0,Fey(4),0     ,0;
             0   ,0     ,0,0     ,Fey(5),0;
             0   ,0     ,0,0     ,Fey(6),0];
     
	Mex=[Fex(1)  ,0     ,0,0     ,0     ,0;
             0   ,Fex(2),0,0     ,0     ,0;
             0   ,Fex(3),0,0     ,0     ,0;
             0   ,0     ,0,Fex(4),0     ,0;
             0   ,0     ,0,0     ,Fex(5),0;
             0   ,0     ,0,0     ,Fex(6),0];
    Meey{e}=Mey;
    Meex{e}=Mex;
    Feyy{e}=T{e}'*Mey*T{e}*direcargay;
    Fexx{e}=T{e}'*Mex*T{e}*direcargax;
    Fe{e}=Fexx{e}+ Feyy{e};
    Kea{e} = T{e}'*Kloc*T{e};    
    K(GLe(e,:),GLe(e,:)) = K(GLe(e,:),GLe(e,:)) + Kea{e}; % ensambla Ke{e} en K global
    f(GLe(e,:))          = f(GLe(e,:))          + Fe{e}; % ensambla fe{e} en f global
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
    %if iteraciones==ite
    %    figure(1)
    %    xr=xyc(tipo_apoyo(ap,2),1);
    %    yr=xyc(tipo_apoyo(ap,2),2);
    %    [x,y,xn,yn]=apoyos(xr,yr,fesc,tipo_apoyo(ap,1));
    %    plot(x,y,'-r',xn,yn,'-r'), hold on;
    %    fill(x,y,'r'),hold on
    %    alpha(.5)
    %end
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
ac = aa; % desplazamientos conocidos en contorno
%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);        % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
% armo los vectores de desplazamientos (a) y fuerzas (q)
a  = zeros(ngdl,1);   a(c) = ac;  a(d) = ad; % desplazamientos
q  = zeros(ngdl,1);   q(c) = qd; % fuerzas nodales de equilibrio

u=a(1:3:nudos*3);
v=a(2:3:nudos*3);
x=[xye(:,1);xye(end,3)];
M1=zeros(elementos(1,1),1);
M2=zeros(elementos(1,1),1);
A1=zeros(elementos(1,1),1);
A2=zeros(elementos(1,1),1);
V1=zeros(elementos(1,1),1);
V2=zeros(elementos(1,1),1);
for e=1:elementos(1,1)   % para cada barra
    qe_glob{e} = T{e}*(Kea{e}*a(GLe(e,:)) - Fe{e});
    qe_loc{e}=T{e}*qe_glob{e};
    A1(e,1)=qe_loc{e}(1);
    V1(e,1)=qe_loc{e}(2);
    M1(e,1)=qe_loc{e}(3);
    A2(e,1)=qe_loc{e}(4);
    V2(e,1)=qe_loc{e}(5);
    M2(e,1)=qe_loc{e}(6);
    %x1 = xye(e,1);  x2 = xye(e,3);
    %y1 = xye(e,2);  y2 = xye(e,4);
    %deformada(AEe(e),EIe(e),Ace(e),Ce(e,:),T{e}*a(GLe(e,:)),Le(e),ang(e),kw(e),P(e),lalb(e,:),esc_def,esc_faxial,esc_V,esc_M,x1,y1,x2,y2,tipo_conti(e,:),iteraciones,ite,puntos_graficas)
    %P(e)=qe_glob{e}(1);
end
M=[M1;-M2(end)];
V=[V1;-V2(end)];
P=[A1;-A2(end)];

%figure
%plot(x,v)
%figure
%plot(x,M)
%figure
%plot(x,V)
%figure
%plot(x,P)

%figure(2); hold on; title('Deformada exagerada');    xlabel('x, m'); ylabel('y, m'); axis equal
%figure(3); hold on; title('Fuerza axial [kN]');      xlabel('x, m'); ylabel('y, m'); axis equal
%figure(4); hold on; title('Fuerza cortante [kN]');   xlabel('x, m'); ylabel('y, m'); axis equal
%figure(5); hold on; title('Momento flector [kN-m]'); xlabel('x, m'); ylabel('y, m'); axis equal


