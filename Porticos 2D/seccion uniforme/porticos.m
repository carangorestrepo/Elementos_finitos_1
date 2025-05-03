function [qe_glob,ad,Kdd,Mdd,GLe,d,elementos,nudos,xye,MMexdd,T,Ke,ang,AEe,EIe,Ace,Le,P,iteraciones,ite,puntos_graficas,gama,Ce]=porticos(EG,mat_elem,xyc,ninf,sec,sec_ele,elemento_carga,nudo_carga,tipo_apoyo,tipo_conti,nudo_resorte,kw,lalb,PD,esc_def,esc_faxial,esc_V,esc_M,fesct,fesc)

% se definen algunas constantes que hacen el codigo mas legible
X = 1;
Y = 2;
G = 3;
XY = 12;
XYG=123;
g=9.8066502;      % aceleracion de la gravedad
%% porticos_calculos(EG,mat_elem,xyc,ninf,sec,sec_ele,elemento_carga,nudo_carga,tipo_apoyo,tipo_conti,X,Y)
elementos = size(ninf,1);
%% numero nudos
nudos = size(xyc,1);
%Grafico nudos
figure(1)
plot(xyc(:,1),xyc(:,2),'o','MarkerSize',6,...
'MarkerEdgeColor','blue',...
'MarkerFaceColor','blue'), hold on
xlabel('x(m)')
ylabel('y(m)')
h = text(xyc(:,1)+0.1,xyc(:,2)+0.1, num2str((1:nudos)'),'FontSize',fesct);
set(h, 'Color','blue');
%% numero grados de livertad
ngdl=nudos*3;
%% separo memoria
xye = zeros(elementos,4);
GLe = zeros(elementos,6);
Le = zeros(elementos,1);
AEe = zeros(elementos,1);
EIe = zeros(elementos,1);
Ace = zeros(elementos,1);
T = cell(elementos,1);
wex = zeros(elementos,2);
wey = zeros(elementos,2);
Wex = cell(elementos,1); 
Wey = cell(elementos,1); 
Cex = cell(elementos,1); 
Cey = cell(elementos,1); 
Ce = cell(elementos,1); 
ang = zeros(elementos,1); 
Fex = cell(elementos,1);
Fey = cell(elementos,1);
Fe = cell(elementos,1);
Ke = cell(elementos,1);
Kloce = cell(elementos,1);
Meex = cell(elementos,1);
Meey = cell(elementos,1);
Myy = cell(elementos,1);
Myyex = cell(elementos,1);
K = zeros(nudos*3);
M = zeros(nudos*3);
MMex = zeros(nudos*3);
My = zeros(nudos*3);
napoyo=size(tipo_apoyo,1);
%apoyo=zeros(napoyo,1);
ca = zeros(nudos*3,1);
f = zeros(nudos*3,1);
fp = zeros(nudos*3,1);%% cargas  y mometos puntuales
P = zeros(elementos,1);
qe_glob = cell(elementos(1,1),1);
qe_loc  = cell(elementos(1,1),1);
puntos_graficas=101; 
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
        fp(nudo_carga(i,2)*3-2,1)=nudo_carga(i,1);
    elseif  nudo_carga(i,3)==Y 
        fp(nudo_carga(i,2)*3-1,1)=nudo_carga(i,1);
        %% matriz de masa asociada a las cargas puntuales
        My(nudo_carga(i,2)*3-2,nudo_carga(i,2)*3-2)=abs(nudo_carga(i,1))/g;
        My(nudo_carga(i,2)*3-1,nudo_carga(i,2)*3-1)=abs(nudo_carga(i,1))/g;
        My(nudo_carga(i,2)*3,nudo_carga(i,2)*3)=abs(nudo_carga(i,1))/g;
    elseif nudo_carga(i,3)==G 
        fp(nudo_carga(i,2)*3,1)=nudo_carga(i,1);
    end
    M=My;
figure(1)
flecha(nudo_carga,xyc(nudo_carga(:,2),:))
end

%% Asignacion de resortes
sizenudo_resorte=size(nudo_resorte,1);
for i=1:sizenudo_resorte
    if nudo_resorte(i,3)==X
        K(nudo_resorte(i,2)*3-2,nudo_resorte(i,2)*3-2)=nudo_resorte(i,1);
    elseif nudo_resorte(i,3)==Y
        K(nudo_resorte(i,2)*3-1,nudo_resorte(i,2)*3-1)=nudo_resorte(i,1);
    elseif nudo_resorte(i,3)==G
        K(nudo_resorte(i,2)*3,nudo_resorte(i,2)*3)=nudo_resorte(i,1);
    end
end
%%
if strcmp(PD,'si')==1
    iteraciones=2;
else
    iteraciones=1;
end
for ite=1:iteraciones
    for e=1:elementos % para cada barra
        % saco los 4 cordenadas de la barra e
        xye(e,:) = [xyc(ninf(e,1),1),xyc(ninf(e,1),2),xyc(ninf(e,2),1),xyc(ninf(e,2),2)];
        GLe(e,:) = [(ninf(e,X)*3-2):(ninf(e,X)*3),(ninf(e,Y)*3-2):(ninf(e,Y)*3)];
        L = hypot(xye(e,3)-xye(e,1),xye(e,4)-xye(e,2));
        Le(e)=L;
        ang(e) = atan2(xye(e,4)-xye(e,2), xye(e,3)-xye(e,1));
        %% Se dibuja la estructura junto con su numeracion
        if iteraciones==ite
            figure(1)
            plot(xye(e,[1,3]),xye(e,[2,4]),'--k','LineWidth', 2),hold on
            %Calculo la posicion del centro de gravedad de la barra
            cgx = mean(xye(e,[1,3]));
            cgy = mean(xye(e,[2,4]));  
            figure(1)
            h = text(cgx+0.1, cgy+0.1, num2str(e),'FontSize',fesct);
            set(h, 'Color', [1 0 0]);
            %% grafico cargas distribuidas
            figure(1)
            cargas(wey(e,1),wey(e,2),xye(e,1),xye(e,3),xye(e,2),xye(e,4),0.1*wey(e,1),0.1*wey(e,2),Y)
            cargas(wex(e,1),wex(e,2),xye(e,1),xye(e,3),xye(e,2),xye(e,4),0.1*wex(e,1),0.1*wex(e,2),X)
        end
        %% materiales y secciones elementos
        AEe(e) = sec(sec_ele(e),1)*EG(mat_elem(e),1); %% AE
        EIe(e) = sec(sec_ele(e),2)*EG(mat_elem(e),1);%% EI
        Ace(e) = sec(sec_ele(e),3)*EG(mat_elem(e),2);%% Ac*G
        EA=AEe(e);
        EI=EIe(e);
        Ac=Ace(e);
        wx = wex(e,:);
        wy = wey(e,:);
        %% seco y coseno direcctor
        s = (xye(e,4)-xye(e,2))/L;   
        c = (xye(e,3)-xye(e,1))/L;
        % matriz de transformacion de coordenadas para la barra e   
        T{e} = [ c  s  0  0  0  0
               -s  c  0  0  0  0
                0  0  1  0  0  0
                0  0  0  c  s  0
                0  0  0 -s  c  0
                0  0  0  0  0  1];
        % matriz de rigidez nudos rigidos local expresada en el sistema de coordenadas locales
        % para la barra e
        [Kloc,~,Mex]=matriz_rigidez(tipo_conti(e,:),EA,Ac,EI,L,wex(e,1),wex(e,2),wex(e,1),wex(e,2),0,0,0,0,0,0,0,0,kw(e),P(e),lalb(e,:),puntos_graficas);
        [~,FEy,Mey]   =matriz_rigidez(tipo_conti(e,:),EA,Ac,EI,L,wey(e,1),wey(e,2),wey(e,1),wey(e,2),0,0,0,0,0,0,0,0,kw(e),P(e),lalb(e,:),puntos_graficas);
        Kloce{e}=Kloc;
        %[Kloc,FE,Me,V,M,fax,v,u]=fuerzas_empotramiento_bwinckler(EA,Ac,EI,L,we(e,1),we(e,2),0,0,0,0,0,0,0)
        %% direccion de carga distribuida
        %if we(e,4)==X   %DIRECCÓN HORIZONTAL
        direcargax=[-1;0;0;-1;0;0];
        Wex{e}=T{e}*[wex(e,1);0;0;wex(e,2);0;0];
        Cex{e}=[wx(1)*c;-wx(1)*s;-wx(1)*s;wx(2)*c;-wx(2)*s;-wx(2)*s];
        %if we(e,4)==Y     %DIRECCÓN VERTICAL
        direcargay=[0;-1;0;0;-1;0]; 
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
    %% calculo de  la matriz de masa concentrada

       Mlocy=  [abs(FEy(2))/g,              0,              0,              0,              0,            0;
                            0,  abs(FEy(2))/g,              0,              0,              0,            0;
                            0,              0,              0,              0,              0,            0;
                            0,              0,              0,  abs(FEy(5))/g,              0,            0;
                            0,              0,              0,              0,  abs(FEy(5))/g,            0;
                            0,              0,              0,              0,              0,           0];
        wa=wey(e,1);
        wb=wey(e,2);
         %% calculo de  la matriz de masa incluye grados de livertad rotacionales
        Myex=[[ (L*(3*wa + wb))/(12*g),                                                                                                                                     0,                                                                                                                                     0,   (L*(wa + wb))/(12*g),                                                                                                                                      0,                                                                                                                                      0];
              [                      0,      (L*(1260*EI^2*wa + 420*EI^2*wb + 10*Ac^2*L^4*wa + 3*Ac^2*L^4*wb + 224*Ac*EI*L^2*wa + 70*Ac*EI*L^2*wb))/(35*g*(Ac*L^2 + 12*EI)^2), (L^2*(1512*EI^2*wa + 1008*EI^2*wb + 15*Ac^2*L^4*wa + 7*Ac^2*L^4*wb + 300*Ac*EI*L^2*wa + 162*Ac*EI*L^2*wb))/(420*g*(Ac*L^2 + 12*EI)^2),                      0,                                                      (3*L*(wa + wb)*(3*Ac^2*L^4 + 84*Ac*EI*L^2 + 560*EI^2))/(140*g*(Ac*L^2 + 12*EI)^2),  -(L^2*(1512*EI^2*wa + 1008*EI^2*wb + 7*Ac^2*L^4*wa + 6*Ac^2*L^4*wb + 216*Ac*EI*L^2*wa + 162*Ac*EI*L^2*wb))/(420*g*(Ac*L^2 + 12*EI)^2)];
              [                      0, (L^2*(1512*EI^2*wa + 1008*EI^2*wb + 15*Ac^2*L^4*wa + 7*Ac^2*L^4*wb + 300*Ac*EI*L^2*wa + 162*Ac*EI*L^2*wb))/(420*g*(Ac*L^2 + 12*EI)^2),      (L^3*(504*EI^2*wa + 504*EI^2*wb + 5*Ac^2*L^4*wa + 3*Ac^2*L^4*wb + 96*Ac*EI*L^2*wa + 72*Ac*EI*L^2*wb))/(840*g*(Ac*L^2 + 12*EI)^2),                      0,   (L^2*(1008*EI^2*wa + 1512*EI^2*wb + 6*Ac^2*L^4*wa + 7*Ac^2*L^4*wb + 162*Ac*EI*L^2*wa + 216*Ac*EI*L^2*wb))/(420*g*(Ac*L^2 + 12*EI)^2),                                                       -(L^3*(wa + wb)*(Ac^2*L^4 + 28*Ac*EI*L^2 + 168*EI^2))/(280*g*(Ac*L^2 + 12*EI)^2)];
              [   (L*(wa + wb))/(12*g),                                                                                                                                     0,                                                                                                                                     0, (L*(wa + 3*wb))/(12*g),                                                                                                                                      0,                                                                                                                                      0];
              [                      0,                                                     (3*L*(wa + wb)*(3*Ac^2*L^4 + 84*Ac*EI*L^2 + 560*EI^2))/(140*g*(Ac*L^2 + 12*EI)^2),  (L^2*(1008*EI^2*wa + 1512*EI^2*wb + 6*Ac^2*L^4*wa + 7*Ac^2*L^4*wb + 162*Ac*EI*L^2*wa + 216*Ac*EI*L^2*wb))/(420*g*(Ac*L^2 + 12*EI)^2),                      0,       (L*(420*EI^2*wa + 1260*EI^2*wb + 3*Ac^2*L^4*wa + 10*Ac^2*L^4*wb + 70*Ac*EI*L^2*wa + 224*Ac*EI*L^2*wb))/(35*g*(Ac*L^2 + 12*EI)^2), -(L^2*(1008*EI^2*wa + 1512*EI^2*wb + 7*Ac^2*L^4*wa + 15*Ac^2*L^4*wb + 162*Ac*EI*L^2*wa + 300*Ac*EI*L^2*wb))/(420*g*(Ac*L^2 + 12*EI)^2)];
              [                      0, -(L^2*(1512*EI^2*wa + 1008*EI^2*wb + 7*Ac^2*L^4*wa + 6*Ac^2*L^4*wb + 216*Ac*EI*L^2*wa + 162*Ac*EI*L^2*wb))/(420*g*(Ac*L^2 + 12*EI)^2),                                                      -(L^3*(wa + wb)*(Ac^2*L^4 + 28*Ac*EI*L^2 + 168*EI^2))/(280*g*(Ac*L^2 + 12*EI)^2),                      0, -(L^2*(1008*EI^2*wa + 1512*EI^2*wb + 7*Ac^2*L^4*wa + 15*Ac^2*L^4*wb + 162*Ac*EI*L^2*wa + 300*Ac*EI*L^2*wb))/(420*g*(Ac*L^2 + 12*EI)^2),       (L^3*(504*EI^2*wa + 504*EI^2*wb + 3*Ac^2*L^4*wa + 5*Ac^2*L^4*wb + 72*Ac*EI*L^2*wa + 96*Ac*EI*L^2*wb))/(840*g*(Ac*L^2 + 12*EI)^2)]];
                 

        %% matriz de rigidez local en coordenadas globales
        Meex{e}=Mex;
        Meey{e}=Mey;
        
        %% matriz de maza trasformada a cordenadas globales
        Myy{e}=T{e}'*Mlocy*T{e};
        Myyex{e}=T{e}'*Myex*T{e};
        
        Fex{e}=T{e}'*Mex*T{e}*direcargax;
        Fey{e}=T{e}'*Mey*T{e}*direcargay;
        Fe{e}=Fex{e}+ Fey{e};
        Ke{e} = T{e}'*Kloc*T{e};    
        K(GLe(e,:),GLe(e,:)) = K(GLe(e,:),GLe(e,:)) + Ke{e}; % ensambla Ke{e} en K global
        if iteraciones==ite
             f(GLe(e,:))          = f(GLe(e,:))-  Fe{e}; % ensambla fe{e} en f global
            f(GLe(e,:))          = f(GLe(e,:))+  Fe{e}; % ensambla fe{e} en f global
        else
            f(GLe(e,:))          = fp(GLe(e,:)) + Fe{e}; % ensambla fe{e} en f global
            M(GLe(e,:),GLe(e,:)) = M(GLe(e,:),GLe(e,:)) + Myy{e}; % ensambla My{e} en M global
            MMex(GLe(e,:),GLe(e,:)) = MMex(GLe(e,:),GLe(e,:)) + Myyex{e}; % ensambla My{e} en M global
        end
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
        if iteraciones==ite
            figure(1)
            xr=xyc(tipo_apoyo(ap,2),1);
            yr=xyc(tipo_apoyo(ap,2),2);
            [x,y,xn,yn]=apoyos(xr,yr,fesc,tipo_apoyo(ap,1));
            plot(x,y,'-r',xn,yn,'-r'), hold on;
            fill(x,y,'r'),hold on
            alpha(.5)
        end
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
    MMexdd=MMex(d,d);
    
    val=zeros(nudos*3,3);
    recorrido1=1:3:nudos*3;
    valones=ones(nudos*3/3,1);
    val(recorrido1(1,:),1)=valones;
    val(recorrido1(1,:)+1,2)=valones;
    val(recorrido1(1,:)+2,3)=valones;
    %% matriz de coeficiones de partricpación pagina 290 luis enrique garcia
    gama=val(d,1:3);%% esta asocioado a los grados de livertad libres
    % armo los vectores de desplazamientos (a) y fuerzas (q)
    a  = zeros(ngdl,1);   a(c) = ac;  a(d) = ad; % desplazamientos
    q  = zeros(ngdl,1);   q(c) = qd; % fuerzas nodales de equilibrio
    for e=1:elementos(1,1)   % para cada barra
        qe_glob{e} = T{e}*(Ke{e}*a(GLe(e,:)) - Fe{e});
        qe_loc{e}=T{e}*qe_glob{e};
        x1 = xye(e,1);  x2 = xye(e,3);
        y1 = xye(e,2);  y2 = xye(e,4);
        deformada(AEe(e),EIe(e),Ace(e),Ce(e,:),T{e}*a(GLe(e,:)),Le(e),ang(e),kw(e),P(e),lalb(e,:),esc_def,esc_faxial,esc_V,esc_M,x1,y1,x2,y2,tipo_conti(e,:),iteraciones,ite,puntos_graficas)
        P(e)=-qe_glob{e}(1);
    end
    if strcmp(PD,'si')==1
        if ite==1
            K = zeros(nudos*3);
            Fe = cell(elementos,1);
            %f = zeros(nudos*3,1);
        end
    end
end
if ite==iteraciones
    figure(2); hold on; title('Deformada exagerada');    xlabel('x, m'); ylabel('y, m'); axis equal
    figure(3); hold on; title('Fuerza axial [kN]');      xlabel('x, m'); ylabel('y, m'); axis equal
    figure(4); hold on; title('Fuerza cortante [kN]');   xlabel('x, m'); ylabel('y, m'); axis equal
    figure(5); hold on; title('Momento flector [kN-m]'); xlabel('x, m'); ylabel('y, m'); axis equal
end

