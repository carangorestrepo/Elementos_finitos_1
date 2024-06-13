function [qe_loc,kmodal,mmodal,gama,GLKM,xye,elementos,nudos,GLe,we,Le,esc_def,esc_faxial,esc_V,esc_M3,esc_M2,matT,T,Ke]= porticos_3D(xyc,ninf,sec,tiposec,elecargas,carga,dcarga,carga_nudo,nudos_carga,nudosapoyos,tipoapoyo,EA,Acy,Acz,EIz,EIy,GJ)
%%escalas de Dibujo la estructura y su deformada
esc_def    = 500;          % escalamiento de la deformada
esc_faxial = 0.01;          % escalamiento del diagrama de axiales
esc_V      = 0.009;           % escalamiento del diagrama de cortantes
esc_M3      = 0.009;           % escalamiento del diagrama de momentos
esc_M2      = 0.01;           % escalamiento del diagrama de momentos
%cantidad elementos
elementos = size(ninf);
%% cantidad de nudos
nudos=size(unique(ninf),1);
hold on
plot3(xyc(:,1),xyc(:,2),xyc(:,3),'.k','MarkerSize',15)
for i=1:nudos
    h=text(xyc(i,1),xyc(i,2),xyc(i,3), num2str((i)'));%,'FontSize',2
    set(h, 'Color', [1 0 1]);
end
size_nudos_carga=size(nudos_carga,1);
%% cantidad de grados de livertad
%CGL = nudos(1,1)*3;
%% separo memoria
xye = zeros(elementos(1,1),6);
%L = zeros(elementos(1,1),1);
GLe = zeros(elementos(1,1),12);
GLep = zeros(size_nudos_carga(1,1),3);
bh = zeros(elementos(1,1),2);
%I = zeros(elementos(1,1),1);
%A = zeros(elementos(1,1),1);
%c = zeros(elementos(1,1),1);
%s = zeros(elementos(1,1),1);
%Ae= zeros(elementos(1,1),1);
%Ixze= zeros(elementos(1,1),1);
%Iyze= zeros(elementos(1,1),1);
Le = zeros(elementos(1,1),1);
ang = zeros(elementos(1,1),1);
T   = cell(elementos(1,1),1);
Ke  = cell(elementos(1,1),1);
Keloc  = cell(elementos(1,1),1);
MMe  = cell(elementos(1,1),1);
MMp  = sparse(nudos*6,nudos*6);
K   = sparse(nudos*6,nudos*6); 
M   = sparse(nudos*6,nudos*6); 
f = zeros(nudos*6,1); 
d = ones(nudos*6,1); 
Fe   = cell(elementos(1,1),1); 
we = zeros(elementos(1,1),2);
qe_glob = cell(elementos(1,1),1);
qe_loc  = cell(elementos(1,1),1);
wed = zeros(elementos(1,1),1);
fp = zeros(nudos*6,1); %%cargas puntuales
%% vector de cargas distribuidas
we(elecargas',1:2) = carga;
%% vector de dirección de cargas
wed(elecargas',1) = dcarga;
%% cargas puntuales 
fp(nudos_carga*6-4,1)=carga_nudo;
%% Matriz de masa cargas puntuales
for i=1:size_nudos_carga
        Mlocp=abs(carga_nudo(i,1))*[1,0,0;
                                    0,1,0;
                                    0,0,1];
       GLep(i,:)=((nudos_carga(i)*6-5):(nudos_carga(i)*6-3));           
       MMp(GLep(i,:),GLep(i,:))=MMp(GLep(i,:),GLep(i,:))+Mlocp;         
end

for e=1:elementos(1,1)
    % saco los 6 cordenadas de la barra e
    xye(e,:) = [xyc(ninf(e,1),1),xyc(ninf(e,1),2),xyc(ninf(e,1),3),xyc(ninf(e,2),1),xyc(ninf(e,2),2),xyc(ninf(e,2),3)];
    % saco los 6 gdls de la barra e
    GLe(e,:) = [(ninf(e,1)*6-5):(ninf(e,1)*6),(ninf(e,2)*6-5):(ninf(e,2)*6)];
    %% Se dibuja la estructura junto con su numeracion
    figure(1)
    plot3(xye(e,[1,4]),xye(e,[2,5]),xye(e,[3,6]),'--k','LineWidth', 2),hold on
    %Calculo la posicion del centro de gravedad de la barra
    cgx = mean(xye(e,[1,4]));
    cgy = mean(xye(e,[2,5]));  
    cgz = mean(xye(e,[3,6]));
    figure(1)
    h = text(cgx,cgy,cgz, num2str(e)); set(h, 'Color', [1 0 0]);
    %L = sqrt((xye(e,4)-xye(e,1))^2+(xye(e,5)-xye(e,2))^2+(xye(e,6)-xye(e,3))^2);
    ang(e) = atan2(xye(e,4)-xye(e,2), xye(e,3)-xye(e,1));
    Le(e) = sqrt((xye(e,4)-xye(e,1))^2+(xye(e,5)-xye(e,2))^2+(xye(e,6)-xye(e,3))^2);
    bh(e,:) = [sec(tiposec(e),1),sec(tiposec(e),2)];% inercias_y  I(m^4)
    %I = sec(tiposec(e),1)*sec(tiposec(e),2)^3/12;% inercias_y  I(m^4)
    %Ixz=bh(e,1).*bh(e,2).^3/12;
    %Iyz=bh(e,2).*bh(e,1).^3/12;   
    %Ixze(e) = Ixz;
    %Iyze(e) = Iyz;
    %A = sec(tiposec(e),1)*sec(tiposec(e),2);% area A(m^2)
    %EIz=Ixz*E;
    %EIy=Iyz*E;
    %EA=E*A;
    %GJ=G*Jc;
    tipo=1;
    %Ae(e) = A;
    %w = we(e);  
    cx=(xye(e,4)-xye(e,1))/Le(e);
    cz=(xye(e,5)-xye(e,2))/Le(e);%%%Y
    cy=(xye(e,6)-xye(e,3))/Le(e);%%%Z
    cxz=sqrt(cx.^2+cz.^2);
    %la = lalb(e,1);
    %lb = lalb(e,2);
    %bj=min([bh(e,1),bh(e,2)]);
    %hj=max([bh(e,1),bh(e,2)]);
    %J=(1/3-0.21*bj/hj*(1-1/12*(bj/hj)^4))*hj*bj^3;  
    % matriz de rigidez nudos rigidos local expresada en el sistema de coordenadas locales
    % para la barra e
    %Kloc1 =[[  (A*E)/L,                                                                                                                                                    0,                                                                                                                                                   0,        0,                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                   0, -(A*E)/L,                                                                                                                                                   0,                                                                                                                                                    0,        0,                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                   0];
    %       [        0,              -(12*A*E*G*Iyz*ks)/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Iyz)),                                                                                                                                                   0,        0,                                                                                                                                                                                                                                   0,                                                                                -(6*A*E*G*Iyz*ks*(L + la - lb))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Iyz)),        0,              (12*A*E*G*Iyz*ks)/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Iyz)),                                                                                                                                                    0,        0,                                                                                                                                                                                                                                   0,                                                                                -(6*A*E*G*Iyz*ks*(L - la + lb))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Iyz))];
    %       [        0,                                                                                                                                                    0,             -(12*A*E*G*Ixz*ks)/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Ixz)),        0,                                                                                 (6*A*E*G*Ixz*ks*(L + la - lb))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Ixz)),                                                                                                                                                                                                                                   0,        0,                                                                                                                                                   0,               (12*A*E*G*Ixz*ks)/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Ixz)),        0,                                                                                 (6*A*E*G*Ixz*ks*(L - la + lb))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Ixz)),                                                                                                                                                                                                                                   0];
    %       [        0,                                                                                                                                                    0,                                                                                                                                                   0,  (G*J)/L,                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                   0,        0,                                                                                                                                                   0,                                                                                                                                                    0, -(G*J)/L,                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                   0];
    %       [        0,                                                                                                                                                    0, (6*A*E*G*Ixz*ks*(L + la - lb))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Ixz)),        0,     -(4*E*Ixz*(A*G*ks*L^2 + A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 - A*G*ks*la*lb + A*G*ks*lb^2 + 3*E*Ixz))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Ixz)),                                                                                                                                                                                                                                   0,        0,                                                                                                                                                   0, -(6*A*E*G*Ixz*ks*(L + la - lb))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Ixz)),        0, -(2*E*Ixz*(A*G*ks*L^2 + A*G*ks*L*la + A*G*ks*L*lb - 2*A*G*ks*la^2 + 2*A*G*ks*la*lb - 2*A*G*ks*lb^2 - 6*E*Ixz))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Ixz)),                                                                                                                                                                                                                                   0];
    %       [        0, -(6*A*E*G*Iyz*ks*(L + la - lb))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Iyz)),                                                                                                                                                   0,        0,                                                                                                                                                                                                                                   0,     -(4*E*Iyz*(A*G*ks*L^2 + A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 - A*G*ks*la*lb + A*G*ks*lb^2 + 3*E*Iyz))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Iyz)),        0, (6*A*E*G*Iyz*ks*(L + la - lb))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Iyz)),                                                                                                                                                    0,        0,                                                                                                                                                                                                                                   0, -(2*E*Iyz*(A*G*ks*L^2 + A*G*ks*L*la + A*G*ks*L*lb - 2*A*G*ks*la^2 + 2*A*G*ks*la*lb - 2*A*G*ks*lb^2 - 6*E*Iyz))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Iyz))];
    %       [ -(A*E)/L,                                                                                                                                                    0,                                                                                                                                                   0,        0,                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                   0,  (A*E)/L,                                                                                                                                                   0,                                                                                                                                                    0,        0,                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                   0];
    %       [        0,               (12*A*E*G*Iyz*ks)/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Iyz)),                                                                                                                                                   0,        0,                                                                                                                                                                                                                                   0,                                                                                 (6*A*E*G*Iyz*ks*(L + la - lb))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Iyz)),        0,             -(12*A*E*G*Iyz*ks)/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Iyz)),                                                                                                                                                    0,        0,                                                                                                                                                                                                                                   0,                                                                                 (6*A*E*G*Iyz*ks*(L - la + lb))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Iyz))];
    %       [        0,                                                                                                                                                    0,              (12*A*E*G*Ixz*ks)/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Ixz)),        0,                                                                                -(6*A*E*G*Ixz*ks*(L + la - lb))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Ixz)),                                                                                                                                                                                                                                   0,        0,                                                                                                                                                   0,              -(12*A*E*G*Ixz*ks)/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Ixz)),        0,                                                                                -(6*A*E*G*Ixz*ks*(L - la + lb))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Ixz)),                                                                                                                                                                                                                                   0];
    %       [        0,                                                                                                                                                    0,                                                                                                                                                   0, -(G*J)/L,                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                   0,        0,                                                                                                                                                   0,                                                                                                                                                    0,  (G*J)/L,                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                   0];
    %       [        0,                                                                                                                                                    0, (6*A*E*G*Ixz*ks*(L - la + lb))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Ixz)),        0, -(2*E*Ixz*(A*G*ks*L^2 + A*G*ks*L*la + A*G*ks*L*lb - 2*A*G*ks*la^2 + 2*A*G*ks*la*lb - 2*A*G*ks*lb^2 - 6*E*Ixz))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Ixz)),                                                                                                                                                                                                                                   0,        0,                                                                                                                                                   0, -(6*A*E*G*Ixz*ks*(L - la + lb))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Ixz)),        0,     -(4*E*Ixz*(A*G*ks*L^2 - 2*A*G*ks*L*la + A*G*ks*L*lb + A*G*ks*la^2 - A*G*ks*la*lb + A*G*ks*lb^2 + 3*E*Ixz))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Ixz)),                                                                                                                                                                                                                                   0];
    %       [        0, -(6*A*E*G*Iyz*ks*(L - la + lb))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Iyz)),                                                                                                                                                   0,        0,                                                                                                                                                                                                                                   0, -(2*E*Iyz*(A*G*ks*L^2 + A*G*ks*L*la + A*G*ks*L*lb - 2*A*G*ks*la^2 + 2*A*G*ks*la*lb - 2*A*G*ks*lb^2 - 6*E*Iyz))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Iyz)),        0, (6*A*E*G*Iyz*ks*(L - la + lb))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Iyz)),                                                                                                                                                    0,        0,                                                                                                                                                                                                                                   0,     -(4*E*Iyz*(A*G*ks*L^2 - 2*A*G*ks*L*la + A*G*ks*L*lb + A*G*ks*la^2 - A*G*ks*la*lb + A*G*ks*lb^2 + 3*E*Iyz))/((la - L + lb)*(A*G*ks*L^2 - 2*A*G*ks*L*la - 2*A*G*ks*L*lb + A*G*ks*la^2 + 2*A*G*ks*la*lb + A*G*ks*lb^2 + 12*E*Iyz))]];

    Kloc=k_viga3d(Le(e),EA(e),GJ(e),Acy(e),Acz(e),EIz(e),EIy(e),tipo);
    Keloc{e}=Kloc;
    % matriz de transformacion de coordenadas para la barra e   
    if  cy~=1
           T{e} = [[           cx,  cy,           cz,            0,   0,            0,            0,   0,            0,            0,   0,            0];
                   [ -(cx*cy)/cxz, cxz, -(cy*cz)/cxz,            0,   0,            0,            0,   0,            0,            0,   0,            0];
                   [      -cz/cxz,   0,       cx/cxz,            0,   0,            0,            0,   0,            0,            0,   0,            0];
                   [            0,   0,            0,           cx,  cy,           cz,            0,   0,            0,            0,   0,            0];
                   [            0,   0,            0, -(cx*cy)/cxz, cxz, -(cy*cz)/cxz,            0,   0,            0,            0,   0,            0];
                   [            0,   0,            0,      -cz/cxz,   0,       cx/cxz,            0,   0,            0,            0,   0,            0];
                   [            0,   0,            0,            0,   0,            0,           cx,  cy,           cz,            0,   0,            0];
                   [            0,   0,            0,            0,   0,            0, -(cx*cy)/cxz, cxz, -(cy*cz)/cxz,            0,   0,            0];
                   [            0,   0,            0,            0,   0,            0,      -cz/cxz,   0,       cx/cxz,            0,   0,            0];
                   [            0,   0,            0,            0,   0,            0,            0,   0,            0,           cx,  cy,           cz];
                   [            0,   0,            0,            0,   0,            0,            0,   0,            0, -(cx*cy)/cxz, cxz, -(cy*cz)/cxz];
                   [            0,   0,            0,            0,   0,            0,            0,   0,            0,      -cz/cxz,   0,       cx/cxz]];
    else
    % matriz de transformacion de coordenadas para la barra e   
     T{e} = [[   0, cy, 0,   0,  0, 0,   0,  0, 0,   0,  0, 0];
             [ -cy,  0, 0,   0,  0, 0,   0,  0, 0,   0,  0, 0];
             [   0,  0, 1,   0,  0, 0,   0,  0, 0,   0,  0, 0];
             [   0,  0, 0,   0, cy, 0,   0,  0, 0,   0,  0, 0];
             [   0,  0, 0, -cy,  0, 0,   0,  0, 0,   0,  0, 0];
             [   0,  0, 0,   0,  0, 1,   0,  0, 0,   0,  0, 0];
             [   0,  0, 0,   0,  0, 0,   0, cy, 0,   0,  0, 0];
             [   0,  0, 0,   0,  0, 0, -cy,  0, 0,   0,  0, 0];
             [   0,  0, 0,   0,  0, 0,   0,  0, 1,   0,  0, 0];
             [   0,  0, 0,   0,  0, 0,   0,  0, 0,   0, cy, 0];
             [   0,  0, 0,   0,  0, 0,   0,  0, 0, -cy,  0, 0];
             [   0,  0, 0,   0,  0, 0,   0,  0, 0,   0,  0, 1]];
    end    
    % matriz de masa nudos local expresada en el sistema de coordenadas locales       
    %gama=w/A;
    Mloc=we(e,1)*Le(e)/2*[1,0,0,0,0,0,0,0,0,0,0,0;
                          0,1,0,0,0,0,0,0,0,0,0,0;
                          0,0,1,0,0,0,0,0,0,0,0,0;
                          0,0,0,0,0,0,0,0,0,0,0,0;
                          0,0,0,0,0,0,0,0,0,0,0,0;
                          0,0,0,0,0,0,0,0,0,0,0,0;
                          0,0,0,0,0,0,1,0,0,0,0,0;
                          0,0,0,0,0,0,0,1,0,0,0,0;
                          0,0,0,0,0,0,0,0,1,0,0,0;
                          0,0,0,0,0,0,0,0,0,0,0,0;
                          0,0,0,0,0,0,0,0,0,0,0,0;
                          0,0,0,0,0,0,0,0,0,0,0,0];
                 
    MMe{e} = T{e}'*Mloc*T{e};   
    M(GLe(e,:),GLe(e,:)) = M(GLe(e,:),GLe(e,:)) + MMe{e}; % ensambla Me{e} en M gl   
    %[Y1,Y2,M1,M2]=reaciones(L,w,1E+11,A,1E+11,1E+11,Ix,1E+11,la,lb,G,E,ks);   
    %% fuerzas de empotramiento distribuidas 
    %Ava=10^11;
    %Avb=A;
    %Avc=10^11;
    %Ixa=10^11;
    %Ixc=10^11;
    %Ixb=I;
    %Y1 =-(w*(12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L^3 - 12*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*la^3 + 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*la^3 + 12*Avb*Avc*E*Ixa*Ixb^2*Ixc^2*la^3 - 12*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*la^3 + 12*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*lb^3 - 12*Ava*Avb*E*Ixa^2*Ixb^2*Ixc*lb^3 - 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*lb^3 + 12*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*lb^3 + Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^5*ks + 12*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*L^2*la - 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L*la^2 - 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L^2*la + 12*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*L*la^2 - 36*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*L*lb^2 + 24*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*L^2*lb + 24*Ava*Avb*E*Ixa^2*Ixb^2*Ixc*L*lb^2 + 36*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L*lb^2 - 36*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L^2*lb - 24*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*L*lb^2 + 12*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*L^2*lb + Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^5 + Ava*Avb*Avc*G*Ixb^2*Ixc^2*ks*la^5 - Ava*Avb*Avc*G*Ixa^2*Ixb^2*ks*lb^5 - Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*lb^5 - 12*Ava*Avb*E*Ixa*Ixb^2*Ixc^2*la*lb^2 + 12*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*la*lb^2 + 12*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*la*lb^2 - 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*la*lb^2 + 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*la^2*lb - 12*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*la^2*lb - 12*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*la^2*lb + 12*Avb*Avc*E*Ixa^2*Ixb^2*Ixc*la^2*lb - 3*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la^4 - 3*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^4*ks*la + 2*Ava*Avb*Avc*G*Ixa^2*Ixb^2*L*ks*lb^4 + 5*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*lb^4 - 5*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^4*ks*lb - 3*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la*lb^4 + 3*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^4*lb + 2*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*la^3 + 2*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^3*ks*la^2 - 10*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*lb^3 + 10*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^3*ks*lb^2 - 2*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^2*lb^3 + 2*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^3*lb^2 - 2*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^5 + 2*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*lb^5 + 24*Ava*Avb*E*Ixa*Ixb^2*Ixc^2*L*la*lb - 24*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*L*la*lb - 24*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*L*la*lb + 24*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L*la*lb + 6*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la^2*lb^2 - 18*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*la*lb^2 - 6*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*la^2*lb + 3*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la^4 + 3*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^4*ks*la - 7*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*lb^4 + 5*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^4*ks*lb + 3*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la*lb^4 - 3*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^4*lb - 3*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la*lb^4 + 3*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la^4*lb + 3*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la*lb^4 - 3*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la^4*lb - 2*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^2*ks*la^3 - 2*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^3*ks*la^2 + 10*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^2*ks*lb^3 - 10*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^3*ks*lb^2 + 2*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^2*lb^3 - 2*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^3*lb^2 - 2*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la^2*lb^3 + 2*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la^3*lb^2 + 2*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la^2*lb^3 - 2*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la^3*lb^2 + 12*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la*lb^3 - 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la^3*lb + 12*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^3*ks*la*lb - 12*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la*lb^3 + 4*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la^3*lb - 12*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^3*ks*la*lb + 12*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L*ks*la*lb^3 - 4*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L*ks*la^3*lb + 12*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L^3*ks*la*lb - 12*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*la*lb^3 + 4*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*la^3*lb - 12*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^3*ks*la*lb - 6*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la^2*lb^2 + 18*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^2*ks*la*lb^2 + 6*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^2*ks*la^2*lb + 6*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L*ks*la^2*lb^2 - 18*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L^2*ks*la*lb^2 - 6*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L^2*ks*la^2*lb - 6*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*la^2*lb^2 + 18*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^2*ks*la*lb^2 + 6*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^2*ks*la^2*lb))/(2*(12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L^2 - 12*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*la^2 + 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*la^2 + 12*Avb*Avc*E*Ixa*Ixb^2*Ixc^2*la^2 - 12*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*la^2 - 12*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*lb^2 + 12*Ava*Avb*E*Ixa^2*Ixb^2*Ixc*lb^2 + 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*lb^2 - 12*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*lb^2 + 12*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*L*la - 24*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L*la + 12*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*L*la + 12*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*L*lb - 24*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L*lb + 12*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*L*lb + 12*Ava*Avb*E*Ixa*Ixb^2*Ixc^2*la*lb - 12*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*la*lb - 12*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*la*lb + 24*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*la*lb - 12*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*la*lb - 12*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*la*lb + 12*Avb*Avc*E*Ixa^2*Ixb^2*Ixc*la*lb + Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^4*ks + Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^4 + Ava*Avb*Avc*G*Ixb^2*Ixc^2*ks*la^4 + Ava*Avb*Avc*G*Ixa^2*Ixb^2*ks*lb^4 + Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*lb^4 - 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la^3 - 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^3*ks*la - 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*lb^3 - 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^3*ks*lb + 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la*lb^3 + 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^3*lb + 6*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*la^2 + 6*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*lb^2 + 6*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^2*lb^2 - 2*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^4 - 2*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*lb^4 + 4*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la^3 + 4*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^3*ks*la + 4*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*lb^3 + 4*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^3*ks*lb - 4*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la*lb^3 - 4*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^3*lb + 4*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la*lb^3 + 4*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la^3*lb - 4*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la*lb^3 - 4*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la^3*lb - 6*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^2*ks*la^2 - 6*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^2*ks*lb^2 - 6*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^2*lb^2 + 6*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la^2*lb^2 - 6*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la^2*lb^2 - 12*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la*lb^2 - 12*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la^2*lb + 12*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*la*lb + 12*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la*lb^2 + 12*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la^2*lb - 12*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^2*ks*la*lb - 12*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L*ks*la*lb^2 - 12*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L*ks*la^2*lb + 12*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L^2*ks*la*lb + 12*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*la*lb^2 + 12*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*la^2*lb - 12*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^2*ks*la*lb));
    %Y2 =(w*(12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L^3 + 12*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*la^3 - 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*la^3 - 12*Avb*Avc*E*Ixa*Ixb^2*Ixc^2*la^3 + 12*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*la^3 - 12*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*lb^3 + 12*Ava*Avb*E*Ixa^2*Ixb^2*Ixc*lb^3 + 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*lb^3 - 12*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*lb^3 + Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^5*ks - 24*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*L*la^2 + 12*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*L^2*la + 36*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L*la^2 - 36*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L^2*la + 24*Avb*Avc*E*Ixa*Ixb^2*Ixc^2*L*la^2 - 36*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*L*la^2 + 24*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*L^2*la + 12*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*L*lb^2 - 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L*lb^2 - 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L^2*lb + 12*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*L^2*lb - Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^5 - Ava*Avb*Avc*G*Ixb^2*Ixc^2*ks*la^5 + Ava*Avb*Avc*G*Ixa^2*Ixb^2*ks*lb^5 + Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*lb^5 + 12*Ava*Avb*E*Ixa*Ixb^2*Ixc^2*la*lb^2 - 12*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*la*lb^2 - 12*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*la*lb^2 + 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*la*lb^2 - 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*la^2*lb + 12*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*la^2*lb + 12*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*la^2*lb - 12*Avb*Avc*E*Ixa^2*Ixb^2*Ixc*la^2*lb + 5*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la^4 - 5*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^4*ks*la + 2*Ava*Avb*Avc*G*Ixb^2*Ixc^2*L*ks*la^4 - 3*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*lb^4 - 3*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^4*ks*lb + 3*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la*lb^4 - 3*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^4*lb - 10*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*la^3 + 10*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^3*ks*la^2 + 2*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*lb^3 + 2*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^3*ks*lb^2 + 2*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^2*lb^3 - 2*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^3*lb^2 + 2*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^5 - 2*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*lb^5 + 24*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L*la*lb - 24*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*L*la*lb - 24*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*L*la*lb + 24*Avb*Avc*E*Ixa^2*Ixb^2*Ixc*L*la*lb + 6*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la^2*lb^2 - 6*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*la*lb^2 - 18*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*la^2*lb - 7*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la^4 + 5*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^4*ks*la + 3*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*lb^4 + 3*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^4*ks*lb - 3*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la*lb^4 + 3*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^4*lb + 3*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la*lb^4 - 3*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la^4*lb - 3*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la*lb^4 + 3*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la^4*lb + 10*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^2*ks*la^3 - 10*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^3*ks*la^2 - 2*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^2*ks*lb^3 - 2*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^3*ks*lb^2 - 2*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^2*lb^3 + 2*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^3*lb^2 + 2*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la^2*lb^3 - 2*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la^3*lb^2 - 2*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la^2*lb^3 + 2*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la^3*lb^2 - 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la*lb^3 + 12*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la^3*lb + 12*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^3*ks*la*lb + 4*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la*lb^3 - 12*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la^3*lb - 12*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^3*ks*la*lb - 4*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L*ks*la*lb^3 + 12*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L*ks*la^3*lb + 12*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L^3*ks*la*lb + 4*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*la*lb^3 - 12*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*la^3*lb - 12*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^3*ks*la*lb - 6*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la^2*lb^2 + 6*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^2*ks*la*lb^2 + 18*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^2*ks*la^2*lb + 6*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L*ks*la^2*lb^2 - 6*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L^2*ks*la*lb^2 - 18*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L^2*ks*la^2*lb - 6*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*la^2*lb^2 + 6*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^2*ks*la*lb^2 + 18*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^2*ks*la^2*lb))/(2*(12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L^2 - 12*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*la^2 + 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*la^2 + 12*Avb*Avc*E*Ixa*Ixb^2*Ixc^2*la^2 - 12*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*la^2 - 12*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*lb^2 + 12*Ava*Avb*E*Ixa^2*Ixb^2*Ixc*lb^2 + 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*lb^2 - 12*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*lb^2 + 12*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*L*la - 24*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L*la + 12*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*L*la + 12*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*L*lb - 24*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L*lb + 12*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*L*lb + 12*Ava*Avb*E*Ixa*Ixb^2*Ixc^2*la*lb - 12*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*la*lb - 12*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*la*lb + 24*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*la*lb - 12*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*la*lb - 12*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*la*lb + 12*Avb*Avc*E*Ixa^2*Ixb^2*Ixc*la*lb + Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^4*ks + Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^4 + Ava*Avb*Avc*G*Ixb^2*Ixc^2*ks*la^4 + Ava*Avb*Avc*G*Ixa^2*Ixb^2*ks*lb^4 + Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*lb^4 - 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la^3 - 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^3*ks*la - 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*lb^3 - 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^3*ks*lb + 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la*lb^3 + 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^3*lb + 6*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*la^2 + 6*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*lb^2 + 6*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^2*lb^2 - 2*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^4 - 2*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*lb^4 + 4*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la^3 + 4*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^3*ks*la + 4*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*lb^3 + 4*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^3*ks*lb - 4*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la*lb^3 - 4*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^3*lb + 4*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la*lb^3 + 4*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la^3*lb - 4*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la*lb^3 - 4*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la^3*lb - 6*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^2*ks*la^2 - 6*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^2*ks*lb^2 - 6*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^2*lb^2 + 6*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la^2*lb^2 - 6*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la^2*lb^2 - 12*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la*lb^2 - 12*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la^2*lb + 12*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*la*lb + 12*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la*lb^2 + 12*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la^2*lb - 12*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^2*ks*la*lb - 12*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L*ks*la*lb^2 - 12*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L*ks*la^2*lb + 12*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L^2*ks*la*lb + 12*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*la*lb^2 + 12*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*la^2*lb - 12*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^2*ks*la*lb));
    %M1 =(w*(12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L^4 - 12*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*la^4 + 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*la^4 + 12*Avb*Avc*E*Ixa*Ixb^2*Ixc^2*la^4 - 12*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*la^4 - 12*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*lb^4 + 12*Ava*Avb*E*Ixa^2*Ixb^2*Ixc*lb^4 + 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*lb^4 - 12*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*lb^4 + 36*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*L^2*la^2 - 72*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L^2*la^2 + 36*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*L^2*la^2 - 108*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*L^2*lb^2 + 72*Ava*Avb*E*Ixa^2*Ixb^2*Ixc*L^2*lb^2 + 72*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L^2*lb^2 - 36*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*L^2*lb^2 - 36*Ava*Avb*E*Ixa*Ixb^2*Ixc^2*la^2*lb^2 + 36*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*la^2*lb^2 + 36*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*la^2*lb^2 - 72*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*la^2*lb^2 + 36*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*la^2*lb^2 + 36*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*la^2*lb^2 - 36*Avb*Avc*E*Ixa^2*Ixb^2*Ixc*la^2*lb^2 + Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^6*ks - 24*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*L*la^3 + 24*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L*la^3 + 24*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L^3*la - 24*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*L^3*la + 72*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*L*lb^3 + 48*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*L^3*lb - 72*Ava*Avb*E*Ixa^2*Ixb^2*Ixc*L*lb^3 - 48*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L*lb^3 - 48*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L^3*lb + 48*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*L*lb^3 + Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^6 + Ava*Avb*Avc*G*Ixb^2*Ixc^2*ks*la^6 + Ava*Avb*Avc*G*Ixa^2*Ixb^2*ks*lb^6 + Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*lb^6 - 24*Ava*Avb*E*Ixa*Ixb^2*Ixc^2*la^3*lb + 24*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*la^3*lb + 24*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*la^3*lb - 24*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*la*lb^3 - 24*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*la^3*lb + 24*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*la*lb^3 + 24*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*la*lb^3 - 24*Avb*Avc*E*Ixa^2*Ixb^2*Ixc*la*lb^3 - 6*Ava*Avb*Avc*G*Ixa^2*Ixb^2*L*ks*lb^5 - 6*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*lb^5 - 6*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^5*ks*lb + 72*Ava*Avb*E*Ixa*Ixb^2*Ixc^2*L*la^2*lb - 72*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*L*la^2*lb - 72*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*L*la^2*lb + 72*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L*la*lb^2 + 144*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L*la^2*lb - 72*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L^2*la*lb - 72*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*L*la*lb^2 - 72*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*L*la^2*lb + 72*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*L^2*la*lb - 72*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*L*la*lb^2 - 72*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*L*la^2*lb + 72*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*L^2*la*lb + 72*Avb*Avc*E*Ixa^2*Ixb^2*Ixc*L*la*lb^2 + 72*Avb*Avc*E*Ixa^2*Ixb^2*Ixc*L*la^2*lb - 72*Avb*Avc*E*Ixa^2*Ixb^2*Ixc*L^2*la*lb - 9*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*la^4 + 16*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^3*ks*la^3 - 9*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^4*ks*la^2 + 6*Ava*Avb*Avc*G*Ixa^2*Ixb^2*L^2*ks*lb^4 + 15*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*lb^4 - 20*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^3*ks*lb^3 + 15*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^4*ks*lb^2 - 9*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^2*lb^4 - 16*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^3*lb^3 - 9*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^4*lb^2 - 2*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^6 - 2*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*lb^6 + 36*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la^2*lb^3 + 48*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la^3*lb^2 - 48*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*la^3*lb + 36*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^3*ks*la^2*lb + 12*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*lb^5 + 6*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^5*ks*lb - 54*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*la^2*lb^2 + 9*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^2*ks*la^4 - 16*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^3*ks*la^3 + 9*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^4*ks*la^2 - 21*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^2*ks*lb^4 + 20*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^3*ks*lb^3 - 15*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^4*ks*lb^2 + 9*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^2*lb^4 + 16*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^3*lb^3 + 9*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^4*lb^2 - 9*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la^2*lb^4 - 16*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la^3*lb^3 - 9*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la^4*lb^2 + 9*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la^2*lb^4 + 16*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la^3*lb^3 + 9*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la^4*lb^2 + 18*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la^4*lb + 54*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^2*ks*la^2*lb^2 - 54*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L^2*ks*la^2*lb^2 + 54*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^2*ks*la^2*lb^2 - 18*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la^4*lb + 18*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L*ks*la^4*lb - 18*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*la^4*lb - 36*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la^2*lb^3 - 48*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la^3*lb^2 + 48*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^2*ks*la^3*lb - 36*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^3*ks*la^2*lb + 36*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L*ks*la^2*lb^3 + 48*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L*ks*la^3*lb^2 - 48*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L^2*ks*la^3*lb + 36*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L^3*ks*la^2*lb - 36*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*la^2*lb^3 - 48*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*la^3*lb^2 + 48*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^2*ks*la^3*lb - 36*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^3*ks*la^2*lb))/(12*(12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L^2 - 12*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*la^2 + 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*la^2 + 12*Avb*Avc*E*Ixa*Ixb^2*Ixc^2*la^2 - 12*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*la^2 - 12*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*lb^2 + 12*Ava*Avb*E*Ixa^2*Ixb^2*Ixc*lb^2 + 12*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*lb^2 - 12*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*lb^2 + 12*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*L*la - 24*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L*la + 12*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*L*la + 12*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*L*lb - 24*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*L*lb + 12*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*L*lb + 12*Ava*Avb*E*Ixa*Ixb^2*Ixc^2*la*lb - 12*Ava*Avb*E*Ixa^2*Ixb*Ixc^2*la*lb - 12*Ava*Avc*E*Ixa*Ixb^2*Ixc^2*la*lb + 24*Ava*Avc*E*Ixa^2*Ixb*Ixc^2*la*lb - 12*Ava*Avc*E*Ixa^2*Ixb^2*Ixc*la*lb - 12*Avb*Avc*E*Ixa^2*Ixb*Ixc^2*la*lb + 12*Avb*Avc*E*Ixa^2*Ixb^2*Ixc*la*lb + Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^4*ks + Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^4 + Ava*Avb*Avc*G*Ixb^2*Ixc^2*ks*la^4 + Ava*Avb*Avc*G*Ixa^2*Ixb^2*ks*lb^4 + Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*lb^4 - 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la^3 - 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^3*ks*la - 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*lb^3 - 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^3*ks*lb + 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la*lb^3 + 4*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^3*lb + 6*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*la^2 + 6*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*lb^2 + 6*Ava*Avb*Avc*G*Ixa^2*Ixc^2*ks*la^2*lb^2 - 2*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^4 - 2*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*lb^4 + 4*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la^3 + 4*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^3*ks*la + 4*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*lb^3 + 4*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^3*ks*lb - 4*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la*lb^3 - 4*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^3*lb + 4*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la*lb^3 + 4*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la^3*lb - 4*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la*lb^3 - 4*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la^3*lb - 6*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^2*ks*la^2 - 6*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^2*ks*lb^2 - 6*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*ks*la^2*lb^2 + 6*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*ks*la^2*lb^2 - 6*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*ks*la^2*lb^2 - 12*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la*lb^2 - 12*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L*ks*la^2*lb + 12*Ava*Avb*Avc*G*Ixa^2*Ixc^2*L^2*ks*la*lb + 12*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la*lb^2 + 12*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L*ks*la^2*lb - 12*Ava*Avb*Avc*G*Ixa*Ixb*Ixc^2*L^2*ks*la*lb - 12*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L*ks*la*lb^2 - 12*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L*ks*la^2*lb + 12*Ava*Avb*Avc*G*Ixa*Ixb^2*Ixc*L^2*ks*la*lb + 12*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*la*lb^2 + 12*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L*ks*la^2*lb - 12*Ava*Avb*Avc*G*Ixa^2*Ixb*Ixc*L^2*ks*la*lb));  
    [Pa,Pb,Ra,Rb,Ma,Mb]=M_V_empo(we(e),we(e),Acy(e),EIz(e),Le(e),EA(e),we(e,1),we(e,2)); 
    %cl=L-la-lb;
    %Pa=w*L/2;
    %Ra=w*cl/2+w*la;
    %Ma=w*cl^2/12+w*cl/2*la+w*la^2/2;
    %Pb=w*L/2;
    %Rb=w*cl/2+w*lb;
    %Mb=w*cl^2/12+w*cl/2*lb+w*lb^2/2;
    %% matriz de mometos de empotramiento
    Me(:,:)=[-Pa,0,0,0,0,0,0,0,0,0,0,0;
               0,-Ra,0,0,0,0,0,0,0,0,0,0;
               0,0,-Ra,0,0,0,0,0,0,0,0,0;
               0,0,0,-Pa,0,0,0,0,0,0,0,0;
               0,0,-Ma,0,0,0,0,0,0,0,0,0;
               0,-Ma,0,0,0,0,0,0,0,0,0,0;
               0,0,0,0,0,0,-Pb,0,0,0,0,0;
               0,0,0,0,0,0,0,-Rb,0,0,0,0;
               0,0,0,0,0,0,0,0,-Rb,0,0,0;
               0,0,0,0,0,0,0,0,0,-Pb,0,0;
               0,0,0,0,0,0,0,0,-Mb,0,0,0;
               0,0,0,0,0,0,0,Mb,0,0,0,0];
    %% direccion de carga
    if wed(e)==1
        direcarga=[0;1;0;0;0;0;0;1;0;0;0;0];
    else
        direcarga=[-1;0;0;0;0;0;-1;0;0;0;0;0];
    end  
    %% matriz de rigidez local en coordenadas globales
    Fe{e}=T{e}'*Me*T{e}*direcarga;
    Ke{e} = T{e}'*Kloc*T{e};    
    K(GLe(e,:),GLe(e,:)) = K(GLe(e,:),GLe(e,:)) + Ke{e}; % ensambla Ke{e} en K global
    f(GLe(e,:))          = f(GLe(e,:))          + Fe{e}; % ensambla fe{e} en f global
end
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
axis equal
grid minor
%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
sizetipoapoyo=size(tipoapoyo);
D = -diag(ones(1,nudos*6));
for i=1:sizetipoapoyo(1,2)
    if tipoapoyo(i)==1
        K(:,nudosapoyos(i)*6-5)=D(:,nudosapoyos(i)*6-5);%% apoyado horizontal
        figure(1)
        plot(x,y,'b'), hold on
        figure(1)
        plot(xyc(nudosapoyos(i),1),xyc(nudosapoyos(i),2),'>','MarkerSize',8,...
       'MarkerEdgeColor','blue'), hold on
        figure(1)
        plot(xyc(nudosapoyos(i),1),xyc(nudosapoyos(i),2),'o','MarkerSize',12,...
       'MarkerEdgeColor','blue'), hold on
        d1 = [0;1;1];
        d((nudosapoyos(i)*6-5):(nudosapoyos(i)*6),1)=d1;
    elseif tipoapoyo(i)==2 %% apoyado vertical
        K(:,nudosapoyos(i)*6-4)=D(:,nudosapoyos(i)*6-4);%% apoyado vertical
        x=[xyc(nudosapoyos(i),1)-0.1,0.1+xyc(nudosapoyos(i),1)];  
        y=[xyc(nudosapoyos(i),2),xyc(nudosapoyos(i),2)]; 
        figure(1)
        plot(xyc(nudosapoyos(i),1),xyc(nudosapoyos(i),2),'^','MarkerSize',8,...
       'MarkerEdgeColor','blue'), hold on
        figure(1)
        plot(xyc(nudosapoyos(i),1),xyc(nudosapoyos(i),2),'o','MarkerSize',12,...
       'MarkerEdgeColor','blue'), hold on
       d1=[1;0;1];
       d((nudosapoyos(i)*6-5):(nudosapoyos(i)*6),1)=d1;
    elseif tipoapoyo(i)==12 %% apoyado simple
       K(:,(nudosapoyos(i)*6-5):(nudosapoyos(i)*5-4))=D(:,(nudosapoyos(i)*6-5):(nudosapoyos(i)*6-4));%horizontal y vertical
       figure(1)
       plot(xyc(nudosapoyos(i),1),xyc(nudosapoyos(i),2),'^','MarkerSize',7,...
       'MarkerEdgeColor','blue',...
       'MarkerFaceColor','blue'), hold on
       d1=[0;0;1];
       d((nudosapoyos(i)*6-5):(nudosapoyos(i)*6),1)=d1;
     elseif tipoapoyo(i)==123  %% empotramiento
       K(:,(nudosapoyos(i)*6-5):(nudosapoyos(i)*6))=D(:,(nudosapoyos(i)*6-5):(nudosapoyos(i)*6));%horizontal y vertica y giro
       figure(1)
       plot3(xyc(nudosapoyos(i),1),xyc(nudosapoyos(i),2),xyc(nudosapoyos(i),3),'-s','MarkerSize',15,...
       'MarkerEdgeColor','blue',...
       'MarkerFaceColor','blue'), hold on
       d1=[0;0;0;0;0;0];
       d((nudosapoyos(i)*6-5):(nudosapoyos(i)*6),1)=d1;
    end 
end
%% armado matriz de rigidez y matriz de masa   para calcular modos de vibración
GLKM=find(d==1);
%% armado de vectores de desplazamientos para caclulo de modos de vibración 
val=zeros(nudos*6,3);
recorrido1=1:6:(nudos*6);
valones=ones(nudos*6/6,1);
val(recorrido1(1,:),1)=valones;
val(recorrido1(1,:)+1,2)=valones;
val(recorrido1(1,:)+2,3)=valones;
%val(recorrido1(1,:)+3,4)=valones;
%val(recorrido1(1,:)+4,5)=valones;
%val(recorrido1(1,:)+5,6)=valones;
kmodal=K(GLKM,GLKM);%% matriz de rigidez para calcular modos de vibración
M=M+MMp;
mmodal=M(GLKM,GLKM);%% matriz de masa para calcular modos de vibración
%% matriz de coeficiones de partricpación pagina 290 luis enrique garcia
gama=val(GLKM,:);
%% resuelvo el sistema de ecuaciones
RD=K\(f+fp);
Def=RD.*d;
%% grafico diagramas para portico
figure(1)
for e=1:elementos(1,1)  
    qe_glob{e} = T{e}*(Ke{e}*Def(GLe(e,:)) - Fe{e});
    qe_loc{e}=T{e}*qe_glob{e};
    x1 = xye(e,1);  x2 = xye(e,4);
    y1 = xye(e,2);  y2 = xye(e,5);
    z1 = xye(e,3);  z2 = xye(e,6);
    matT = T{e}(1:3,1:3)';
%    deformada(-we(e,:),T{e}*Def(GLe(e,:)),Le(e),esc_def,esc_faxial,esc_V,esc_M3,esc_M2,x1,y1,x2,y2,z1,z2,qe_loc{e},matT,EA(e),GJ(e),Acy(e),Acz(e),EIz(e),EIy(e));
end
a=1


