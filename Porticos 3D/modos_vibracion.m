%% LUIS ENRIQUE GARCIA Figura 14-15 - Ejemplo 14-4 PAG 267
%k=10^3*[31.6850000000000,0,0,-36.1350000000000,0,0,4.71800000000000,0,0;0,29.7670000000000,0,0,-35.8580000000000,107.570000000000,0,6.59800000000000,-19.7940000000000;0,0,656.040000000000,0,0,-765.380000000000,0,0,117.180000000000;-36.1350000000000,0,0,90.2790000000000,0,0,-60.5530000000000,0,0;0,-35.8580000000000,0,0,88.2390000000000,-131.030000000000,0,-60.2860000000000,25.4180000000000;0,107.570000000000,-765.380000000000,0,-131.030000000000,2961.60000000000,0,25.4180000000000,-2137.80000000000;4.71800000000000,0,0,-60.5530000000000,0,0,96.7250000000000,0,0;0,6.59800000000000,0,0,-60.2860000000000,25.4180000000000,0,94.8450000000000,-6.10700000000000;0,-19.7940000000000,117.180000000000,0,25.4180000000000,-2137.80000000000,0,-6.10700000000000,3449];
%m=[29.4000000000000,0,0,0,0,0,0,0,0;0,29.4000000000000,0,0,0,0,0,0,0;0,0,208.250000000000,0,0,0,0,0,0;0,0,0,58.8000000000000,0,0,0,0,0;0,0,0,0,58.8000000000000,0,0,0,0;0,0,0,0,0,945.700000000000,0,0,0;0,0,0,0,0,0,58.8000000000000,0,0;0,0,0,0,0,0,0,58.8000000000000,0;0,0,0,0,0,0,0,0,945.700000000000];
%gama=[1,0,0;0,1,0;0,0,1;1,0,0;0,1,0;0,0,1;1,0,0;0,1,0;0,0,1];
%% defino parametros de espectro de respuesta NSR-10
%% parematrios para espectro medellin
Aa=0.15;%% aceleración pico efectiva
Av=0.2;%% velocidad pico efectiva
I=1;%% importancia
suelo=5;% suelo tipo D

%% matriz de rigidez entrada
k=kmodal; %%% kN/m
%% matriz de masa entrada
m=mmodal/9.8065174;
sizem=size(m,1);
%% grados de livertad por nudo
GL=6;
%% numero de modos de vibración 
n=50;
OPT.maxit = 10000;
OPT.tol   = 1e-32;
%OPT.issym = 1;  % La matriz Kdd es simetrica
[Phi,Omega2,flag] = eigs(k,m,n,'SM',OPT);% autovalores y autovectores
iteracionesss=80;
[Phi1,Omega1,T1]=estodola(k,m,iteracionesss,n);
if flag ~= 0
   warning('El calculo de valores y vectores propios no convergio');
end
Omega = sqrt(diag(Omega2));  % frecuencias angulares
[Omega,In] = sort(Omega);     % ordena las frecuencias
Te = 2*pi./Omega;             % seg - periodo de vibracion 
Phi = Phi(:,In); % ordena los modos segun el orden de las frecuencias
Phi = Phi./repmat(sqrt(diag(Phi'*m*Phi)'),size(Phi,1),1); % normaliza los modos
%% función que caclula la espectro de respuesta
%iteraciones=200;
%[Phi,Omega1,T1]=estodola(k,m,iteraciones,n);
sa=espectro(Aa,Av,I,suelo,Te);
def=sa*9.8067./(Omega).^2;
phim=(Phi'*m)';
meff=(Phi'*m*gama).^2;% claculo de masa modal efectiva 
alfa=(Phi'*m*gama);%coeficiente de participación 

sm=meff./sum(meff);%% participacín modal

meff1=(Phi1'*m*gama).^2;% claculo de masa modal efectiva
sm1=meff1./sum(meff1);%% participacín modal

alfa=((Phi'*m*gama).*def);%coeficiente de participación 


eta1=Phi*diag(alfa(:,1));

eta2=Phi*diag(alfa(:,2));
eta3=Phi*diag(alfa(:,3));

sizexy=size(xyc,1);
fuer=zeros(sizexy*GL,1);
%sizem=size(m,1);

defx=zeros(sizem,1);
defz=zeros(sizem,1);
defy=zeros(sizem,1);
for i=1:sizem
    defx(i,1)=rssq(eta1(i,:));
    defz(i,1)=rssq(eta2(i,:));
    defy(i,1)=rssq(eta3(i,:));
    %fuer(i,1)=rssq(d(i,:));
end
def1=zeros(sizexy(1,1)*GL,n);
def1(GLKM,1:n)=eta3;
qe_glob = cell(elementos(1,1),1);
qe_loc  = cell(elementos(1,1),1);
qe_loc_env  = cell(elementos(1,1),1);
%Ce = cell(elementos,1); 

for e=1:elementos(1,1)   % para cada barra
    qe_glob{e} = T{e}*(Ke{e}*def1(GLe(e,:),:));
    qe_loc{e}=T{e}*qe_glob{e};
    x1 = xye(e,1);  x2 = xye(e,4);
    y1 = xye(e,2);  y2 = xye(e,5);
    z1 = xye(e,3);  z2 = xye(e,6);
    qe_loc_env{e}=rssq(T{e}*qe_glob{e},2);
    %deformada(-we(e,:),T{e}*def1(GLe(e,:)),Le(e),esc_def,esc_faxial,esc_V,esc_M3,esc_M2,x1,y1,x2,y2,z1,z2,qe_loc{e},matT,EA(e),GJ(e),Acy(e),Acz(e),EIz(e),EIy(e));
end


%figure(2); hold on; title('Deformada exagerada');    xlabel('x, m'); ylabel('y, m'); axis equal
%figure(3); hold on; title('Fuerza axial [kN]');      xlabel('x, m'); ylabel('y, m'); axis equal
%figure(4); hold on; title('Fuerza cortante [kN]');   xlabel('x, m'); ylabel('y, m'); axis equal
%figure(5); hold on; title('Momento flector [kN-m]'); xlabel('x, m'); ylabel('y, m'); axis equal

