k=[216760,-306770,105490,-19561,4282.20000000000,-510.880000000000;-306770,668240,-475140,137940,-29375,5385.70000000000;105490,-475140,731370,-493230,159600,-29327;-19561,137940,-493230,749020,-494470,145710;4282.20000000000,-29375,159600,-494470,738110,-515900;-510.880000000000,5385.70000000000,-29375,145710,-515900,889940];
m=[256,0,0,0,0,0;0,256,0,0,0,0;0,0,256,0,0,0;0,0,0,256,0,0;0,0,0,0,256,0;0,0,0,0,0,256];


%% defino parametros de espectro de respuesta NSR-10
%% parematrios para espectro medellin
Aa=0.15;%% aceleración pico efectiva
Av=0.2;%% velocidad pico efectiva
I=1;%% importancia
suelo=5;% suelo tipo D

%% matriz de rigidez entrada
%k=kmodal; %%% kN/m
%% matriz de masa entrada
%m=mmodal/9.8065174;
sizem=size(m);
%% grados de livertad por nudo
GL=3;
%% numero de modos de vibración 
n=6;
OPT.maxit = 10000;
OPT.tol   = 1e-32;
%OPT.issym = 1;  % La matriz Kdd es simetrica
[Phi,Omega2,flag] = eigs(k,m,n,'SM',OPT);% autovalores y autovectores
if flag ~= 0
   warning('El calculo de valores y vectores propios no convergio');
end
Omega = sqrt(diag(Omega2));  % frecuencias angulares
[Omega,In] = sort(Omega);     % ordena las frecuencias
T = 2*pi./Omega;             % seg - periodo de vibracion 
Phi = Phi(:,In); % ordena los modos segun el orden de las frecuencias
Phi = Phi./repmat(sqrt(diag(Phi'*m*Phi)'),size(Phi,1),1); % normaliza los modos
%% función que caclula la espectro de respuesta
%iteraciones=200;
%[Phi,Omega1,T1]=estodola(k,m,iteraciones,n);
%sa=espectro(Aa,Av,I,suelo,T);
gama=ones(6,1);
phim=(Phi'*m)';
meff=(Phi'*m*gama).^2;% claculo de masa modal efectiva 
alfa=(Phi'*m*gama);%coeficiente de participación 
sm=meff./sum(meff);%% participacín modal
F=zeros(sizem(1,1),n);
%% (15-8)coeficuiente de participaciòn amplñificado pagina 286 luis enrique garcia
%r=diag(sa*9.806./(2*pi./T).^2);
%valalfa=gama*alfa'*r;
%Umod1=zeros(sizem(1,1),n);
%Fmod=zeros(sizem(1,1),n);
%for i=1:n
%    Umod1(:,i)=Phi(:,i).*valalfa(:,i);%%(15.10) desplazamietos  dinamicos modales 286 luis enrique garcia pagina 288
%    Fmod(:,i)=k*Phi(:,i).*valalfa(:,i);%%(15.11) fuerzas modales pagina 286 luis enrique garcia pagina 288
    %Umod=Phi*diag(alfa.*sa*9.806./(2*pi./T).^2);
    %F1=k*Phi*diag(alfa.*sa*9.806./(2*pi./T).^2);
%end
%Ugtt=zeros(nudos(1,1)*GL,n);
%Phia=zeros(nudos(1,1)*GL,n);
%kgttm=zeros(nudos(1,1)*GL,n);
%Ugtt(GLKM,:)=Umod1;% desplazamietos modales
%kgttm(GLKM,GLKM)=k;% desplazamietos modales
%Fb=(kgttm)*Ugtt; %fuerzas sismicas
%Phia(GLKM,:)=Phi;

n=[-3.5459;
-0.15329;
-0.001903;
-0.0035796;
-0.0012279;
-0.00035675];

U=(Phi*diag(n));
F=k*U

