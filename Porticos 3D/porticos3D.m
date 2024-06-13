%escamilla 11.3
% coordenadas de cada nodo [x, y]
xy=[0,3,3;5,3,3;0,3,0;0,0,3];
% nudos con restricciones
ap=[3,2,4];
w=[24,35,0];
tipap=ones(1,3)*3;
bh=[0.3,0.4;0.4,0.25;0.3,0.4];
nudosele=[1,2;4,1;3,1];


%xy=[0,0,0;0,0,4;6,0,0;6,0,4;0,5,0;0,5,4;6,5,0;6,5,4];
%ap=[1,3,5,7];
%w=[0,0,0,0,10,10,10,10];
%tipap=ones(1,4)*3;
%bh=ones(8,2)*0.4;
%nudosele=[1,2;3,4;5,6;7,8;2,4;4,8;8,6;6,2];


%xy=[0,0,0;0,0,3;0,0,6;6,0,0;6,0,3;6,0,6;12,0,0;12,0,3;12,0,6;0,4,0;0,4,3;0,4,6;6,4,0;6,4,3;6,4,6;12,4,0;12,4,3;12,4,6];
%ap=[1,4,7,10,13,16];
%w=[0;0;0;0;0;0;0;0;0;0;0;0;30;30;30;30;30;30;30;30;30;30;30;30;30;30]';
%tipap=ones(1,6)*3;
%bh=ones(26,2)*0.4;
%nudosele=[1,2;2,3;4,5;5,6;7,8;8,9;10,11;11,12;13,14;14,15;16,17;17,18;2,5;5,8;11,14;14,17;2,11;5,14;8,17;3,12;6,15;9,18;3,6;6,9;12,15;15,18];
%xy=[0,0,0;0,0,4;0,0,8;0,4,0;0,4,4;0,4,8;0,13,0;0,13,4;0,13,8;6,0,0;6,0,4;6,0,8;6,4,0;6,4,4;6,4,8;6,13,0;6,13,4;6,13,8];



%escamilla 11.3
%xy=[0,3,3;5,3,3;0,3,0;0,0,3];
%ap=[2,3,4];
%w=[24,35,0];
%tipap=[3,3,3];
%bh=[0.3,0.4;0.4,0.25;0.3,0.4];
%nudosele=[1,2;4,1;3,1];
%%% portico en Y
%xy=[0,3,0;0,3,3;0,0,3];
%ap=[1,3];
%w=[0,35];
%tipap=[3,3];
%bh=[0.3,0.4;0.4,0.25];
%nudosele=[1,2;2,3];
%%% portico en Xclc
%xy=[0,3,0;0,3,3;3,3,3];
%ap=[1,3];
%w=[0,35];
%tipap=[3,3];
%bh=[0.4,0.3;0.4,0.25];
%nudosele=[1,2;2,3];
sizexy=size(xy);
nudos=1:sizexy(1,1);
GL=nudosele*6-5;
elemetos=size(nudosele);
xye=zeros(elemetos(1,1),6);
L=zeros(elemetos(1,1),1);
cx=zeros(elemetos(1,1),1);
cy=zeros(elemetos(1,1),1);
cz=zeros(elemetos(1,1),1);
kgt=zeros(sizexy(1,1)*6);
mgt=zeros(sizexy(1,1)*6,sizexy(1,1)*6);
femt=zeros(sizexy(1,1)*6,1);
ke=zeros(elemetos(1,1)*12,12);
Te=zeros(elemetos(1,1)*12,12);
Femt=zeros(elemetos(1,1)*12,1);
recorrido=1:12:elemetos(1,1)*12;
GLele=zeros(elemetos(1,1)*12,1);
glm=1:(sizexy(1,1)*6);
E = 22*1000^2;
%np=[
G=8.5*1000^2;           
ks=5/6;
reacciones=diag(ones(sizexy(1,1)*6,1));
for i=1:elemetos(1,1)
    xye(i,:)=[xy(nudosele(i,1),1),xy(nudosele(i,2),1),xy(nudosele(i,1),2),xy(nudosele(i,2),2),xy(nudosele(i,1),3),xy(nudosele(i,2),3)];
    L(i,1)=((xye(i,2)-xye(i,1))^2+(xye(i,4)-xye(i,3))^2+(xye(i,5)-xye(i,6))^2)^(1/2);
    cx(i,1)=(xye(i,2)-xye(i,1))/L(i);
    cz(i,1)=(xye(i,4)-xye(i,3))/L(i);%%%Y
    cy(i,1)=(xye(i,6)-xye(i,5))/L(i);%%%Z
    Ra=w(i)*L(i)/2;
    Ma=w(i)*L(i)^2/12;
    Rb=w(i)*L(i)/2;
    Mb=w(i)*L(i)^2/12;
    Pa=w(i)*L(i)/2;
    Pb=w(i)*L(i)/2;
    cxz=sqrt(cx.^2+cz.^2);
    kel=zeros(12);
    k=zeros(12);
    T=zeros(12);
    C=zeros(12);
    mc=zeros(12);
    Fem=zeros(12,1);
    Iz=bh(i,1).*bh(i,2).^3/12;
    Iy=bh(i,2).*bh(i,1).^3/12;
    A=bh(i,2).*bh(i,1); 
    pz=12*E*Iz/(ks*A*G*L(i)^2);
    z12=12*E*Iz/(L(i)^3*(1+pz));
    z6=6*E*Iz/(L(i)^2*(1+pz));
    z4=(4+pz)*E*Iz/(L(i)*(1+pz));
    z2=(2-pz)*E*Iz/(L(i)*(1+pz));
    py=12*E*Iy/(ks*A*G*L(i)^2);
    y12=12*E*Iy/(L(i)^3*(1+py));
    y6=6*E*Iy/(L(i)^2*(1+py));
    y4=(4+py)*E*Iy/(L(i)*(1+py));
    y2=(2-py)*E*Iy/(L(i)*(1+py));
    AE=A*E;
    bj=min([bh(i,1),bh(i,2)]);
    hj=max([bh(i,1),bh(i,2)]);
    J=(1/3-0.21*bj/hj*(1-1/12*(bj/hj)^4))*hj*bj^3;
    GJ=J*G;
    gama=w(i)/A;
   % matriz de rigidez local expresada en el sistema de coordenadas locales
   % para la barra e
    k(:,:)= [AE/L(i),0     ,0     ,0       ,0    ,0       ,-AE/L(i),0      ,0    ,0       ,0  ,0  ;
                0   ,z12   ,0     ,0       ,0    ,z6      ,0       ,-z12   ,0    ,0       ,0  ,z6 ;
                0   ,0     ,y12   ,0       ,-y6  ,0       ,0       ,0      ,-y12 ,0       ,-y6,0  ;
                0   ,0     ,0     ,GJ/L(i) ,0    ,0       ,0       ,0      ,0    ,-GJ/L(i),0  ,0  ;
                0   ,0     ,-y6   ,0       ,y4   ,0       ,0       ,0      ,y6   ,0       ,y2 ,0  ;
                0   ,z6    ,0     ,0       ,0    ,z4      ,0       ,-z6    ,0    ,0       ,0  ,z2 ;
            -AE/L(i),0     ,0     ,0       ,0    ,0       ,AE/L(i) ,0      ,0    ,0       ,0  ,0  ;
               0    ,-z12  ,0      ,0       ,0   ,-z6     ,0       ,z12    ,0    ,0       ,0  ,-z6;
               0    ,0     ,-y12   ,0       ,y6  ,0       ,0       ,0      ,y12  ,0       ,y6  ,0 ;
               0    ,0     ,0      ,-GJ/L(i),0   ,0       ,0       ,0      ,0    ,GJ/L(i) ,0  ,0  ;
               0    ,0     ,-y6    ,0       ,y2  ,0       ,0       ,0      ,y6   ,0       ,y4 ,0  ;
               0    ,z6    ,0      ,0       ,0   ,z2      ,0       ,-z6    ,0    ,0       ,0  ,z4];
             
           
        if isnan(-cx(i)*cy(i)/cxz(i))==0
         T(:,:)=[cx(i)            ,cy(i) ,cz(i)              ,0                  ,0     ,0                  ,0                  ,0     ,0                  ,0                  ,0     ,0;
               -cx(i)*cy(i)/cxz(i),cxz(i),-cy(i)*cz(i)/cxz(i),0                  ,0     ,0                  ,0                  ,0     ,0                  ,0                  ,0     ,0;
               -cz(i)/cxz(i)      ,0     ,cx(i)/cxz(i)       ,0                  ,0     ,0                  ,0                  ,0     ,0                  ,0                  ,0     ,0;
                 0                ,0     ,0                  ,cx(i)              ,cy(i) ,cz(i)              ,0                  ,0     ,0                  ,0                  ,0     ,0;                   
                 0                ,0     ,0                  ,-cx(i)*cy(i)/cxz(i),cxz(i),-cy(i)*cz(i)/cxz(i),0                  ,0     ,0                  ,0                  ,0     ,0; 
                 0                ,0     ,0                  ,-cz(i)/cxz(i)      ,0     ,cx(i)/cxz(i)       ,0                  ,0     ,0                  ,0                  ,0     ,0;
                 0                ,0     ,0                  ,0                  ,0     ,0                  ,cx(i)              ,cy(i) ,cz(i)              ,0                  ,0     ,0;
                 0                ,0     ,0                  ,0                  ,0     ,0                  ,-cx(i)*cy(i)/cxz(i),cxz(i),-cy(i)*cz(i)/cxz(i),0                  ,0     ,0;
                 0                ,0     ,0                  ,0                  ,0     ,0                  ,-cz(i)/cxz(i)      ,0     ,cx(i)/cxz(i)       ,0                  ,0     ,0;
                 0                ,0     ,0                  ,0                  ,0     ,0                  ,0                  ,0     ,0                  ,cx(i)              ,cy(i) ,cz(i)              ;
                 0                ,0     ,0                  ,0                  ,0     ,0                  ,0                  ,0     ,0                  ,-cx(i)*cy(i)/cxz(i),cxz(i),-cy(i)*cz(i)/cxz(i);
                 0                ,0     ,0                  ,0                  ,0     ,0                  ,0                  ,0     ,0                  ,-cz(i)/cxz(i)      ,0     ,cx(i)/cxz(i)       ] ; 
        end    
        if cy(i)==1
         T(:,:)=[0,cy(i),0,0,0,0,0,0,0,0,0,0;
                -cy(i)  ,0,0,0,0,0,0,0,0,0,0,0;
                 0   ,0,1,0,0,0,0,0,0,0,0,0;
                 0   ,0,0,0,cy(i),0,0,0,0,0,0,0;
                 0   ,0,0,-cy(i),0,0,0,0,0,0,0,0;
                 0   ,0,0,0,0,1,0,0,0,0,0,0;
                 0   ,0,0,0,0,0,0,cy(i),0,0,0,0;
                 0   ,0,0,0,0,0,-cy(i),0,0,0,0,0;
                 0   ,0,0,0,0,0,0,0,1,0,0,0;
                 0   ,0,0,0,0,0,0,0,0,0,cy(i),0;
                 0   ,0,0,0,0,0,0,0,0,-cy(i),0,0;
                 0   ,0,0,0,0,0,0,0,0,0,0,1];
        end           
C(:,:)=[-Pa,0,0,0,0,0,0,0,0,0,0,0;
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
   
mc(:,:)=gama*L(i)*A/2*[1,0,0,0,0,0,0,0,0,0,0,0;
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
                   
    we=[0;-1;0;0;0;0;0;-1;0;0;0;0];
    kel(:,:)=T'*k*T;    
    Fem(:,:)=T'*C*T*we;
    mc(:)=T'*mc*T;

    kg=zeros(sizexy(1,1)*6);  
    %kg([GL(i,1):GL(i,1)+5,GL(i,2):GL(i,2)+5],[GL(i,1):GL(i,1)+5,GL(i,2):GL(i,2)+5])=k;   
    mg=zeros(sizexy(1,1)*6,sizexy(1,1)*6);
    mgt([GL(i,1):GL(i,1)+5,GL(i,2):GL(i,2)+5],[GL(i,1):GL(i,1)+5,GL(i,2):GL(i,2)+5])=mc;
    kg(GL(i,1):GL(i,1)+5,GL(i,1):GL(i,1)+5)=kel(1:6,1:6);  
    kg(GL(i,1):GL(i,1)+5,GL(i,2):GL(i,2)+5)=kel(1:6,7:12);  
    kg(GL(i,2):GL(i,2)+5,GL(i,1):GL(i,1)+5)=kel(7:12,1:6);  
    kg(GL(i,2):GL(i,2)+5,GL(i,2):GL(i,2)+5)=kel(7:12,7:12); 
    fg=zeros(sizexy(1,1)*6,1);
    fg(GL(i,1):GL(i,1)+5,:)=Fem(1:6,:);  
    fg(GL(i,2):GL(i,2)+5,:)=Fem(7:12,:);   
    kgt(:,:)=kg+kgt(:,:);
    mgt(:,:)=mg(:,:)+mgt(:,:);
    femt(:,:)=fg+femt(:,:);
    ke(recorrido(i):(recorrido(i)+11),:)=kel;
    Te(recorrido(i):(recorrido(i)+11),:)=T;
    Femt(recorrido(i):(recorrido(i)+11),1)=Fem;
    GLele(recorrido(i):(recorrido(i)+11),1)=[GL(i,1):GL(i,1)+5,GL(i,2):GL(i,2)+5]';
    %plot((xye(i,2)+xye(i,1))/2,(xye(i,4)+xye(i,3))/2,'o','MarkerSize',17,...
    %    'MarkerEdgeColor','k',...
     %   'MarkerFaceColor','k'), hold on
     x=[xy(nudosele(i,1),1),xy(nudosele(i,2),1)];
     y=[xy(nudosele(i,1),2),xy(nudosele(i,2),2)];
     z=[xy(nudosele(i,1),3),xy(nudosele(i,2),3)];
     if w(i)>0    
         plot3(x,y,z,'r'), hold on,grid on;
     else
        plot3(x,y,z,'b'), hold on,grid on;
     end
    text(xye(i,1),xye(i,3),xye(i,5),num2str(nudosele(i,1)),'Color','r','HorizontalAlignment','left','FontSize',14), hold on;%% grafica nudos elementos    
    text(xye(i,2),xye(i,4),xye(i,6),num2str(nudosele(i,2)),'Color','r','HorizontalAlignment','left','FontSize',14), hold on;%% grafica nudos elementos    
end
a=1;
set(gca,'DataAspectRatio',[1 1 1],...
        'PlotBoxAspectRatio',[1 1 1]);
val=zeros(sizexy(1,1)*6,1);
sizenp=size(ap);
for i=1:sizenp(1,2)
    if tipap(i)==3
        kgt(:,ap(i)*6-5)=val;
        kgt(:,ap(i)*6-5+1)=val;
        kgt(:,ap(i)*6-5+2)=val;
        kgt(:,ap(i)*6-5+3)=val;
        kgt(:,ap(i)*6-5+4)=val;
        kgt(:,ap(i)*6-5+5)=val;
        kgt(ap(i)*6-5,ap(i)*6-5)=-1;
        kgt(ap(i)*6-5+1,ap(i)*6-5+1)=-1;
        kgt(ap(i)*6-5+2,ap(i)*6-5+2)=-1;
        kgt(ap(i)*6-5+3,ap(i)*6-5+3)=-1;
        kgt(ap(i)*6-5+4,ap(i)*6-5+4)=-1;
        kgt(ap(i)*6-5+5,ap(i)*6-5+5)=-1;
        glm((ap(i)*6-5):(ap(i)*6))=zeros(6,1);
        plot3(xy(ap(i),1),xy(ap(i),2),xy(ap(i),3),'-s','MarkerSize',15,...
        'MarkerEdgeColor','blue',...
        'MarkerFaceColor','blue'), hold on
    end      
end
F=zeros(sizexy(1,1)*2,1);
%F(np*2-1+1,1)=p;
f=find(glm~=0);
glm=f;
kmodal=full(kgt(f,f));%% matriz de rigidez modal kN/m % matriz de rigidez timoshenko
mmodal=full(mgt(f,f));%% matriz de masa modal kN % matriz de masa de timoshenko
RD=kgt\femt;
D=RD;
for i=1:sizenp(1,2)
    if  tipap(i)==3
        D((ap(i)*6-5):(ap(i)*6),1)=zeros(6,1);
    end 
end
def=zeros(elemetos(1,1)*12,1);
femp=zeros(elemetos(1,1)*12,1);
for i=1:elemetos(1,1)
    def(recorrido(i):recorrido(i)+11,:)=D(GLele(recorrido(i):recorrido(i)+11),1);
    femp(recorrido(i):recorrido(i)+11,:)=Te(recorrido(i):recorrido(i)+11,:)*(ke(recorrido(i):recorrido(i)+11,:)*def(recorrido(i):recorrido(i)+11,:)-Femt(recorrido(i):recorrido(i)+11,:));
    text((xye(i,2)+xye(i,1))/2,(xye(i,4)+xye(i,3))/2,(xye(i,5)+xye(i,6))/2,num2str(i),'Color','M','HorizontalAlignment','center'), hold on;%% elemetos 
end    