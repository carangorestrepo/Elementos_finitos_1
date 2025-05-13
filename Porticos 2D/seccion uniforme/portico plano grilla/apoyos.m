function [x,y,xn,yn]=apoyos(xr,yr,f,t)
%xr coordenada x apoyo
%yr coordenada y apoyo
%f escala

if t==12 %% apoyado
    %xr=0;
    %yr=-f;
    %f=0.3;
    n = 3;
    phi=360/(n*2)*pi/180;
    R = (f/2)/(cos(phi));
    tita = (0:(2*pi/n):2*pi)+pi/2; 
    x = R*cos(tita)+xr;
    y = R*sin(tita)+yr-f;
    xn=0;
    yn=0;
elseif t==2 %% apoyado vertical
%plot(x,y)
    ang=(0:360)*pi/180;
    x1=cos(ang)*f/2+xr;
    y1=sin(ang)*f/2+yr;
    x=x1;
    y=y1-f/2;
    xn=[-f/2,f/2]+xr;
    yn=[-f,-f]+yr;
elseif t==1 %% apoyado horizontal
    ang=(0:360)*pi/180;
    x1=cos(ang)*f/2+xr;
    y1=sin(ang)*f/2+yr;
    x=x1+f/2;
    y=y1;
    xn=[0,0]+xr;
    yn=[f/2,-f/2]+yr; 
elseif t==123% empotrado
    x=[-f/2,-f/2,f/2,f/2,-f/2]+xr;
    y=[-f,0,0,-f,-f]+yr;
    xn=0;
    yn=0;
 end

%plot(x,y,xn,yn)