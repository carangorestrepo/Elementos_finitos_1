function [Axh,Ixh,As2h]=secciones(sec,b1,b2,nb,h1,h2,nh,q1,q2,nq,G,E,tf,tw,a,b,L)
syms x
bfx=(b1 - (x^nb*(b1 - b2))/L^nb);
dx =(h1 - (x^nh*(h1 - h2))/L^nh);
%seccion C ==1
%Seccion rectangular  ==2
%Circular solida ==3 
%Seccion I ==4
%Seccion circular hueca ==5
%Tubular hueca ==6

if sec==1 
    %seccion C
    Ax=bfx*dx-(bfx-tw)*(dx-tf*2);
    Ix=bfx*dx^3/12-(bfx-tw)*(dx-tf*2)^3/12;
    As2x=dx*tw*G;
elseif sec==2
%Seccion rectangular 
    Ax=dx*bfx;
    Ix=bfx*dx^3/12;
    As2x=bfx*dx*5/6*G;
elseif sec==3
%Circular solida
    Ax=dx^2*pi/4;
    Ix=1/4*pi*(dx/2)^4;
    As2x=dx^2*pi/4*9/10*G;
elseif sec==4
    %Seccion I
    Ax=bfx*dx-(bfx-tw)*(dx-tf*2);
    Ix=bfx*dx^3/12-(bfx-tw)*(dx-tf*2)^3/12;
    As2x=tw*dx*G;
elseif sec==5
%Seccion circular hueca
    Ax=dx^2*pi/4-(dx-tf*2)^2*pi/4;
    Ix=1/4*pi*(dx/2)^4-1/4*pi*(dx/2-tf)^4;
    As2x=(pi*(2*dx/2-tf)*tf)*(0.9-0.4*(dx/2-tf)/(dx/2))*G;
elseif sec==6
    %Tubular hueca
    Ax=bfx*dx-(bfx-tw*2)*(dx-tf*2);
    Ix=bfx*dx^3/12-(bfx-tw*2)*(dx-tf*2)^3/12;
    As2x=dx*tw*2*G;
elseif sec==7
    %Perfil T
    Ax=bfx*tf+(dx-tf)*tw;
    Ix=bfx*tf^3/12+bfx*tf*((dx-tf/2)-((bfx*tf*(dx-tf/2)+(dx-tf)*tw*(dx-tf)/2)/(bfx*tf+(dx-tf)*tw)))^2+tw*(dx-tf)^3/12+tw*(dx-tf)*((dx-tf)/2-((bfx*tf*(dx-tf/2)+(dx-tf)*tw*(dx-tf)/2)/(bfx*tf+(dx-tf)*tw)))^2;
    As2x=dx*tw*2*G;
end
Axh=matlabFunction(Ax,"Vars",{x});
Ixh=matlabFunction(Ix,"Vars",{x});
As2h=matlabFunction(As2x,"Vars",{x});