function [A,Ix,Iy,As2,As3,J]=propiedades_geometrica_perfiles_v1(d,bf,tw,tf)
%https://calcresource.com/cross-section-doubletee-unsym.html
%d=300;
%bf=150;
%tw=6;
%tf=9;
%Channel
%Rectangular
%Circle
%I/Wide Flange
%Pipe             Seccion circular hueca
%Box/Tube         Tubular hueca
%Angle            Perfil L 
%Tee              Perfil T
A=[bf*d-(bf-tw)*(d-tf*2);
   d*bf; 
   d^2*pi/4;
   bf*d-(bf-tw)*(d-tf*2);
   d^2*pi/4-(d-tf*2)^2*pi/4;
   bf*d-(bf-tw*2)*(d-tf*2);
   d*tf+tf*(d-tf);
   bf*tf+(d-tf)*tw];
Ix=[bf*d^3/12-(bf-tw)*(d-tf*2)^3/12;
    bf*d^3/12;
    1/4*pi*(d/2)^4;
    bf*d^3/12-(bf-tw)*(d-tf*2)^3/12;
    1/4*pi*(d/2)^4-1/4*pi*(d/2-tf)^4;
    bf*d^3/12-(bf-tw*2)*(d-tf*2)^3/12;
   (tf*(5*d^4 - 10*d^3*tf + 11*d^2*tf^2 - 6*d*tf^3 + tf^4))/(12*(2*d - tf));
   bf*tf^3/12+bf*tf*((d-tf/2)-((bf*tf*(d-tf/2)+(d-tf)*tw*(d-tf)/2)/(bf*tf+(d-tf)*tw)))^2+tw*(d-tf)^3/12+tw*(d-tf)*((d-tf)/2-((bf*tf*(d-tf/2)+(d-tf)*tw*(d-tf)/2)/(bf*tf+(d-tf)*tw)))^2];
Iy=[(bf^3*tf)/6 + (tw^3*(d - 2*tf))/12 + 2*bf*tf*(bf/2 - (bf^2*tf + (tw^2*(d - 2*tf))/2)/(2*bf*tf + tw*(d - 2*tf)))^2 + tw*(d - 2*tf)*(tw/2 - (bf^2*tf + (tw^2*(d - 2*tf))/2)/(2*bf*tf + tw*(d - 2*tf)))^2;
    d*bf^3/12;
    1/4*pi*(d/2)^4;
    tf*bf^3/12*2+(d-tf*2)*tw^3/12;
    1/4*pi*(d/2)^4-1/4*pi*(d/2-tf)^4;
    bf^3*d/12-(bf-tw*2)^3*(d-tf*2)/12;
    bf^3*tf/12+tw^3*(d-tf)/12;
    (tf*(5*d^4 - 10*d^3*tf + 11*d^2*tf^2 - 6*d*tf^3 + tf^4))/(12*(2*d - tf))];
As3=[bf*tf*2;
    bf*d*5/6;
    d^2*pi/4*9/10;
    5/3*bf*tf;
    (pi*(2*d/2-tf)*tf)*(0.9-0.4*(d/2-tf)/(d/2));
    2*tf*bf;
    d*tf;
    bf*tf*5/6];
As2=[d*tw;
    bf*d*5/6;
    d^2*pi/4*9/10;
    tw*d;
    (pi*(2*d/2-tf)*tf)*(0.9-0.4*(d/2-tf)/(d/2));
    d*tw*2;
    d*tf;
    d*tw];
bj=min(bf,d);
dj=max(bf,d);
Jrec=(1/3-0.21*bj/dj*(1-1/12*(bj/dj)^4))*dj*bj^3;
J=[1/3*((d-tf)*tw^3+2*(bf-tw/2)*tf^3);
   Jrec;
   1/2*pi*(d/2)^4;
   1/3*((d-tf*2)*tw^3+2*(bf-tw)*tf^3);
   1/2*pi*((d/2)^4-((d-tf*2)/2)^4);
   2*((bf-tw)*(d-tf))^2/((bf-tw)/tf+(d-tf)/tw);
   1/3*((d-tf/2)*tf^3+(bf-tf/2)*tf^3);
   1/3*((d-tf/2)*tw^3+bf*tf^3)];
%Ai=A(i);
%Ixi=Ix(i);
%As3i=As3(i);
       