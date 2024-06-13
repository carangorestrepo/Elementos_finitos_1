function [Kloc,FE,Me,V,M,fax,v,u]=matriz_rigidez(tipo_conti,EA,Ac,EI,L,b1,b2,q1,q2,ana,v1,v2,t1,t2,u1,u2,x,k,P,lalb,puntos_graficas)

% Se definen algunas iniciales
V=0;
M=0;
v=0;
fax=0;
u=0;
Kloc=0;
FE=0;
Me=0;
%syms x  
%% caso EE
if strcmp(tipo_conti,'EE')==1
   if ana==0 
        % matriz de rigidez local expresada en el sistema de coordenadas locales
        % empotrado empotrado
        Kloc =[[  EA/L,                                0,                                              0, -EA/L,                                0,                                              0];
              [     0,  (12.*Ac.*EI)/(L.*(Ac.*L.^2 + 12.*EI)),                     (6.*Ac.*EI)/(Ac.*L.^2 + 12.*EI),     0, -(12.*Ac.*EI)/(L.*(Ac.*L.^2 + 12.*EI)),                     (6.*Ac.*EI)/(Ac.*L.^2 + 12.*EI)];
              [     0,       (6.*Ac.*EI)/(Ac.*L.^2 + 12.*EI),    (4.*EI.*(Ac.*L.^2 + 3.*EI))/(L.*(Ac.*L.^2 + 12.*EI)),     0,      -(6.*Ac.*EI)/(Ac.*L.^2 + 12.*EI), -(2.*EI.*(- Ac.*L.^2 + 6.*EI))/(L.*(Ac.*L.^2 + 12.*EI))];
              [ -EA/L,                                0,                                              0,  EA/L,                                0,                                              0];
              [     0, -(12.*Ac.*EI)/(L.*(Ac.*L.^2 + 12.*EI)),                    -(6.*Ac.*EI)/(Ac.*L.^2 + 12.*EI),     0,  (12.*Ac.*EI)/(L.*(Ac.*L.^2 + 12.*EI)),                    -(6.*Ac.*EI)/(Ac.*L.^2 + 12.*EI)];
              [     0,       (6.*Ac.*EI)/(Ac.*L.^2 + 12.*EI), -(2.*EI.*(- Ac.*L.^2 + 6.*EI))/(L.*(Ac.*L.^2 + 12.*EI)),     0,      -(6.*Ac.*EI)/(Ac.*L.^2 + 12.*EI),    (4.*EI.*(Ac.*L.^2 + 3.*EI))/(L.*(Ac.*L.^2 + 12.*EI))]];
       % momentos y reaccciones de empotramiento
        x=0;
        X1 = ((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L);
        Y1 = q1.*x - (240.*Ac.*EI.*v2 - 240.*Ac.*EI.*v1 + 7.*Ac.*L.^4.*q1 + 3.*Ac.*L.^4.*q2 + 80.*EI.*L.^2.*q1 + 40.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t1 - 120.*Ac.*EI.*L.*t2)/(20.*L.*(Ac.*L.^2 + 12.*EI)) - (x.^2.*(q1 - q2))/(2.*L);
        M1 =(720.*EI.^2.*t2 - 720.*EI.^2.*t1 + 3.*Ac.*L.^5.*q1 + 2.*Ac.*L.^5.*q2 + 30.*EI.*L.^3.*q1 + 30.*EI.*L.^3.*q2 - 240.*Ac.*EI.*L.^2.*t1 - 120.*Ac.*EI.*L.^2.*t2 - 360.*Ac.*EI.*L.*v1 + 360.*Ac.*EI.*L.*v2)/(60.*L.*(Ac.*L.^2 + 12.*EI)) - (q1.*x.^2)/2 - x.*(((q1 - q2).*x.^2)/(2.*L) - q1.*x + (240.*Ac.*EI.*v2 - 240.*Ac.*EI.*v1 + 7.*Ac.*L.^4.*q1 + 3.*Ac.*L.^4.*q2 + 80.*EI.*L.^2.*q1 + 40.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t1 - 120.*Ac.*EI.*L.*t2)/(20.*L.*(Ac.*L.^2 + 12.*EI))) + (x.^3.*(q1 - q2))/(3.*L); 
        x=L;
        X2 = ((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L);
        Y2 = q1.*x - (240.*Ac.*EI.*v2 - 240.*Ac.*EI.*v1 + 7.*Ac.*L.^4.*q1 + 3.*Ac.*L.^4.*q2 + 80.*EI.*L.^2.*q1 + 40.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t1 - 120.*Ac.*EI.*L.*t2)/(20.*L.*(Ac.*L.^2 + 12.*EI)) - (x.^2.*(q1 - q2))/(2.*L);
        M2 = (720.*EI.^2.*t2 - 720.*EI.^2.*t1 + 3.*Ac.*L.^5.*q1 + 2.*Ac.*L.^5.*q2 + 30.*EI.*L.^3.*q1 + 30.*EI.*L.^3.*q2 - 240.*Ac.*EI.*L.^2.*t1 - 120.*Ac.*EI.*L.^2.*t2 - 360.*Ac.*EI.*L.*v1 + 360.*Ac.*EI.*L.*v2)/(60.*L.*(Ac.*L.^2 + 12.*EI)) - (q1.*x.^2)/2 - x.*(((q1 - q2).*x.^2)/(2.*L) - q1.*x + (240.*Ac.*EI.*v2 - 240.*Ac.*EI.*v1 + 7.*Ac.*L.^4.*q1 + 3.*Ac.*L.^4.*q2 + 80.*EI.*L.^2.*q1 + 40.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t1 - 120.*Ac.*EI.*L.*t2)/(20.*L.*(Ac.*L.^2 + 12.*EI))) + (x.^3.*(q1 - q2))/(3.*L); 
        FE=[X1;-Y1;M1;X2;Y2;-M2];
    elseif ana==1
        % momentos, cortantes,axiales y deformadas
        V =q1.*x - (240.*Ac.*EI.*v2 - 240.*Ac.*EI.*v1 + 7.*Ac.*L.^4.*q1 + 3.*Ac.*L.^4.*q2 + 80.*EI.*L.^2.*q1 + 40.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t1 - 120.*Ac.*EI.*L.*t2)/(20.*L.*(Ac.*L.^2 + 12.*EI)) - (x.^2.*(q1 - q2))/(2.*L);
        M =(720.*EI.^2.*t2 - 720.*EI.^2.*t1 + 3.*Ac.*L.^5.*q1 + 2.*Ac.*L.^5.*q2 + 30.*EI.*L.^3.*q1 + 30.*EI.*L.^3.*q2 - 240.*Ac.*EI.*L.^2.*t1 - 120.*Ac.*EI.*L.^2.*t2 - 360.*Ac.*EI.*L.*v1 + 360.*Ac.*EI.*L.*v2)/(60.*L.*(Ac.*L.^2 + 12.*EI)) - (q1.*x.^2)/2 - x.*(((q1 - q2).*x.^2)/(2.*L) - q1.*x + (240.*Ac.*EI.*v2 - 240.*Ac.*EI.*v1 + 7.*Ac.*L.^4.*q1 + 3.*Ac.*L.^4.*q2 + 80.*EI.*L.^2.*q1 + 40.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t1 - 120.*Ac.*EI.*L.*t2)/(20.*L.*(Ac.*L.^2 + 12.*EI))) + (x.^3.*(q1 - q2))/(3.*L); 
        v =(((q1 - q2).*x.^5)/(30.*L) - (q1.*x.^4)/24 + (EI.*(720.*EI.^2.*t2 - 720.*EI.^2.*t1 + 60.*Ac.^2.*L.^3.*v1 + 3.*Ac.*L.^5.*q1 + 2.*Ac.*L.^5.*q2 + 30.*EI.*L.^3.*q1 + 30.*EI.*L.^3.*q2 - 240.*Ac.*EI.*L.^2.*t1 - 120.*Ac.*EI.*L.^2.*t2 + 360.*Ac.*EI.*L.*v1 + 360.*Ac.*EI.*L.*v2))/(60.*Ac.*L.*(Ac.*L.^2 + 12.*EI)))/EI - (1/Ac - x.^2/(2.*EI)).*(((q1 - q2).*x.^3)/(3.*L) - (q1.*x.^2)/2 + (720.*EI.^2.*t2 - 720.*EI.^2.*t1 + 3.*Ac.*L.^5.*q1 + 2.*Ac.*L.^5.*q2 + 30.*EI.*L.^3.*q1 + 30.*EI.*L.^3.*q2 - 240.*Ac.*EI.*L.^2.*t1 - 120.*Ac.*EI.*L.^2.*t2 - 360.*Ac.*EI.*L.*v1 + 360.*Ac.*EI.*L.*v2)/(60.*L.*(Ac.*L.^2 + 12.*EI))) + (- x.^3/(6.*EI) + x/Ac).*(((q1 - q2).*x.^2)/(2.*L) - q1.*x + (240.*Ac.*EI.*v2 - 240.*Ac.*EI.*v1 + 7.*Ac.*L.^4.*q1 + 3.*Ac.*L.^4.*q2 + 80.*EI.*L.^2.*q1 + 40.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t1 - 120.*Ac.*EI.*L.*t2)/(20.*L.*(Ac.*L.^2 + 12.*EI))) + (x.*(- ((q1 - q2).*x.^4)/(8.*L) + (q1.*x.^3)/6 + EI.*t1))/EI;
        u =(- ((b1 - b2).*x.^3)/(3.*L) + (b1.*x.^2)/2 + EA.*u1)/EA + (x.*(((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L)))/EA;
        fax =((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L);
   
   end
   
%% caso RE
elseif strcmp(tipo_conti,'RE')==1
    if ana==0
        % matriz de rigidez local expresada en el sistema de coordenadas locales
        % articulado empotrado
        Kloc = [[  EA/L,                                    0,0,-EA/L,                              0,                           0];
                [     0,  (3.*Ac.*EI)/(L.*(Ac.*L.^2 + 3.*EI)),0,0    , -(3.*Ac.*EI)/(L.*(Ac.*L.^2 + 3.*EI)),   (3.*Ac.*EI)/(Ac.*L.^2 + 3.*EI)];
                [     0,                                    0,0,0    ,                              0,                           0];
                [ -EA/L,                                    0,0,EA/L ,                              0,                           0];
                [     0, -(3.*Ac.*EI)/(L.*(Ac.*L.^2 + 3.*EI)),0,    0,  (3.*Ac.*EI)/(L.*(Ac.*L.^2 + 3.*EI)),  -(3.*Ac.*EI)/(Ac.*L.^2 + 3.*EI)];
                [     0,       (3.*Ac.*EI)/(Ac.*L.^2 + 3.*EI),0,    0,     -(3.*Ac.*EI)/(Ac.*L.^2 + 3.*EI), (3.*Ac.*EI.*L)/(Ac.*L.^2 + 3.*EI)]];   
        % momentos y reaccciones de empotramiento
        x=0;
        X1 = ((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L);
        Y1 = q1.*x - (120.*Ac.*EI.*v2 - 120.*Ac.*EI.*v1 + 11.*Ac.*L.^4.*q1 + 4.*Ac.*L.^4.*q2 + 40.*EI.*L.^2.*q1 + 20.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t2)/(40.*L.*(Ac.*L.^2 + 3.*EI)) - (x.^2.*(q1 - q2))/(2.*L);
        M1 = 0;
        x=L;
        X2 = ((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L);
        Y2 = q1.*x - (120.*Ac.*EI.*v2 - 120.*Ac.*EI.*v1 + 11.*Ac.*L.^4.*q1 + 4.*Ac.*L.^4.*q2 + 40.*EI.*L.^2.*q1 + 20.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t2)/(40.*L.*(Ac.*L.^2 + 3.*EI)) - (x.^2.*(q1 - q2))/(2.*L);
        M2 = (x.^3.*(q1 - q2))/(3.*L) - (q1.*x.^2)/2 - x.*(((q1 - q2).*x.^2)/(2.*L) - q1.*x + (120.*Ac.*EI.*v2 - 120.*Ac.*EI.*v1 + 11.*Ac.*L.^4.*q1 + 4.*Ac.*L.^4.*q2 + 40.*EI.*L.^2.*q1 + 20.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t2)/(40.*L.*(Ac.*L.^2 + 3.*EI)));  
        FE=[X1;Y1;M1;X2;Y2;M2];
    elseif ana==1 
        % cortantes, momentos, deformadas verticales, deformadas axiales, fuerzas axiales 
        V =q1.*x - (120.*Ac.*EI.*v2 - 120.*Ac.*EI.*v1 + 11.*Ac.*L.^4.*q1 + 4.*Ac.*L.^4.*q2 + 40.*EI.*L.^2.*q1 + 20.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t2)/(40.*L.*(Ac.*L.^2 + 3.*EI)) - (x.^2.*(q1 - q2))/(2.*L);
        M =(x.^3.*(q1 - q2))/(3.*L) - (q1.*x.^2)/2 - x.*(((q1 - q2).*x.^2)/(2.*L) - q1.*x + (120.*Ac.*EI.*v2 - 120.*Ac.*EI.*v1 + 11.*Ac.*L.^4.*q1 + 4.*Ac.*L.^4.*q2 + 40.*EI.*L.^2.*q1 + 20.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t2)/(40.*L.*(Ac.*L.^2 + 3.*EI))) ;
        v =(((q1 - q2).*x.^5)/(30.*L) - (q1.*x.^4)/24 + EI.*v1)/EI + (- ((q1 - q2).*x.^3)/(3.*L) + (q1.*x.^2)/2).*(1/Ac - x.^2/(2.*EI)) + (- x.^3/(6.*EI) + x/Ac).*(((q1 - q2).*x.^2)/(2.*L) - q1.*x + (120.*Ac.*EI.*v2 - 120.*Ac.*EI.*v1 + 11.*Ac.*L.^4.*q1 + 4.*Ac.*L.^4.*q2 + 40.*EI.*L.^2.*q1 + 20.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t2)/(40.*L.*(Ac.*L.^2 + 3.*EI))) + (x.*(- ((q1 - q2).*x.^4)/(8.*L) + (q1.*x.^3)/6 + (720.*EI.^2.*t2 + 3.*Ac.*L.^5.*q1 + 2.*Ac.*L.^5.*q2 + 30.*EI.*L.^3.*q1 + 30.*EI.*L.^3.*q2 - 120.*Ac.*EI.*L.^2.*t2 - 360.*Ac.*EI.*L.*v1 + 360.*Ac.*EI.*L.*v2)/(240.*Ac.*L.^2 + 720.*EI)))/EI;
        u =(- ((b1 - b2).*x.^3)/(3.*L) + (b1.*x.^2)/2 + EA.*u1)/EA + (x.*(((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L)))/EA;
        fax =((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L);
 
    end
%% caso ER
elseif strcmp(tipo_conti,'ER')==1 
    if ana==0
        % matriz de rigidez local expresada en el sistema de coordenadas locales
        % empotrado articulado 
        Kloc = [[  EA/L,                              0,                           0, -EA/L,                              0, 0];
                [     0,  (3.*Ac.*EI)/(L.*(Ac.*L.^2 + 3.*EI)),   (3.*Ac.*EI)/(Ac.*L.^2 + 3.*EI),     0, -(3.*Ac.*EI)/(L.*(Ac.*L.^2 + 3.*EI)), 0];
                [     0,      (3.*Ac.*EI)/(Ac.*L.^2 + 3.*EI), (3.*Ac.*EI.*L)/(Ac.*L.^2 + 3.*EI),     0,     -(3.*Ac.*EI)/(Ac.*L.^2 + 3.*EI), 0];
                [ -EA/L,                              0,                           0,  EA/L,                              0, 0];
                [     0, -(3.*Ac.*EI)/(L.*(Ac.*L.^2 + 3.*EI)),  -(3.*Ac.*EI)/(Ac.*L.^2 + 3.*EI),     0,  (3.*Ac.*EI)/(L.*(Ac.*L.^2 + 3.*EI)), 0];
                [     0,                              0,                           0,     0,                              0, 0]];
       % momentos y reaccciones de empotramiento
        x=0;
        X1 = ((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L);
        Y1 = q1.*x - (120.*Ac.*EI.*v2 - 120.*Ac.*EI.*v1 + 16.*Ac.*L.^4.*q1 + 9.*Ac.*L.^4.*q2 + 40.*EI.*L.^2.*q1 + 20.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t2)/(40.*L.*(Ac.*L.^2 + 3.*EI)) - (x.^2.*(q1 - q2))/(2.*L);      
        M1 = (360.*Ac.*EI.*v2 - 360.*Ac.*EI.*v1 + 8.*Ac.*L.^4.*q1 + 7.*Ac.*L.^4.*q2 - 360.*Ac.*EI.*L.*t2)/(120.*(Ac.*L.^2 + 3.*EI)) - x.*(((q1 - q2).*x.^2)/(2.*L) - q1.*x + (120.*Ac.*EI.*v2 - 120.*Ac.*EI.*v1 + 16.*Ac.*L.^4.*q1 + 9.*Ac.*L.^4.*q2 + 40.*EI.*L.^2.*q1 + 20.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t2)/(40.*L.*(Ac.*L.^2 + 3.*EI))) - (q1.*x.^2)/2 + (x.^3.*(q1 - q2))/(3.*L);
        x=L;
        X2 = ((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L);
        Y2 = q1.*x - (120.*Ac.*EI.*v2 - 120.*Ac.*EI.*v1 + 16.*Ac.*L.^4.*q1 + 9.*Ac.*L.^4.*q2 + 40.*EI.*L.^2.*q1 + 20.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t2)/(40.*L.*(Ac.*L.^2 + 3.*EI)) - (x.^2.*(q1 - q2))/(2.*L);
        M2 = 0;
        FE=[X1;Y1;M1;X2;Y2;M2];
    elseif ana==1
        % cortantes, momentos, deformadas verticales, deformadas axiales, fuerzas axiales 
        V =q1.*x - (120.*Ac.*EI.*v2 - 120.*Ac.*EI.*v1 + 16.*Ac.*L.^4.*q1 + 9.*Ac.*L.^4.*q2 + 40.*EI.*L.^2.*q1 + 20.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t2)/(40.*L.*(Ac.*L.^2 + 3.*EI)) - (x.^2.*(q1 - q2))/(2.*L);
        M =(360.*Ac.*EI.*v2 - 360.*Ac.*EI.*v1 + 8.*Ac.*L.^4.*q1 + 7.*Ac.*L.^4.*q2 - 360.*Ac.*EI.*L.*t2)/(120.*(Ac.*L.^2 + 3.*EI)) - x.*(((q1 - q2).*x.^2)/(2.*L) - q1.*x + (120.*Ac.*EI.*v2 - 120.*Ac.*EI.*v1 + 16.*Ac.*L.^4.*q1 + 9.*Ac.*L.^4.*q2 + 40.*EI.*L.^2.*q1 + 20.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t2)/(40.*L.*(Ac.*L.^2 + 3.*EI))) - (q1.*x.^2)/2 + (x.^3.*(q1 - q2))/(3.*L);
        v =(((q1 - q2).*x.^5)/(30.*L) - (q1.*x.^4)/24 + (EI.*(360.*EI.*v2 + 8.*L.^4.*q1 + 7.*L.^4.*q2 - 360.*EI.*L.*t2 + 120.*Ac.*L.^2.*v1))/(120.*Ac.*L.^2 + 360.*EI))/EI - (1/Ac - x.^2/(2.*EI)).*(((q1 - q2).*x.^3)/(3.*L) - (q1.*x.^2)/2 + (360.*Ac.*EI.*v2 - 360.*Ac.*EI.*v1 + 8.*Ac.*L.^4.*q1 + 7.*Ac.*L.^4.*q2 - 360.*Ac.*EI.*L.*t2)/(120.*Ac.*L.^2 + 360.*EI)) + (- x.^3/(6.*EI) + x/Ac).*(((q1 - q2).*x.^2)/(2.*L) - q1.*x + (120.*Ac.*EI.*v2 - 120.*Ac.*EI.*v1 + 16.*Ac.*L.^4.*q1 + 9.*Ac.*L.^4.*q2 + 40.*EI.*L.^2.*q1 + 20.*EI.*L.^2.*q2 - 120.*Ac.*EI.*L.*t2)/(40.*L.*(Ac.*L.^2 + 3.*EI))) + (x.*(- ((q1 - q2).*x.^4)/(8.*L) + (q1.*x.^3)/6 + EI.*t2))/EI; 
        u =(- ((b1 - b2).*x.^3)/(3.*L) + (b1.*x.^2)/2 + EA.*u1)/EA + (x.*(((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L)))/EA;
        fax =((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L);
    end
elseif strcmp(tipo_conti,'RR')==1
   if ana==0
        % matriz de rigidez local expresada en el sistema de coordenadas locales
        % articulado  articulado 
        Kloc =[[  EA/L,0,0,-EA/L,0,0];
              [     0,0,0,0    ,0,0];
              [     0,0,0,0    ,0,0];
              [ -EA/L,0,0,EA/L ,0,0];
              [     0,0,0,0    ,0,0];
              [     0,0,0,0    ,0,0]];
        % momentos y reaccciones de empotramiento
        x=0;
        X1 = ((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L);
        Y1 = q1.*x - (L.*q2)/6 - (L.*q1)/3 - (x.^2.*(q1 - q2))/(2.*L);
        M1 = 0;
        x=L;
        X2 = ((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L);
        Y2 = q1.*x - (L.*q2)/6 - (L.*q1)/3 - (x.^2.*(q1 - q2))/(2.*L);
        M2 = 0;
        FE=[X1;Y1;M1;X2;Y2;M2];
   elseif ana==1 
       % cortantes, momentos, deformadas verticales, deformadas axiales, fuerzas axiales 
        V =q1.*x - (L.*q2)/6 - (L.*q1)/3 - (x.^2.*(q1 - q2))/(2.*L);
        M = (x.^3.*(q1 - q2))/(3.*L) - (q1.*x.^2)/2 - x.*(((q1 - q2).*x.^2)/(2.*L) - q1.*x + (L.*q1)/3 + (L.*q2)/6); 
        v =(- x.^3/(6.*EI) + x/Ac).*(((q1 - q2).*x.^2)/(2.*L) - q1.*x + (L.*q1)/3 + (L.*q2)/6) + (((q1 - q2).*x.^5)/(30.*L) - (q1.*x.^4)/24 + EI.*v1)/EI + (- ((q1 - q2).*x.^3)/(3.*L) + (q1.*x.^2)/2).*(1/Ac - x.^2/(2.*EI)) + (x.*(- ((q1 - q2).*x.^4)/(8.*L) + (q1.*x.^3)/6 + (360.*EI.*v2 - 360.*EI.*v1 + 8.*L.^4.*q1 + 7.*L.^4.*q2)/(360.*L)))/EI; 
        u =(- ((b1 - b2).*x.^3)/(3.*L) + (b1.*x.^2)/2 + EA.*u1)/EA + (x.*(((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L)))/EA;
        fax =((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L);
   end
elseif strcmp(tipo_conti,'kw')==1
   if ana==0
        % matriz de rigidez local expresada en el sistema de coordenadas locales
        % viga wincler
        P=0;
        Kloc = Ke_ec_dif(L, EI, EA,Ac,k,P,puntos_graficas);
        % momentos y reaccciones de empotramiento
        [X1,Y1,M1,X2,Y2,M2] = fe_ec_dif(L, EI, Ac, EA, q1,q2,q1,q2,k,P,puntos_graficas);
        FE=[X1;Y1;M1;X2;Y2;M2];
        X1=abs(X1);
        Y1=abs(Y1);
        M1=-abs(M1);
        M2=-abs(M2);
   elseif ana==1 
       % cortantes, momentos, deformadas verticales, deformadas axiales, fuerzas axiales 
       [fax,V,M,u,v] = def_ec_dif(L, EI, Ac, EA, q1,q1, q1,q2,k,P,v1,v2,t1,t2,u1,u2,puntos_graficas);  
   end
elseif strcmp(tipo_conti,'NR')==1%%nudo rigido
   if ana==0
        % matriz de rigidez local expresada en el sistema de coordenadas locales
        % nudos rigidos viga continua
        la=lalb(1);
        lb=lalb(2);
        Kloc =[[  EA/L,                                                                                                                  0,                                                                                                                                                                          0, -EA/L,                                                                                                                 0,                                                                                                                                                                          0];
               [     0,              -(12*Ac*EI)/((la - L + lb)*(Ac*L^2 - 2*Ac*L*la - 2*Ac*L*lb + Ac*la^2 + 2*Ac*la*lb + Ac*lb^2 + 12*EI)),                                                         -(6*Ac*EI*(L + la - lb))/((la - L + lb)*(Ac*L^2 - 2*Ac*L*la - 2*Ac*L*lb + Ac*la^2 + 2*Ac*la*lb + Ac*lb^2 + 12*EI)),     0,              (12*Ac*EI)/((la - L + lb)*(Ac*L^2 - 2*Ac*L*la - 2*Ac*L*lb + Ac*la^2 + 2*Ac*la*lb + Ac*lb^2 + 12*EI)),                                                         -(6*Ac*EI*(L - la + lb))/((la - L + lb)*(Ac*L^2 - 2*Ac*L*la - 2*Ac*L*lb + Ac*la^2 + 2*Ac*la*lb + Ac*lb^2 + 12*EI))];
               [     0, -(6*Ac*EI*(L + la - lb))/((la - L + lb)*(Ac*L^2 - 2*Ac*L*la - 2*Ac*L*lb + Ac*la^2 + 2*Ac*la*lb + Ac*lb^2 + 12*EI)),     -(4*EI*(Ac*L^2 + Ac*L*la - 2*Ac*L*lb + Ac*la^2 - Ac*la*lb + Ac*lb^2 + 3*EI))/((la - L + lb)*(Ac*L^2 - 2*Ac*L*la - 2*Ac*L*lb + Ac*la^2 + 2*Ac*la*lb + Ac*lb^2 + 12*EI)),     0, (6*Ac*EI*(L + la - lb))/((la - L + lb)*(Ac*L^2 - 2*Ac*L*la - 2*Ac*L*lb + Ac*la^2 + 2*Ac*la*lb + Ac*lb^2 + 12*EI)), -(2*EI*(Ac*L^2 + Ac*L*la + Ac*L*lb - 2*Ac*la^2 + 2*Ac*la*lb - 2*Ac*lb^2 - 6*EI))/((la - L + lb)*(Ac*L^2 - 2*Ac*L*la - 2*Ac*L*lb + Ac*la^2 + 2*Ac*la*lb + Ac*lb^2 + 12*EI))];
               [ -EA/L,                                                                                                                  0,                                                                                                                                                                          0,  EA/L,                                                                                                                 0,                                                                                                                                                                          0];
               [     0,               (12*Ac*EI)/((la - L + lb)*(Ac*L^2 - 2*Ac*L*la - 2*Ac*L*lb + Ac*la^2 + 2*Ac*la*lb + Ac*lb^2 + 12*EI)),                                                          (6*Ac*EI*(L + la - lb))/((la - L + lb)*(Ac*L^2 - 2*Ac*L*la - 2*Ac*L*lb + Ac*la^2 + 2*Ac*la*lb + Ac*lb^2 + 12*EI)),     0,             -(12*Ac*EI)/((la - L + lb)*(Ac*L^2 - 2*Ac*L*la - 2*Ac*L*lb + Ac*la^2 + 2*Ac*la*lb + Ac*lb^2 + 12*EI)),                                                          (6*Ac*EI*(L - la + lb))/((la - L + lb)*(Ac*L^2 - 2*Ac*L*la - 2*Ac*L*lb + Ac*la^2 + 2*Ac*la*lb + Ac*lb^2 + 12*EI))];
               [     0, -(6*Ac*EI*(L - la + lb))/((la - L + lb)*(Ac*L^2 - 2*Ac*L*la - 2*Ac*L*lb + Ac*la^2 + 2*Ac*la*lb + Ac*lb^2 + 12*EI)), -(2*EI*(Ac*L^2 + Ac*L*la + Ac*L*lb - 2*Ac*la^2 + 2*Ac*la*lb - 2*Ac*lb^2 - 6*EI))/((la - L + lb)*(Ac*L^2 - 2*Ac*L*la - 2*Ac*L*lb + Ac*la^2 + 2*Ac*la*lb + Ac*lb^2 + 12*EI)),     0, (6*Ac*EI*(L - la + lb))/((la - L + lb)*(Ac*L^2 - 2*Ac*L*la - 2*Ac*L*lb + Ac*la^2 + 2*Ac*la*lb + Ac*lb^2 + 12*EI)),     -(4*EI*(Ac*L^2 - 2*Ac*L*la + Ac*L*lb + Ac*la^2 - Ac*la*lb + Ac*lb^2 + 3*EI))/((la - L + lb)*(Ac*L^2 - 2*Ac*L*la - 2*Ac*L*lb + Ac*la^2 + 2*Ac*la*lb + Ac*lb^2 + 12*EI))]];
        % momentos y reaccciones de empotramiento
        EIa=10^10;
        Aca=10^10;
        %[Y1,M1,Y2,M2,~,~,~,~]=ecuacion_diferencial_viga(EIa,EI,EIc,Asa,Ac,Asc,q1,q2,L,la,lb,0,0,0,0,x);
        [Y1,M1,Y2,M2,~,~,~]=fe_nudos_rigidos(EIa,EI,Aca,Ac,q1,q2,L,la,lb,t1,t2,v1,v2,x);
        x=0;
        X1 = ((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L);
        x=L;
        X2 = ((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L);
        FE=[X1;Y1;M1;X2;Y2;M2];
        X1=abs(X1);
        Y1=abs(Y1);
        M1=-abs(M1);
        M2=-abs(M2);
   elseif ana==1 
       EIa=10^10;
       Aca=10^10;
       la=lalb(1);
       lb=lalb(2);
       % cortantes, momentos, deformadas verticales, deformadas axiales, fuerzas axiales 
       xlalb=x;
       %xlalb=x.*((x>la)&(x<=(L-lb)));
       f=find(x<=la);
       sizef=size(f,2);
       xlalb(f)=la*ones(1,sizef);
       f=find(x>=(L-lb));
       sizef=size(f,2);
       xlalb(f)=((L-lb)*ones(1,sizef));
       [~,~,~,~,V,M,v]=fe_nudos_rigidos(EIa,EI,Aca,Ac,q1,q2,L,la,lb,t1,t2,v1,v2,xlalb); 
       %[fax,V,M,u,v] = def_ec_dif(L, EI, Ac, EA, q1,q1, q1,q2,k,v1,v2,t1,t2,u1,u2);
       u =(- ((b1 - b2).*x.^3)/(3.*L) + (b1.*x.^2)/2 + EA.*u1)/EA + (x.*(((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L)))/EA;
       fax =((b1 - b2).*x.^2)/(2.*L) - b1.*x + (6.*EA.*u2 - 6.*EA.*u1 + 2.*L.^2.*b1 + L.^2.*b2)/(6.*L);  
   end
elseif strcmp(tipo_conti,'PD')==1
   if ana==0
        % matriz de rigidez local expresada en el sistema de coordenadas locales
        % analsis de  esbetez
        k=0;
        Kloc = Ke_ec_dif(L, EI, EA,Ac,k,P,puntos_graficas);
        % momentos y reaccciones de empotramiento
        [X1,Y1,M1,X2,Y2,M2] = fe_ec_dif(L, EI, Ac, EA, q1,q2,q1,q2,k,P,puntos_graficas);
        FE=[X1;Y1;M1;X2;Y2;M2];
        X1=abs(X1);
        Y1=abs(Y1);
        M1=-abs(M1);
        M2=-abs(M2);
   elseif ana==1 
       % cortantes, momentos, deformadas verticales, deformadas axiales, fuerzas axiales 
       [fax,V,M,u,v] = def_ec_dif(L, EI, Ac, EA, q1,q1, q1,q2,k,P,v1,v2,t1,t2,u1,u2,puntos_graficas);
        
   end 
end
if ana==0
    % matriz de mometos de empotramiento
    Me=[-X1 ,0  ,0,0  ,0  ,0;
         0  ,Y1,0,0   ,0  ,0;
         0  ,-M1,0,0  ,0  ,0;
         0  ,0  ,0,X2 ,0  ,0;
         0  ,0  ,0,0  ,-Y2,0;
         0  ,0  ,0,0  ,M2 ,0];
end