clear 
clc
syms m a b alfam Xe Ye D nu q
%simplemente apoyada Pagina 114

%w=(1/m^5*(1-(m*pi*b/(2*a)*tanh(alfam)+2)/(2*cosh(alfam))*cosh(2*alfam*Ye/b)+1/(2*cosh(alfam))*m*pi*Ye/a*sinh(2*alfam*Ye/b))*sin(m*pi*Xe/a));

%empotrado superior y inferior y apoyado en caras laterales  %pagina 185
w=1/m^5*sin(m*pi*Xe/a)*(1-(alfam*tanh(alfam)+2)/(2*cosh(alfam))...
    *cosh(m*pi*Ye/a)+1/(2*cosh(alfam))*m*pi*Ye/a*sinh(m*pi*Ye/a));

%empotrado en 4 caras 

%w=(-1)^((m-1)/2)/m^5*cos(m*pi*Xe/a)*(1-(alfam*tanh(alfam)+2)/(2*cosh(alfam))*cosh(m*pi*Ye/a)+1/(2*cosh(alfam))*m*pi*Ye/a*senh(m*pi*Ye/a));


Mx=simplify(-D*(diff(w,Xe,2)+nu*diff(w,Ye,2)));
My=simplify(-D*(diff(w,Ye,2)+nu*diff(w,Xe,2)));
Mxy=simplify(D*(1-nu)*diff(diff(w,Xe),Ye));
Myx=-simplify(D*(1-nu)*diff(diff(w,Xe),Ye));

Qx=-D*diff((diff(w,Xe,2)+diff(w,Ye,2)),Xe);
Qy=-D*diff((diff(w,Xe,2)+diff(w,Ye,2)),Ye);


