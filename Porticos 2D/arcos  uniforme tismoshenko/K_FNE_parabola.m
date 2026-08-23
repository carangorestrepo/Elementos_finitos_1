function [K,R]=K_FNE_parabola(x1,x2,x3,y1,y2,y3,E,G,A,I,alpha,L,h,qb,qa)
dato = 100;

a=-(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)/((x1 - x2)*(x1*x2 - x1*x3 - x2*x3 + x3^2));
b=(x1^2*y2 - x2^2*y1 - x1^2*y3 + x3^2*y1 + x2^2*y3 - x3^2*y2)/((x1 - x2)*(x1*x2 - x1*x3 - x2*x3 + x3^2));
c=-(- y3*x1^2*x2 + y2*x1^2*x3 + y3*x1*x2^2 - y2*x1*x3^2 - y1*x2^2*x3 + y1*x2*x3^2)/((x1 - x2)*(x1*x2 - x1*x3 - x2*x3 + x3^2));

%a = p(1); b = p(2); c = p(3);
% Geometría del arco en cada punto x
x = linspace(x1, x3, dato)-x1;
y = a*(x+x1).^2 + b*(x+x1) + c-y1;                 % coordenada y
tan_theta = 2*a*x + b;               % pendiente (dy/dx)

[RR,R]=parabola_4(x,y,tan_theta,E,G,A,I,alpha,dato,L,h,qa,qb);

x = linspace(-x3, x1, dato)+x1;
y = a*(x+x3).^2 + b*(x+x3) + c-y3;                % coordenada y
tan_theta = b + a*(2*x + 2*x3);               % pendiente (dy/dx)
L=-x3+x1;
h=y3-y1;
[RR1,~]=parabola_4(x,y,tan_theta,E,G,A,I,alpha,dato,L,h,qb,qa);

K=[RR,RR1];