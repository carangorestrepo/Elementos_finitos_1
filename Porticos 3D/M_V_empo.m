function [Ra,Rb,Va,Vb,Ma,Mb]=M_V_empo(q1,q2,Ac,EI,L,EA,b1,b2)  

u1=0;
u2=0;
v1=0;
v2=0;
t1=0;
t2=0;

x=0;
Va =-(q1*x - (240*Ac*EI*v2 - 240*Ac*EI*v1 + 7*Ac*L^4*q1 + 3*Ac*L^4*q2 + 80*EI*L^2*q1 + 40*EI*L^2*q2 - 120*Ac*EI*L*t1 - 120*Ac*EI*L*t2)./(20*L*(Ac*L^2 + 12*EI)) - (x.^2*(q1 - q2))/(2*L));
Ma =(720*EI^2*t2 - 720*EI^2*t1 + 3*Ac*L^5*q1 + 2*Ac*L^5*q2 + 30*EI*L^3*q1 + 30*EI*L^3*q2 - 240*Ac*EI*L^2*t1 - 120*Ac*EI*L^2*t2 - 360*Ac*EI*L*v1 + 360*Ac*EI*L*v2)./(60*L*(Ac*L^2 + 12*EI)) - (q1*x.^2)/2 - x.*(((q1 - q2)*x.^2)/(2*L) - q1*x + (240*Ac*EI*v2 - 240*Ac*EI*v1 + 7*Ac*L^4*q1 + 3*Ac*L^4*q2 + 80*EI*L^2*q1 + 40*EI*L^2*q2 - 120*Ac*EI*L*t1 - 120*Ac*EI*L*t2)./(20*L*(Ac*L^2 + 12*EI))) + (x.^3*(q1 - q2))./(3*L); 
Ra =((b1 - b2)*x.^2)./(2*L) - b1*x + (6*EA*u2 - 6*EA*u1 + 2*L^2*b1 + L^2*b2)/(6*L);

x=L;
Vb =q1*x - (240*Ac*EI*v2 - 240*Ac*EI*v1 + 7*Ac*L^4*q1 + 3*Ac*L^4*q2 + 80*EI*L^2*q1 + 40*EI*L^2*q2 - 120*Ac*EI*L*t1 - 120*Ac*EI*L*t2)./(20*L*(Ac*L^2 + 12*EI)) - (x.^2*(q1 - q2))/(2*L);
Mb =(720*EI^2*t2 - 720*EI^2*t1 + 3*Ac*L^5*q1 + 2*Ac*L^5*q2 + 30*EI*L^3*q1 + 30*EI*L^3*q2 - 240*Ac*EI*L^2*t1 - 120*Ac*EI*L^2*t2 - 360*Ac*EI*L*v1 + 360*Ac*EI*L*v2)./(60*L*(Ac*L^2 + 12*EI)) - (q1*x.^2)/2 - x.*(((q1 - q2)*x.^2)/(2*L) - q1*x + (240*Ac*EI*v2 - 240*Ac*EI*v1 + 7*Ac*L^4*q1 + 3*Ac*L^4*q2 + 80*EI*L^2*q1 + 40*EI*L^2*q2 - 120*Ac*EI*L*t1 - 120*Ac*EI*L*t2)./(20*L*(Ac*L^2 + 12*EI))) + (x.^3*(q1 - q2))./(3*L);
Rb =-(((b1 - b2)*x.^2)./(2*L) - b1*x + (6*EA*u2 - 6*EA*u1 + 2*L^2*b1 + L^2*b2)/(6*L));

