q1=20;
q2=30;
d=0.3;
bf=0.15;
tf=0.009;
tw=0.006;

A=b*h;
I=bf*d^3/12-(bf-tw)*(d-tf*2)^3/12;
L=2.5;
b1=15;
b2=30;
E=200000*1000;
v=0.3;
G=E/(2*(1+v));
EI=E*I;
EA=E*A;

Ac=tw*d*G;
t1=0;
t2=0;
v1=0;
v2=0;
u1=0;
u2=0;
x=linspace(0,L,200);



%[V,M,v,u,fax,Va,Vb,Ma,Mb,Ra,Rb]=cor_mome_def(q1,q2,Ac,EI,L,t1,t2,v1,v2,EA,b1,b2,u1,u2,x);
[Ra,Rb,Va,Vb,Ma,Mb]=M_V_empo(q1,q2,Ac,EI,L,EA,b1,b2)  