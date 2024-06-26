function  Sa=espectro(Aa,Av,I,suelo,t)



Fat=[0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5;0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8;1 1 1 1 1 1 1 1 1;1.2 1.2 1.2 1.15 1.1 1.05 1 1 1;1.6 1.5 1.4 1.3 1.2 1.15 1.1 1.05 1;2.5 2.1 1.7 1.45 1.2 1.05 0.9 0.9 0.9];

Fvt=[0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5;0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8;1 1 1 1 1 1 1 1 1;1.7 1.65 1.6 1.55 1.5 1.45 1.4 1.35 1.3;2.4 2.2 2 1.9 1.8 1.7 1.6 1.55 1.5;3.5 3.35 3.2 3 2.8 2.6 2.4 2.4 2.4];

fa=Fat(1,:)==Aa;
fv=Fat(1,:)==Av;

Fa=Fat(suelo,fa);
Fv=Fvt(suelo,fv);
Ir=[1,1.1,1.25,1.5];
I=Ir(I);

Sa=2.5*Aa*Fa*I.*(t<=0.48*Av*Fv/(Aa*Fa))+1.2*Av*Fv*I./t.*((t>0.48*Av*Fv/(Aa*Fa))&(t<=2.4*Fv))+1.2*Av*Fv*2.4*Fv*I./t.^2.*((t>2.4*Fv)&(t<=2.4*Fv*2));