
t=linspace(0,10,2.2)';
P=0.2*ones(100,1);%[0,5,8.6603,10,8.6603,5,0,0,0,0,0]'

dt=0.2;
m=1;
x=0.05;
w=3.1416;
c1=w^2;
c2=2*xi*w+w^2*dt;
c3=w*dt*(xi+w*dt/3);
c4=1+xi*w*dt+(w*dt)^2/6;
c5=dt/2;
c6=dt;
c7=dt^2/3;
c8=dt^2/6;
sizeP=size(P,1);
v=zeros(sizeP,1);
a=zeros(sizeP,1);
u=zeros(sizeP,1);

for i=1:(sizeP-1)
    a(i+1)=P(i+1)/m-c1*u(i)-c2*v(i)-c3*a(i);
    v(i+1)=v(i)+c5*(a(i)+a(i+1));
    u(i+1)=u(i)+c6*v(i)+c7*a(i)+c8*a(i+1);
end


