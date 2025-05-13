function [DDD,AAA]= transformation(xe)
% matrice de transformation  u=C.alpha
....... Polynomial xy coordinates ........
    
 for j=1:3
if j==1
i=3;
else
i=j-1;    
end
% si sx(ij)=xi-xj;
sx(i)=xe(i,1)-xe(j,1);
sy(i)=xe(i,2)-xe(j,2);
rl(i)=sqrt(sx(i)^2+sy(i)^2);
c(i)=-sy(i)/rl(i);
s(i)=sx(i)/rl(i);
% si sx(ij)=xj-xi;
%sx(j)=xe(j,1)-xe(i,1);
%sy(j)=xe(j,2)-xe(i,2);
%c(j)=sy(j)/rl(j);
%s(j)=-sx(j)/rl(j);
 end
 
%L'inverse de la matrice [C]
x2=xe(2,1);
x3=xe(3,1);
y3=xe(3,2);
DDD(1,:)=[1,0,0,0,0,0,0,0,0];
DDD(2,:)=[-1/x2,0,0,1/x2,0,0,0,0,0];
DDD(3,:)=[-(x2-x3)/(x2*y3),0,0,-(x3)/(x2*y3),0,0,1/y3,0,0];
DDD(4,:)=[0,1,0,0,0,0,0,0,0];
DDD(5,:)=[0,-1/x2,0,0,1/x2,0,0,0,0];
DDD(6,:)=[0,-(x2-x3)/(x2*y3),0,0,-(x3)/(x2*y3),0,0,1/y3,0];
DDD(7,:)=[0,0,-rl(1)/6,0,0,rl(1)/6,0,0,0];
DDD(8,:)=[0,0,0,0,0,-rl(2)/6,0,0,rl(2)/6];
DDD(9,:)=[0,0,rl(3)/6,0,0,0,0,0,-rl(3)/6];

.......................................................
% Spurious zero energy mode stabilisation
phiR=(-0.5*DDD(3,:)+0.5*DDD(5,:));
RRR=-phiR+[0,0,1/3,0,0,1/3,0,0,1/3];
DDD(7,:)=DDD(7,:)+0*rl(1)*RRR;
DDD(8,:)=DDD(8,:)+0*rl(2)*RRR;
DDD(9,:)=DDD(9,:)+0*rl(3)*RRR;

AAA=zeros(3,9);
SSSS=1;%rl(1)/3+rl(2)/3+rl(3)/3;
AAA(1,:)=1*rl(1)*RRR/SSSS;
AAA(2,:)=1*rl(2)*RRR/SSSS;
AAA(3,:)=1*rl(3)*RRR/SSSS;

.......................................................
% matrice pénalité
%x1=xe(1,1);
%y1=xe(1,2);
%y2=xe(2,2);

%air=abs(0.5*(sx(1)*sy(2)-sx(2)*sy(1)));
%AAA=[(x3-x2)/(4*air),-(y2-y3)/(4*air),1/3,(x1-x3)/(4*air),-(y3-y1)/(4*air),1/3,(x2-x1)/(4*air),-(y1-y2)/(4*air),1/3];
end
