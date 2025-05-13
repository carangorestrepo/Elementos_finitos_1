function [bm] = bmatCTMTDR(xe,xsi,eta)
    
%bm = zeros(3,9);
fnT3=zeros(1,3);
% les fonctions d'interpolation géométrique
fnT3(1)=1-xsi-eta;
fnT3(2)=xsi;
fnT3(3)=eta;
x = fnT3*xe(:,1); % x = N1 x1 + N2 x2 + N3 x3
y = fnT3*xe(:,2); % y = N1 y1 + N2 y2 + N3 y3
sx=zeros(1,3);
sy=zeros(1,3);
rl=zeros(1,3);
c=zeros(1,3);
s=zeros(1,3);
ij=zeros(3,2);
for j=1:3
    if j==1
        i=3;
    else
        i=j-1;    
    end
    ij(i,:)=[i,j];
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

 a2=4/xe(2,1);
 a3=-(4*xe(3,1))/(xe(2,1)*xe(3,2));
 a4=-4/xe(2,1)^2;
 a5=-(4*xe(2,1)-8*xe(3,1))/(xe(2,1)^2*xe(3,2));
 a6=-(4*xe(3,1)^2-4*xe(2,1)*xe(3,1))/(xe(2,1)^2*xe(3,2)^2);
 b5=4/(xe(2,1)*xe(3,2));
 b6=-(4*xe(3,1))/(xe(2,1)*xe(3,2)^2);
 d3=4/xe(3,2);
 d5=-4/(xe(2,1)*xe(3,2));
 d6=-(4*xe(2,1)-4*xe(3,1))/(xe(2,1)*xe(3,2)^2);
 
 abd=[a2,a3,a4,a5,a6,b5,b6,d3,d5,d6];
.........................................................................
bm =[0,1,0,0,0,0,c(1)*a2+2*c(1)*a4*x+c(1)*a5*y                              ,c(2)*b5*y                      ,c(3)*d5*y;
     0,0,0,0,0,1,s(1)*a3+s(1)*a5*x+2*s(1)*a6*y                              ,s(2)*b5*x+2*s(2)*b6*y          ,s(3)*d3+s(3)*d5*x+2*s(3)*d6*y;
     0,0,1,0,1,0,c(1)*a3+s(1)*a2+c(1)*a5*x+2*c(1)*a6*y+2*s(1)*a4*x+s(1)*a5*y,c(2)*b5*x+2*c(2)*b6*y+s(2)*b5*y,c(3)*d3+c(3)*d5*x+2*c(3)*d6*y+s(3)*d5*y];

end