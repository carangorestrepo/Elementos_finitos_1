function [rog,xe] = rotaT3(coor,connex,iel)
% *************  matrice rotation  **************
%  membrane triangulair aved ddl rotation nodale
nne = 3;

x=zeros(nne,1);
y=zeros(nne,1);
for i=1:nne 
nn(i)=connex(iel,i);           % extract connected node for (iel element)
x(i)=coor(nn(i),1);            % extract x coordinate
y(i)=coor(nn(i),2);            % extract y coordinate
end

 xi(1)=x(1);
 xj(1)=x(2);
 xk(1)=x(3);
 xi(2)=y(1);
 xj(2)=y(2);
 xk(2)=y(3);
 
 xji=xj(1)-xi(1);
 yji=xj(2)-xi(2);
 xki=xk(1)-xi(1);
 yki=xk(2)-xi(2);
 cij=sqrt(xji^2+yji^2);
 ro(1,1)=xji/cij;
 ro(1,2)=yji/cij;
 c=xji*yki-xki*yji;
 ro(2,1)=-ro(1,2);
 ro(2,2)=ro(1,1);
 ch=ro(2,1)*xki+ro(2,2)*yki;
 if ch < 0
 ch=-ch;
 for j=1:2
 ro(2,j)=-ro(2,j);
 end
 end
 chi=-ro(1,1)*xki-ro(1,2)*yki;
 chj=cij+chi;
 xe(1,1)=chi;
 xe(1,2)=0;
 xe(2,1)=chj;
 xe(2,2)=0;
 xe(3,1)=0;
 xe(3,2)=ch;
 
%  coor locales / axe x-x ---> 1-2 avec (0,0) sur le  1
 xe(2,1)=xe(2,1)-xe(1,1);
 xe(3,1)=xe(3,1)-xe(1,1);
 xe(1,1)=0;
 
 rog = zeros(nne*3,nne*3);
 for ij=1:3:nne*3
 ji=ij-1;
 for k=1:2
 for l=1:2
 rog(ji+k,ji+l)=ro(k,l);
 end
 end
 rog(ji+3,ji+3)=1;
 end     
 
end