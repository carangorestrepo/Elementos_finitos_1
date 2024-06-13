function [KG,FG]=cond_lim(KG,FG,climit,NTDDL,ddln)
for i=1:size(climit,1)
nn=climit(i,1);
for j=1:ddln
if climit(i,j+1)==123456
else
ii=(nn-1)*ddln+(j);
for jj=1:NTDDL
if jj==ii
%FG(ii)=climit(i,j+1);
FG(ii)=KG(ii,ii)*climit(i,j+1);
else
FG(jj)=FG(jj)-KG(ii,jj)*climit(i,j+1);
end
end
K=KG(ii,ii);
KG(ii,:)=zeros(1,NTDDL);
KG(:,ii)=zeros(NTDDL,1);
%KG(ii,ii)=1; 
KG(ii,ii)=K; 
end
end
end