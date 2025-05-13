function KG = assemblage(KG,connex,keg,nne,ddln,iel)

NN=connex(iel,:);
nddle=nne*ddln;
ddle=zeros(1,nddle);
k=0;
for i=1:nne
start=(NN(i)-1)*ddln;
for j=1:ddln
k=k+1;
ddle(k)=start+j;
end
end

for i=1:nddle
ii=ddle(i);
for j=1:nddle
jj=ddle(j);
KG(ii,jj)=KG(ii,jj)+keg(i,j);
end
end

end