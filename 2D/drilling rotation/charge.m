function [FG]=charge (FG,charg,ddln)
nchex=size(charg,1);
for inoeud = 1:nchex
nn=charg(inoeud,1);    
for i=1:ddln

FG((nn-1)*ddln+i)= charg(inoeud,i+1);   
end    
end

end