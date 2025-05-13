function [depl] = local_disp(DG,connex,ro,iel)

global nne ddln

depg = zeros(nne*ddln,1);

for in=1:nne
 nn=connex(iel,in);
 ii=(nn-1)*ddln;
for i=1:ddln
 depg((in-1)*ddln+i,1)=DG(ii+i,1);
end
end
    
 depl=ro*depg;

end