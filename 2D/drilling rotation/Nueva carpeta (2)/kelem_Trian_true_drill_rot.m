function [kel] = kelem_Trian_true_drill_rot(xe) 

global nne E nu
kel=zeros(9,9);
nne=3;
..................................................
   
xlocal=xe(:,1);
ylocal=xe(:,2);
% compute constitutive matrix
[Dm]=mater;   

% sampling / quadrature / integration Points and weights
nphamer=3;
[xsi,eta,gw] = hamer(nphamer);

for k=1:nphamer
[detj,invj] = jacob3(xlocal,ylocal);
...........................................................
[bm] = bmatCTMTDR(xe,xsi(k),eta(k));  % Compatible Traingulat Membrane with True Drilling Rotation
...........................................................
kel=kel+bm'*Dm*bm*gw(k)*detj;
end
[DDD,AAA]=transformation(xe);
kel = transpose(DDD)*kel*DDD;
%kel = kel + 0.0001*detj/2*AAA'*E/(2+2*nu)*AAA;
end