function [forcm]= efforts_elem(xe,depl)

global nne

..................................................
   
xlocal=xe(:,1);   % 
ylocal=xe(:,2);   % 
% compute constitutive matrix
[Dm]=mater;   

% les points de calcul des contraintes généralisées
% aux points d'integration
 xsi(1)= 0.5;
 xsi(2)= 0.5;
 xsi(3)= 0;      
 eta(1)= 0;
 eta(2)= 0.5;
 eta(3)= 0.5;

[DDD,AAA]=transformation(xe);
for k=1:3         
[detj,invj] = jacob3(xlocal,ylocal);
.....................................................
[bm] = bmatCTMTDR(xe,xsi(k),eta(k));  % Compatible Traingulat Membrane with True Drilling Rotation
fm = Dm*((bm*DDD)*depl);
....................................................

forcm(k,:)=fm';     % aux points d'integration

....................................................
% extrapolation des contraintes aux noeuds (interpolation linéaire)
% les fonctions d'interpolation géométrique
fnT3(1) = 1-xsi(k)-eta(k);
fnT3(2) = xsi(k);
fnT3(3) = eta(k);
xi(k) = fnT3*xe(:,1); % x = N1 x1 + N2 x2 + N3 x3
yi(k) = fnT3*xe(:,2); % y = N1 y1 + N2 y2 + N3 y3
end

X= [1 xi(1) yi(1);
    1 xi(2) yi(2);
    1 xi(3) yi(3)];

XX=[1 xe(1,1) xe(1,2);
    1 xe(2,1) xe(2,2);
    1 xe(3,1) xe(3,2)];

for j=1:3
forcm(:,j)=XX/X*forcm(:,j);     % aux noeuds
end

end