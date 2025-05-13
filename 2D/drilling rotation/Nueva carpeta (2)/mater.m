function [dm] = mater
%*******************    materiaux    *******************
 .....................  isotrope  ......................

global E nu Gxy th

   E; nuyx = nu; nuxy = nu;
   Gxy=E/2/(1+nuxy);
   
    dm = zeros(3,3);

	a = 0;
	if a == 0
%   contrainte plane
      dm(1,1)=1;
      dm(1,2)=nuxy;
      dm(2,1)=nuyx;
      dm(2,2)=1;
      dm(3,3)=(1-nuxy)/2;
      for i=1:3
      for j=1:3
      dm(i,j)=dm(i,j)*E*th/(1-nuyx*nuxy);
      end
      end
    end
    
    if a == 1
%   deformation plane
 	  un = 1-nuxy;
	  d = E*th/((1+nuxy)*(1-2*nuxy));
	  dm(1,1)=un*d;
      dm(1,2)=nuxy*d;
      dm(2,1)=nuyx*d;
      dm(2,2)=un*d;
      dm(3,3)=0.5*d*(1-2*nuxy);
    end
      
end
