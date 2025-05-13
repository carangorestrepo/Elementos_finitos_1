function [detj,InvJ] = jacob3(xlocal,ylocal)
%.......................................................
%*****************    jacobian matrix     ******************
nne = 3;

     jacob=zeros(2,2);

     [dndxsi,dndeta]=T3_deriv;
     
     for i=1: nne
     jacob(1,1)=jacob(1,1)+dndxsi(i)*xlocal(i);
     jacob(1,2)=jacob(1,2)+dndxsi(i)*ylocal(i);
     jacob(2,1)=jacob(2,1)+dndeta(i)*xlocal(i);
     jacob(2,2)=jacob(2,2)+dndeta(i)*ylocal(i);
     end
 
     
      detj=jacob(1,1)*jacob(2,2)-jacob(2,1)*jacob(1,2);

      InvJ(1,1)=jacob(2,2)/detj;
      InvJ(1,2)=-jacob(1,2)/detj;
      InvJ(2,1)=-jacob(2,1)/detj;
      InvJ(2,2)=jacob(1,1)/detj;
      
end

function [dndxsi,dndeta] = T3_deriv

    dndxsi(1)=-1;
	dndxsi(2)=1;
	dndxsi(3)=0;

	dndeta(1)=-1;
	dndeta(2)=0;
	dndeta(3)=1;
end