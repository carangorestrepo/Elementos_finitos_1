function [effort_nod]=effort_nod_moy_T3(coor,connex,efforts,NTE)

%------------------------------------------------------------------------
%  Calcul des contraintes (g�n�ralis�es) moyennes aux noeuds commains
%                       aux plusieurs �l�ments
%  efforts = matrices des efforts nodaux pour chaque �l�ment
%  nadjele_nod = number d'�l�ments adjacent autour de chaue noeud
%------------------------------------------------------------------------

global nne

NTN = size(coor,1);        
effort_nod=zeros(1,NTN);
nadjele_nod=zeros(NTN,1);
kk=1;

for iel=1:NTE   % loop for the total number of element
    
    for ii=1:nne
    effort_elem(ii) = efforts(kk);
    kk=kk+1;
    end
    
    for i=1:nne
    nn(i)=connex(iel,i);   % extract connected node for (iel element)
    end          
        
    for i=1:nne  % loop for 4 nodes of (iel)-th element              
        
        effort_nod(nn(i)) = effort_nod(nn(i)) + effort_elem(i);
        nadjele_nod(nn(i)) = nadjele_nod(nn(i)) + 1;              
    end
end

%------------------------------------------------------------------------
% contrainte g�n�ralis�e moyenne aux noeuds � partir d'�l�ments adjacent
%------------------------------------------------------------------------
for i=1:NTN
    effort_nod(i)=effort_nod(i)/nadjele_nod(i);
end

end