function [] = plot_stress_T3(coor,connex,NTE,effort)
% **********  contour des contraintes généralisées   ***********

%plot des contraintes    
    
[effort_nod]=effort_nod_moy_T3(coor,connex,effort,NTE);
figure (2)
trisurf(connex,coor(:,1),coor(:,2),effort_nod);
view(0,90);
axis equal;
%axis ([0 1 0 1])
shading interp   % flat, faceted, interp
colormap;
colorbar('location','southoutside');
title('contraintes généralisées')
grid off;

end
