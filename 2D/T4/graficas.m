figure
contourf(Xe,Ye,sxx,12)
ylabel('\sigma_x (Pa)','FontSize',10); axis equal tight; colorbar;

figure
contourf(Xe,Ye,syy,12)
ylabel('\sigma_y (Pa)','FontSize',10); axis equal tight; colorbar;

figure
contourf(Xe,Ye,txyy,12)
ylabel('\tau_{xy} (Pa)','FontSize',10); axis equal tight; colorbar;



syms xo xi yo yi
xy=[xo,0;
    xo,yo;
    xi,yi;
    xi,0;
    xo,0]
[AREA,XCEN,YCEN]=propiedadesgeo(xy)
area=zeros(Nx-1,1);
xcen=zeros(Nx-1,1);
%M=zeros(Ny,1);
for i=1:(Nx-1)
    yi=sxx(1,i+1);
    yo=sxx(1,i);
    xi=Xe(1,i+1);
    xo=Xe(1,i);

    xcen(i)=((xo^2*yo)/2 - (xi^2*yi)/2 + (yi/8 - yo/8)*((xi - xo)^2/3 + (xi + xo)^2))/(xo*yo - xi*yi + ((xi + xo)*(yi - yo))/2);
    area(i)=xi*yi - xo*yo - ((xi + xo)*(yi - yo))/2;
  
end
M=sum(xcen.*area)



