function cargas(carga1,carga2,xi,xf,yi,yf,h1,h2,dir)

%carga=20;
%xi=0;
%xf=5;
%yi=0;
%yf=3;
%h1=0.3;
%h2=1;
%dir=19;

if dir==2%% dirección vertical
    ha = annotation('arrow');
    ha.Parent = gca;           % associate the arrow the the current axes
    ha.X = [xi, xi];          %
    ha.Y = [yi-h1, yi];

    ha = annotation('arrow');
    ha.Parent = gca;           % associate the arrow the the current axes
    ha.X = [xf, xf];          %
    ha.Y = [yf-h2, yf];

    ha = annotation('arrow');
    ha.Parent = gca;           % associate the arrow the the current axes
    ha.X = [(xf+xi)/2, (xf+xi)/2];          %
    ha.Y = [(yf-h2+yi-h1)/2, (yf+yi)/2];
    
    
    xy=[xi,yi-h1;xf,yf-h2;xf,yf;xi,yi;xi,yi-h1];
    line(xy(:,1),xy(:,2))
    
    fill(xy(:,1),xy(:,2),'b')
    alpha(0.2)
    
    text(xi,yi-h1, num2str(carga1),'FontSize',14,'Color','m');
    text(xf,yf-h2, num2str(carga2),'FontSize',14,'Color','m');
    
elseif dir==29%% dirección vertical proyectada  
    
    ha = annotation('arrow');
    ha.Parent = gca;           % associate the arrow the the current axes
    yn=max([yi,yf]);
    
    ha.X = [xi, xi];          %
    ha.Y = [yn+h1, yn];

    ha = annotation('arrow');
    ha.Parent = gca;           % associate the arrow the the current axes
    ha.X = [xf, xf];          %
    ha.Y = [yn+h2, yn];
    
    ha = annotation('arrow');
    ha.Parent = gca;           % associate the arrow the the current axes
    ha.X = [(xf+xi)/2, (xf+xi)/2];          %
    ha.Y = [(yn+h2+yn+h1)/2, (yn+yn)/2];

    xy=[xi,yn+h1;xf,yn+h2;xf,yn;xi,yn;xi,yn+h1];
    line(xy(:,1),xy(:,2))
    fill(xy(:,1),xy(:,2),'b')
    alpha(0.2)
    
    text(xi,yn-h1, num2str(carga1),'FontSize',14,'Color','m');
    text(xf,yn-h2, num2str(carga2),'FontSize',14,'Color','m');
    
elseif dir==1 %% dirección horizontal
    ha = annotation('arrow');
    ha.Parent = gca;           % associate the arrow the the current axes
    ha.X = [xi-h1, xi];          %
    ha.Y = [yi, yi];

    ha = annotation('arrow');
    ha.Parent = gca;           % associate the arrow the the current axes
    ha.X = [xf-h2, xf];          %
    ha.Y = [yf, yf];

    ha = annotation('arrow');
    ha.Parent = gca;           % associate the arrow the the current axes
    ha.X = [(xf-h2+xi-h1)/2,(xf+xi)/2];          %
    ha.Y = [(yf+yi)/2, (yf+yi)/2];
    
    xy=[xi-h1,yi;xf-h2,yf;xf,yf;xi,yi;xi-h1,yi];
    line(xy(:,1),xy(:,2))
    fill(xy(:,1),xy(:,2),'b')
    alpha(0.2)
    
    text(xi-h1,yi, num2str(carga1),'FontSize',14,'Color','m');
    text(xf-h2,yf, num2str(carga2),'FontSize',14,'Color','m');
    
elseif dir==19%% dirección vertical proyectada
    ha = annotation('arrow');
    ha.Parent = gca;           % associate the arrow the the current axes
    xn=min([xi,xf]);
    ha.X = [xn-h1, xn];
    ha.Y = [yi, yi];

    ha = annotation('arrow');
    ha.Parent = gca;           % associate the arrow the the current axes
    ha.X = [xn-h2, xn];          %
    ha.Y = [yf, yf];
    
    ha = annotation('arrow');
    ha.Parent = gca;           % associate the arrow the the current axes
    ha.X = [(xn-h1+xn-h2)/2, (xn+xn)/2];          %
    ha.Y = [(yf+yi)/2,(yf+yi)/2];

    xy=[xn-h1,yi;xn-h2,yf;xn,yf;xn,yi;xn-h1,yi];
    line(xy(:,1),xy(:,2))   
    fill(xy(:,1),xy(:,2),'b')
    alpha(0.2)
    
    text(xn-h1,yi, num2str(carga1),'FontSize',14,'Color','m');
    text(xn-h2,yf, num2str(carga2),'FontSize',14,'Color','m')
    
end