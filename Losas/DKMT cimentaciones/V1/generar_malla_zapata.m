function [xnod,LaG,in_col,info] = generar_malla_zapata(geo,malla)
% Malla triangular estructurada y conforme con la huella del pedestal.

validar_geometria(geo);

x1 = geo.cx0;
x2 = geo.cx0+geo.lcx;
y1 = geo.cy0;
y2 = geo.cy0+geo.lcy;

xvec = unir_segmentos( ...
    segmento(0,x1,malla.h_fuera), ...
    segmento(x1,x2,malla.h_huella), ...
    segmento(x2,geo.Lx,malla.h_fuera));

yvec = unir_segmentos( ...
    segmento(0,y1,malla.h_fuera), ...
    segmento(y1,y2,malla.h_huella), ...
    segmento(y2,geo.Ly,malla.h_fuera));

[X,Y] = meshgrid(xvec,yvec);
xnod = [X(:),Y(:)];

nx = length(xvec);
ny = length(yvec);

LaG = zeros(2*(nx-1)*(ny-1),3);
e = 0;

for j = 1:ny-1
    for i = 1:nx-1
        n1 = sub2ind([ny,nx],j,i);
        n2 = sub2ind([ny,nx],j,i+1);
        n3 = sub2ind([ny,nx],j+1,i+1);
        n4 = sub2ind([ny,nx],j+1,i);

        if isfield(malla,'diagonal_alternada') && ...
                malla.diagonal_alternada && mod(i+j,2)==0
            e=e+1; LaG(e,:)=[n1,n2,n4];
            e=e+1; LaG(e,:)=[n2,n3,n4];
        else
            e=e+1; LaG(e,:)=[n1,n2,n3];
            e=e+1; LaG(e,:)=[n1,n3,n4];
        end
    end
end

for e = 1:size(LaG,1)
    nod = LaG(e,:);
    xe = xnod(nod,1);
    ye = xnod(nod,2);
    detJ = (xe(2)-xe(1))*(ye(3)-ye(1)) ...
         - (xe(3)-xe(1))*(ye(2)-ye(1));
    if detJ < 0
        LaG(e,[2 3])=LaG(e,[3 2]);
    elseif detJ == 0
        error('Se genero un elemento de area nula.');
    end
end

cent = (xnod(LaG(:,1),:)+xnod(LaG(:,2),:)+xnod(LaG(:,3),:))/3;
tol = 1e-10;

in_col = cent(:,1) >= x1-tol & cent(:,1) <= x2+tol & ...
         cent(:,2) >= y1-tol & cent(:,2) <= y2+tol;

Ae = areas_elementos(xnod,LaG);

info.xvec = xvec;
info.yvec = yvec;
info.area_elemento = Ae;
info.area_total = sum(Ae);
info.area_huella = sum(Ae(in_col));

end

function validar_geometria(geo)
if geo.Lx<=0 || geo.Ly<=0 || geo.hz<=0
    error('Las dimensiones de la zapata deben ser positivas.');
end
if geo.lcx<=0 || geo.lcy<=0
    error('Las dimensiones de la huella deben ser positivas.');
end
if geo.cx0<0 || geo.cy0<0 || ...
   geo.cx0+geo.lcx>geo.Lx || geo.cy0+geo.lcy>geo.Ly
    error('La huella debe estar contenida en la zapata.');
end
if geo.pd<geo.hz
    error('pd debe ser mayor o igual que hz.');
end
end

function v = segmento(a,b,h)
if b<a
    error('Intervalo geometrico invalido.');
end
if abs(b-a)<eps
    v=a;
    return
end
n=max(1,ceil((b-a)/h));
v=linspace(a,b,n+1);
end

function v = unir_segmentos(varargin)
v=[];
for k=1:nargin
    vk=varargin{k};
    if isempty(v)
        v=vk;
    else
        v=[v,vk(2:end)]; %#ok<AGROW>
    end
end
v=unique(v,'stable');
end

function Ae = areas_elementos(xnod,LaG)
xe1=xnod(LaG(:,1),1); ye1=xnod(LaG(:,1),2);
xe2=xnod(LaG(:,2),1); ye2=xnod(LaG(:,2),2);
xe3=xnod(LaG(:,3),1); ye3=xnod(LaG(:,3),2);

Ae=0.5*abs((xe2-xe1).*(ye3-ye1)-(xe3-xe1).*(ye2-ye1));
end
