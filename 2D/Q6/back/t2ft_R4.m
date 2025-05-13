function fte=t2ft_R4(xnod, lado, carga, espesor)

    %'''Función que convierte las fuerzas superficiales aplicadas a un elemento
    %finito rectangular de 4 nodos a sus correspondientes cargas nodales 
    %equivalentes ft
    
    %Recibe:
    %    xnod:  coordenadas nodales del elemento finito de 4 nodos
    %      xnod = [ x1e y1e
    %               x2e y2e
    %               x3e y3e
    %               x4e y4e ]

    %   lado:  arista en la que se aplica la carga, puede tomar los siguientes
    %           valores: 12, 23, 34, 41

    %    carga: fuerza distribuida en los nodos
        
    %           [ t1x t1y t2x t2y ]; % si carga se aplica sobre lado 12
    %           [ t2x t2y t3x t3y ]; % si carga se aplica sobre lado 23
    %           [ t3x t3y t4x t4y ]; % si carga se aplica sobre lado 34you must be very proud of her

    %           [ t4x t4y t1x t1y ]; % si carga se aplica sobre lado 41
    
    %     espesor: espesor del elemento
    %'''
%% Se definen algunas constantes
X = 1; Y = 2;    
 %se definen los indices de los lados
if   lado == 12
    idx=[1,2];
    %odenen=[1, 2, 3, 4];
elseif lado == 23
    idx=[2,3];
    %odenen=[3, 4, 5, 6];
elseif lado == 34
    idx=[3,4];
    %odenen=[5, 6, 7, 8];
elseif lado == 41
    idx=[4,1];
    %odenen=[7, 8, 1, 2];
else 
    disp('Únicamente se permiten los lados 12, 23, 34 o 41')
end
nno = size(xnod,1);
if nno ~= 4    
    disp( 'Solo para elementos rectangulares de 4 nodos');
end


%% Cuadratura de Gauss-Legendre
%NOTA: se asumirá aquí el mismo orden de la cuadratura tanto en la dirección
%       de xi como en la dirección de eta
n_gl         = 5; % orden de la cuadratura de Gauss-Legendre
[x_gl, w_gl] = gausslegendre_quad(n_gl);

% se definen las funciones de forma unidimensionales y sus derivadas
NN      = @(xi)[ (1-xi)/2,...
                (1+xi)/2 ];
dNN_dxi = @(xi)[     -1/2,...    
                     1/2 ];

% se calcula el vector de fuerzas distribuidas en los nodos
te = zeros(1,2*nno);
te(reshape([2*idx-1; 2*idx],4,1)) = carga(:);

% cálculo de la integral:
suma   = zeros(2*nno, 2*nno);
N      = zeros(nno,1);
dN_dxi = zeros(nno,1);

for p = 1:n_gl
   N(idx)      = NN(x_gl(p));
   
   for i =1:nno
        Nijk(:,[2*i-1, 2*i]) = [[N(i), 0   ];
                               [0,    N(i)]];
   end
 
   dN_dxi(idx) = dNN_dxi(x_gl(p));
   
   dx_dxi = dN_dxi'*xnod(:,X);
   dy_dxi = dN_dxi'*xnod(:,Y);
   ds_dxi = hypot(dx_dxi, dy_dxi);

   suma = suma + Nijk'*Nijk*ds_dxi*w_gl(p);
end

fte = espesor*suma*te';


