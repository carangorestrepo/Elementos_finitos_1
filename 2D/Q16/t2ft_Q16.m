function ft = t2ft_Q16(xnod, LaG, nodos_ijk, carga, espesor)
% Esta función convierte las fuerzas superficiales aplicadas a un elemento
% finito rectangular de 8 (serendipito) nodos a sus correspondientes cargas 
% nodales equivalentes ft
%
% SERENDIPITO 8          
% xnod = [ x1e y1e       
%          x2e y2e       
%          ... ...       
%          x8e y8e ];
%
% lado = 123, 345, 567, 781
%
% carga = [ t1x t1y t2x t2y t3x t3y ]; % si carga se aplica sobre lado 123
%         [ t3x t3y t4x t4y t5x t5y ]; % si carga se aplica sobre lado 345
%         [ t5x t5y t6x t6y t7x t7y ]; % si carga se aplica sobre lado 567
%         [ t7x t7y t8x t8y t1x t1y ]; % si carga se aplica sobre lado 781

%% Se definen algunas constantes
X = 1; Y = 2;

%% Parametros de la cuadratura de Gauss-Legendre
n_gl = 3;                        % orden de la cuadratura de Gauss-Legendre
[x_gl, w_gl] = gausslegendre_quad(n_gl);

%% Se definen las funciones de forma unidimensionales y sus derivadas

NN = @(xi) [- (9*xi^3)/16 + (9*xi^2)/16 + xi/16 - 1/16;
           (27*xi^3)/16 - (9*xi^2)/16 - (27*xi)/16 + 9/16;
           (27*xi)/16 - (9*xi^2)/16 - (27*xi^3)/16 + 9/16;
           (9*xi^3)/16 + (9*xi^2)/16 - xi/16 - 1/16];

dNN_dxi = @(xi)[(9*xi)/8 - (27*xi^2)/16 + 1/16;
                (81*xi^2)/16 - (9*xi)/8 - 27/16;
                27/16 - (81*xi^2)/16 - (9*xi)/8;
                (27*xi^2)/16 + (9*xi)/8 - 1/16];

%% Se definen las coordenadas locales de los lados
if     isequal(nodos_ijk, LaG([1 2 3 4]))
    idx = [ 1 2 3 4 ];
elseif isequal(nodos_ijk, LaG([4 5 6 7]))
    idx = [ 4 5 6 7 ];
elseif isequal(nodos_ijk, LaG([7 8 9 10]))
    idx = [ 7 8 9 10];
elseif isequal(nodos_ijk, LaG([10 11 12 1]))
    idx = [ 10 11 12 1 ];
else 
    error('Lado para la carga incorrectamente especificado.')
end

%% Se calcula el vector de fuerzas distribuidas en los nodos
te = zeros(16*2,1);
te(reshape([2*idx-1; 2*idx],8,1)) = carga(:);%% 8 grados de libertad ppor lado de elemento finito

%% Se calcula la integral
suma   = zeros(16*2);
N      = zeros(1,16);
dN_dxi = zeros(1,16);
for p = 1:n_gl
   N(idx)      = NN(x_gl(p));
   
   Nijk = [ N(1) 0    N(2) 0    N(3) 0    N(4) 0    N(5) 0    N(6) 0    N(7) 0    N(8) 0    N(9) 0    N(10) 0    N(11) 0    N(12) 0     N(13) 0     N(14) 0     N(15) 0     N(16) 0; 
            0    N(1) 0    N(2) 0    N(3) 0    N(4) 0    N(5) 0    N(6) 0    N(7) 0    N(8)   0  N(9) 0     N(10) 0    N(11) 0    N(12) 0     N(13) 0     N(14) 0     N(15) 0     N(16)];
 
   dN_dxi(idx) = dNN_dxi(x_gl(p));
   
   dx_dxi = dN_dxi*xnod(:,X);
   dy_dxi = dN_dxi*xnod(:,Y);
   ds_dxi = hypot(dx_dxi, dy_dxi);

   suma = suma + Nijk'*Nijk*ds_dxi*w_gl(p);
end

ft = espesor*suma*te;

%%
return;







