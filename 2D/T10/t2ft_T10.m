function   fte=t2ft_T10(xnod, lado, carga, espesor)
    %'''Función que convierte las fuerzas superficiales aplicadas a un elemento
    %finito triangular de 10 nodos a sus correspondientes cargas nodales equiva-
    %lentes ft

    %Recibe:
    %    xnod:  coordenadas nodales del elemento finito

    %    xnod = [ x1e y1e
    %            x2e  y2e
    %             ... ...
    %            x8e  y8e ]

    %    lado:  arista en la que se aplica la carga, puede tomar los siguientes
    %           valores: 1452, 2673, 3891

    %    carga: fuerza distribuida en los nodos
        
    
    %    espesor: espesor del elemento
    %'''
    
    %# se definen los indices de los lados
    if   lado == 1452 
        idx = [ 1, 4, 5, 2 ];
    elseif lado == 2673
         idx = [ 2, 6, 7, 3 ]; 
    elseif lado == 3891 
        idx = [ 3, 8, 9, 1 ]; 
    else 
        raise Exception('Únicamente se permiten los lados 1452, 2673 o 3891')
    end

    nno = size(xnod,1);
    if nno ~= 10    
        disp( 'Solo para elementos rectangulares de 4 nodos');
    end
    %# parámetros para mejorar la lectura del código
    X=1;
    Y=2;
   
    %# se define el número de puntos de la cuadratura y se obtienen los puntos
    %# de evaluación y los pesos de la cuadratura

    n_gl       = 4;
    %xi_gl, w_gl = leggauss(n_gl)
    
    xi_gl=[-0.86113631 -0.33998104  0.33998104  0.86113631];
    w_gl=[0.34785485 0.65214515 0.65214515 0.34785485];
    %# se definen las funciones de forma unidimensionales y sus derivadas
    NN     = @(xi)[-9*xi^3/16 + 9*xi^2/16 + xi/16 - 1/16
                                   27*xi^3/16 - 9*xi^2/16 - 27*xi/16 + 9/16
                                  -27*xi^3/16 - 9*xi^2/16 + 27*xi/16 + 9/16
                                    9*xi^3/16 + 9*xi^2/16 - xi/16 - 1/16];

    dNN_dxi = @(xi)[-27*xi^2/16 + 9*xi/8 + 1/16
                                     81*xi^2/16 - 9*xi/8 - 27/16
                                    -81*xi^2/16 - 9*xi/8 + 27/16
                                     27*xi^2/16 + 9*xi/8 - 1/16];

    %# se calcula el vector de fuerzas distribuidas en los nodos
    te = zeros(1,2*nno);
    te(reshape([2*idx-1; 2*idx],8,1)) = carga(:);

    %# cálculo de la integral:
    suma   = zeros(2*nno, 2*nno);
    N      = zeros(nno,1);
    dN_dxi = zeros(nno,1);
    for p = 1:n_gl
        %# se evalúan las funciones de forma
        N(idx) = NN(xi_gl(p));
        %matN = np.empty((2,2*nno))
        for i =1:nno
            Nijk(:,[2*i-1, 2*i]) = [N(i), 0   ;
                                    0    , N(i)];
        end

        %# se calcula el jacobiano
        dN_dxi(idx) = dNN_dxi(xi_gl(p));
        dx_dxi      = dN_dxi'*xnod(:,X);
        dy_dxi      = dN_dxi'*xnod(:,Y);
        ds_dxi      = hypot(dx_dxi, dy_dxi);

        %# y se calcula la sumatoria
        suma = suma + Nijk'*Nijk*ds_dxi*w_gl(p);
    end
        

    %# finalmente, se retorna el vector de fuerzas nodales equivalentes
 fte = espesor*suma*te';