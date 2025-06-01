function [Phi, Omega, T] = estodola_1(k, m, iteraciones, modos)
% ESTODOLA Calcula valores y vectores propios de un sistema din�mico
% 
% Entradas:
%   k: Matriz de rigidez del sistema (matriz cuadrada)
%   m: Matriz de masa del sistema (matriz cuadrada del mismo tama�o que k)
%   iteraciones: N�mero de iteraciones para convergencia de cada modo
%   modos: N�mero de modos de vibraci�n a calcular
%
% Salidas:
%   Phi: Matriz de modos de vibraci�n (vectores propios, columnas normalizadas)
%   Omega: Vector de frecuencias naturales (valores propios, orden ascendente)
%   T: Vector de periodos de vibraci�n correspondientes

    % Obtener tama�o del sistema (n�mero de grados de libertad)
    sizek = size(k, 1);
    
    % Inicializar matriz de deflaci�n como matriz identidad
    Bo = diag(ones(1, sizek));
    
    % Prealocaci�n de memoria para almacenamiento de resultados intermedios
    x1n = cell(modos, 1);    % Celdas para historial de vectores iterados
    xm1n = cell(modos, 1);   % Celdas para historial de vectores transformados
    Phi = zeros(sizek, modos); % Matriz para modos de vibraci�n finales
    Omega = zeros(modos, 1);   % Vector para frecuencias naturales
    
    % Bucle principal para calcular cada modo
    for n = 1:modos
        % Matriz de transformaci�n para el modo actual (con deflaci�n)
        Fb = k^-1 * m * Bo;
        
        % Inicializar matrices para el proceso iterativo
        x1 = ones(sizek, iteraciones-1);   % Vectores iterados
        xm1 = zeros(sizek, iteraciones-1); % Vectores transformados
        
        % Proceso iterativo para converger al modo actual
        for i = 2:iteraciones
            % Aplicar transformaci�n y normalizar
            xm1(:, i-1) = Fb * x1(:, i-1);              % Transformaci�n
            x1(:, i) = xm1(:, i-1) / max(abs(xm1(:, i-1))); % Normalizaci�n
        end
        
        % Almacenar historial de iteraciones (�til para depuraci�n)
        x1n{n} = x1;
        xm1n{n} = xm1;
        
        % Actualizar matriz de deflaci�n para eliminar modo encontrado
        Bo = Bo - x1(:, i) * x1(:, i)' * m / (x1(:, i)' * m * x1(:, i));
        
        % Almacenar modo convergido (�ltima iteraci�n)
        Phi(:, n) = x1(:, i-1);
        
        % Calcular frecuencia natural correspondiente (? = sqrt(?))
        Omega(n) = sqrt(abs(1 / max(abs(xm1(:, i-1)))));
    end
    
    % Ordenar resultados por frecuencia ascendente
    [Omega, In] = sort(Omega);     % Ordena frecuencias y obtiene �ndices
    T = 2 * pi ./ Omega;           % Calcula periodos de vibraci�n (T = 2?/?)
    Phi = Phi(:, In);              % Reordena modos seg�n frecuencias
    
    % Normalizar modos con respecto a la matriz de masa (??M? = I)
    Phi = Phi ./ repmat(sqrt(diag(Phi' * m * Phi)'), size(Phi, 1), 1);
end