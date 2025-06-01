function [eigenvalues, eigenvectors] = jacobi_method(A, tol, max_iter)
    % Verificar que la matriz es sim�trica
    if ~isequal(A, A')
        error('La matriz debe ser sim�trica.');
    end
    
    % Dimensi�n de la matriz
    n = size(A, 1);
    
    % Inicializaci�n de la matriz de autovectores
    V = eye(n);
    
    % Iteraci�n del m�todo de Jacobi
    for iter = 1:max_iter
        % Encontrar el elemento fuera de la diagonal con mayor valor absoluto
        [max_val, p, q] = max_off_diagonal(A);
        
        % Criterio de convergencia
        if max_val < tol
            break;
        end
        
        % Calcular el �ngulo de rotaci�n
        theta = 0.5 * atan2(2 * A(p, q), A(p, p) - A(q, q));
        
        % Construir la matriz de rotaci�n
        R = eye(n);
        R([p, q], [p, q]) = [cos(theta), -sin(theta); sin(theta), cos(theta)];
        
        % Actualizar la matriz A y los autovectores
        A = R' * A * R;
        V = V * R;
    end
    
    % Extraer autovalores y autovectores
    eigenvalues = diag(A);
    eigenvectors = V;
end

function [max_val, p, q] = max_off_diagonal(A)
    % Encuentra el mayor valor absoluto fuera de la diagonal
    [n, ~] = size(A);
    max_val = 0;
    p = 1; q = 2;
    for i = 1:n-1
        for j = i+1:n
            if abs(A(i, j)) > max_val
                max_val = abs(A(i, j));
                p = i; q = j;
            end
        end
    end
end