function [Phi, Omega, T] = estodola_2(k, m, iteraciones, modos)
% ESTODOLA_1
%
% Calcula modos de vibración, frecuencias naturales y periodos
% mediante iteración inversa de Stodola con ortogonalización modal.
%
% El problema dinámico es:
%
%               K*phi = lambda*M*phi
%
% donde:
%
%               lambda = omega^2
%
% Por tanto:
%
%               omega = sqrt(lambda)     [rad/s]
%               f     = omega/(2*pi)     [Hz]
%               T     = 2*pi/omega       [s]
%
% -------------------------------------------------------------------------
% ENTRADAS
% -------------------------------------------------------------------------
%
%   k           : Matriz de rigidez K
%   m           : Matriz de masa M
%   iteraciones : Número máximo de iteraciones para cada modo
%   modos       : Número de modos que se desean calcular
%
% -------------------------------------------------------------------------
% SALIDAS
% -------------------------------------------------------------------------
%
%   Phi         : Modos de vibración normalizados respecto a M
%
%                       Phi'*M*Phi = I
%
%   Omega       : Frecuencias circulares naturales [rad/s]
%
%   T           : Periodos naturales [s]
%
% -------------------------------------------------------------------------

%% ========================================================================
% 1. TAMAÑO DEL SISTEMA
% =========================================================================

sizek = size(k,1);

%% ========================================================================
% 2. COMPROBACIONES
% =========================================================================

if size(k,1) ~= size(k,2)
    error('La matriz de rigidez K debe ser cuadrada.')
end

if size(m,1) ~= size(m,2)
    error('La matriz de masa M debe ser cuadrada.')
end

if size(k,1) ~= size(m,1)
    error('Las matrices K y M deben tener el mismo tamaño.')
end

if modos > sizek
    error('El número de modos no puede superar el número de grados de libertad.')
end

% Eliminar pequeñas asimetrías numéricas
k = 0.5*(k+k');
m = 0.5*(m+m');

%% ========================================================================
% 3. INICIALIZACIÓN
% =========================================================================

% Matriz de modos
Phi = zeros(sizek,modos);

% Frecuencias naturales
Omega = zeros(modos,1);

% Valores propios
Lambda = zeros(modos,1);

% Historial de iteraciones
x1n  = cell(modos,1);
xm1n = cell(modos,1);

% Tolerancia de convergencia
tolerancia = 1e-10;

%% ========================================================================
% 4. CÁLCULO DE LOS MODOS
% =========================================================================

for n = 1:modos

    %% --------------------------------------------------------------------
    % Vector inicial
    % ---------------------------------------------------------------------

    x1 = ones(sizek,iteraciones);

    % Se modifica ligeramente el vector inicial para evitar problemas
    % cuando ciertos modos son ortogonales al vector [1 1 ... 1].
    x1(:,1) = ones(sizek,1);

    % Para modos superiores puede ser conveniente variar el vector inicial
    if n > 1
        x1(:,1) = (1:sizek)' + 0.1*n;
    end

    xm1 = zeros(sizek,iteraciones-1);

    %% --------------------------------------------------------------------
    % Ortogonalizar el vector inicial respecto de modos anteriores
    % ---------------------------------------------------------------------

    if n > 1

        for j = 1:n-1

            x1(:,1) = x1(:,1) ...
                - Phi(:,j)*(Phi(:,j)'*m*x1(:,1));

        end

    end

    % Normalización respecto a la matriz de masa
    normaM = sqrt(x1(:,1)'*m*x1(:,1));

    if normaM < eps
        error('Vector inicial degenerado para el modo %d.',n)
    end

    x1(:,1) = x1(:,1)/normaM;

    %% --------------------------------------------------------------------
    % Iteración inversa
    % ---------------------------------------------------------------------

    lambda_anterior = inf;

    for i = 2:iteraciones

        % -------------------------------------------------------------
        % Iteración inversa:
        %
        %       K*y = M*x
        %
        % equivalente a:
        %
        %       y = inv(K)*M*x
        %
        % pero SIN calcular explícitamente inv(K).
        % -------------------------------------------------------------

        xm1(:,i-1) = k \ (m*x1(:,i-1));

        %% -------------------------------------------------------------
        % Deflación / ortogonalización respecto a los modos anteriores
        %
        % Se impone:
        %
        %       phi_j' * M * phi_n = 0
        % -------------------------------------------------------------

        if n > 1

            for j = 1:n-1

                xm1(:,i-1) = xm1(:,i-1) ...
                    - Phi(:,j)*(Phi(:,j)'*m*xm1(:,i-1));

            end

        end

        %% -------------------------------------------------------------
        % Normalización respecto a M
        % -------------------------------------------------------------

        normaM = sqrt(xm1(:,i-1)'*m*xm1(:,i-1));

        if normaM < eps
            error('Se obtuvo un vector nulo durante la iteración.')
        end

        x1(:,i) = xm1(:,i-1)/normaM;

        %% -------------------------------------------------------------
        % Cociente de Rayleigh
        %
        %       lambda = (phi' K phi)/(phi' M phi)
        %
        % Como phi está normalizado respecto a M:
        %
        %       phi'Mphi = 1
        % -------------------------------------------------------------

        lambda_actual = ...
            (x1(:,i)'*k*x1(:,i)) / ...
            (x1(:,i)'*m*x1(:,i));

        %% -------------------------------------------------------------
        % Verificación de convergencia
        % -------------------------------------------------------------

        if abs(lambda_actual-lambda_anterior) ...
                <= tolerancia*max(1,abs(lambda_actual))

            break

        end

        lambda_anterior = lambda_actual;

    end

    %% --------------------------------------------------------------------
    % Guardar historial
    % ---------------------------------------------------------------------

    x1n{n}  = x1(:,1:i);
    xm1n{n} = xm1(:,1:i-1);

    %% --------------------------------------------------------------------
    % Modo convergido
    % ---------------------------------------------------------------------

    phi = x1(:,i);

    %% --------------------------------------------------------------------
    % Nueva ortogonalización para reducir errores acumulados
    % ---------------------------------------------------------------------

    if n > 1

        for j = 1:n-1

            phi = phi ...
                - Phi(:,j)*(Phi(:,j)'*m*phi);

        end

    end

    %% --------------------------------------------------------------------
    % Normalización modal respecto a M
    % ---------------------------------------------------------------------

    phi = phi/sqrt(phi'*m*phi);

    %% --------------------------------------------------------------------
    % Convención de signo
    %
    % El signo de un vector propio es arbitrario.
    % Se hace positivo el componente de mayor valor absoluto.
    % ---------------------------------------------------------------------

    [~,imax] = max(abs(phi));

    if phi(imax) < 0
        phi = -phi;
    end

    %% --------------------------------------------------------------------
    % Valor propio mediante cociente de Rayleigh
    % ---------------------------------------------------------------------

    lambda = (phi'*k*phi)/(phi'*m*phi);

    %% --------------------------------------------------------------------
    % Almacenar resultados
    % ---------------------------------------------------------------------

    Phi(:,n) = phi;

    Lambda(n) = lambda;

    Omega(n) = sqrt(abs(lambda));

end

%% ========================================================================
% 5. ORDENAR LOS MODOS DE MENOR A MAYOR FRECUENCIA
% =========================================================================

[Omega,In] = sort(Omega);

Phi = Phi(:,In);

Lambda = Lambda(In);

%% ========================================================================
% 6. PERIODOS
% =========================================================================

T = 2*pi./Omega;

%% ========================================================================
% 7. NORMALIZACIÓN FINAL RESPECTO A LA MATRIZ DE MASA
% =========================================================================

for n = 1:modos

    Phi(:,n) = Phi(:,n)/sqrt(Phi(:,n)'*m*Phi(:,n));

end

%% ========================================================================
% 8. INFORMACIÓN DE CONTROL
% =========================================================================
%
% Debería cumplirse aproximadamente:
%
%       Phi' * M * Phi = I
%
% y
%
%       K*Phi(:,i) = Omega(i)^2*M*Phi(:,i)
%
% -------------------------------------------------------------------------

end