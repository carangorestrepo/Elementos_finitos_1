clear
clc
format long g

%% ============================================================
% MATRICES DEL PROBLEMA
% ============================================================

K = 1e3*[ ...
     76.8690,  -99.6910,   25.5830,   -3.7470;
    -99.6910,  209.1400, -136.0200,   31.1080;
     25.5830, -136.0200,  221.7600, -142.1100;
     -3.4747,   31.1080, -142.1100,  252.1000];

M = [ ...
    60,  0,  0,  0;
     0, 60,  0,  0;
     0,  0, 60,  0;
     0,  0,  0, 60];

%% ============================================================
% REVISAR SIMETRIA
% ============================================================

fprintf('Error de simetria de K = %.12e\n',norm(K-K.','fro'));

% Para este ejemplo se fuerza simetria
K = 0.5*(K + K.');

fprintf('Error despues de simetrizar = %.12e\n\n', ...
    norm(K-K.','fro'));

%% ============================================================
% NUMERO DE VECTORES DE LANCZOS
% ============================================================

n = size(K,1);

% Para este ejemplo pequeño podemos usar todo el espacio
mLanczos = n;

%% ============================================================
% VECTOR INICIAL
% ============================================================

f = ones(n,1);

%% ============================================================
% LANCZOS GENERALIZADO
%
% Problema:
%
%       K*phi = lambda*M*phi
%
% Usamos el operador:
%
%       A = K^(-1)*M
%
% cuyos autovalores son:
%
%       mu = 1/lambda
%
% ============================================================

[Q,T,alpha,beta] = lanczos_KinvM(K,M,f,mLanczos);

%% ============================================================
% PROBLEMA PROPIO REDUCIDO
% ============================================================

[S,Dmu] = eig(T);

mu = diag(Dmu);

% Como:
%
%       mu = 1/lambda
%
lambdaLanczos = 1./mu;

%% Ordenar de menor a mayor lambda

[lambdaLanczos,ind] = sort(lambdaLanczos);

S = S(:,ind);

%% ============================================================
% RECUPERAR VECTORES PROPIOS
%
%       Phi = Q*S
% ============================================================

PhiLanczos = Q*S;

%% ============================================================
% NORMALIZACION RESPECTO A LA MASA
%
%       phi_i' * M * phi_i = 1
% ============================================================

for i = 1:size(PhiLanczos,2)

    masaModal = PhiLanczos(:,i)'*M*PhiLanczos(:,i);

    PhiLanczos(:,i) = ...
        PhiLanczos(:,i)/sqrt(masaModal);

end

%% ============================================================
% FRECUENCIAS
% ============================================================

omegaLanczos = sqrt(lambdaLanczos);

freqLanczos = omegaLanczos/(2*pi);

%% ============================================================
% SOLUCION DIRECTA CON MATLAB
% PARA COMPARACION
% ============================================================

[PhiExacta,DExacta] = eig(K,M);

lambdaExacta = diag(DExacta);

[lambdaExacta,ind] = sort(lambdaExacta);

PhiExacta = PhiExacta(:,ind);

%% Normalización respecto a M

for i = 1:size(PhiExacta,2)

    masaModal = PhiExacta(:,i)'*M*PhiExacta(:,i);

    PhiExacta(:,i) = ...
        PhiExacta(:,i)/sqrt(masaModal);

end

omegaExacta = sqrt(lambdaExacta);

freqExacta = omegaExacta/(2*pi);

%% ============================================================
% RESULTADOS
% ============================================================

fprintf('\n====================================================\n');
fprintf('          VALORES PROPIOS\n');
fprintf('====================================================\n\n');

fprintf('Modo       Lanczos              eig(K,M)\n');

for i = 1:length(lambdaLanczos)

    fprintf('%3d   %18.10f   %18.10f\n', ...
        i,lambdaLanczos(i),lambdaExacta(i));

end

fprintf('\n====================================================\n');
fprintf('          FRECUENCIAS NATURALES [Hz]\n');
fprintf('====================================================\n\n');

fprintf('Modo       Lanczos              Exacta\n');

for i = 1:length(freqLanczos)

    fprintf('%3d   %18.10f   %18.10f\n', ...
        i,freqLanczos(i),freqExacta(i));

end

%% ============================================================
% COMPROBAR ORTOGONALIDAD MODAL
% ============================================================

Mmodal = PhiLanczos'*M*PhiLanczos;

Kmodal = PhiLanczos'*K*PhiLanczos;

fprintf('\nPhi''*M*Phi =\n');
disp(Mmodal)

fprintf('\nPhi''*K*Phi =\n');
disp(Kmodal)

%% ============================================================
% RESIDUOS
%
%       r_i = K*phi_i - lambda_i*M*phi_i
%
% ============================================================

fprintf('\n====================================================\n');
fprintf('          RESIDUOS LANCZOS\n');
fprintf('====================================================\n\n');

for i = 1:length(lambdaLanczos)

    r = K*PhiLanczos(:,i) ...
        - lambdaLanczos(i)*M*PhiLanczos(:,i);

    errorRelativo = norm(r)/( ...
        norm(K*PhiLanczos(:,i)) + eps);

    fprintf('Modo %d   residuo relativo = %.12e\n', ...
        i,errorRelativo);

end

%% ============================================================
% MOSTRAR MATRIZ TRIDIAGONAL DE LANCZOS
% ============================================================

fprintf('\nMatriz T de Lanczos:\n');
disp(T)

fprintf('\nCoeficientes alpha:\n');
disp(alpha)

fprintf('\nCoeficientes beta:\n');
disp(beta)


%% ============================================================
% FUNCION LANCZOS
% ============================================================

function [Q,T,alpha,beta] = lanczos_KinvM(K,M,f,m)
%LANCZOS_KINVM
%
% Implementación del método de Lanczos para:
%
%           K*phi = lambda*M*phi
%
% utilizando:
%
%           A = inv(K)*M
%
% sin formar explícitamente inv(K):
%
%           A*q = K\(M*q)
%
% Los autovalores del operador son:
%
%           mu = 1/lambda
%
% ------------------------------------------------------------
% ENTRADAS
%
% K : matriz de rigidez
% M : matriz de masa
% f : vector inicial
% m : número máximo de vectores de Lanczos
%
% ------------------------------------------------------------
% SALIDAS
%
% Q     : base de Lanczos
% T     : matriz tridiagonal proyectada
% alpha : diagonal de T
% beta  : sub/superdiagonal de T
%
% ------------------------------------------------------------

n = size(K,1);

m = min(m,n);

Q = zeros(n,m);

alpha = zeros(m,1);

beta = zeros(max(m-1,1),1);

%% ============================================================
% NORMALIZACION DEL VECTOR INICIAL
%
% Para A=K^-1 M, el operador es auto-adjunto respecto
% al producto interno asociado a M:
%
%          <x,y>_M = x' M y
%
% ============================================================

norma = sqrt(f'*M*f);

if norma < eps
    error('El vector inicial f no puede ser nulo.');
end

q = f/norma;

qAnterior = zeros(n,1);

betaAnterior = 0;

%% ============================================================
% ITERACION DE LANCZOS
% ============================================================

mReal = m;

for j = 1:m

    Q(:,j) = q;

    % --------------------------------------------------------
    % Aplicación del operador:
    %
    %       z = K^-1 M q
    %
    % Nunca utilizar:
    %
    %       inv(K)*M*q
    %
    % --------------------------------------------------------

    z = K\(M*q);

    % --------------------------------------------------------
    % Eliminar contribución del vector anterior
    % --------------------------------------------------------

    if j > 1

        z = z - betaAnterior*qAnterior;

    end

    % --------------------------------------------------------
    % Coeficiente diagonal
    %
    % Producto interno M
    % --------------------------------------------------------

    alpha(j) = q'*M*z;

    % --------------------------------------------------------
    % Ortogonalización respecto al vector actual
    % --------------------------------------------------------

    z = z - alpha(j)*q;

    % --------------------------------------------------------
    % Reortogonalización completa
    %
    % Numéricamente es recomendable aunque Lanczos teórico
    % utiliza solamente la recurrencia de tres términos.
    % --------------------------------------------------------

    for k = 1:j

        correccion = Q(:,k)'*M*z;

        z = z - Q(:,k)*correccion;

    end

    %% Última iteración

    if j == m
        break
    end

    % --------------------------------------------------------
    % Coeficiente beta
    % --------------------------------------------------------

    beta(j) = sqrt(z'*M*z);

    % --------------------------------------------------------
    % Verificar ruptura de Lanczos
    % --------------------------------------------------------

    if beta(j) < 1e-12

        mReal = j;

        break

    end

    % --------------------------------------------------------
    % Preparar siguiente iteración
    % --------------------------------------------------------

    qAnterior = q;

    q = z/beta(j);

    betaAnterior = beta(j);

end

%% ============================================================
% RECORTAR SI HUBO CONVERGENCIA PREMATURA
% ============================================================

Q = Q(:,1:mReal);

alpha = alpha(1:mReal);

if mReal > 1
    beta = beta(1:mReal-1);
else
    beta = [];
end

%% ============================================================
% MATRIZ TRIDIAGONAL
% ============================================================

T = diag(alpha);

if mReal > 1

    T = T ...
        + diag(beta,1) ...
        + diag(beta,-1);

end

end