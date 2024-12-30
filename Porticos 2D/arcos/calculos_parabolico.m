 function R=calculos_parabolico(y,EI,M,a,b,L,h,qa,qb,DX,DY,GM)
syms x

%% Longitud funcion de arco  
dxdym=sqrt(1+diff(y,x)^2);

%% C�lculo de deformaciones
% Deformaci�n vertical por cargas externas
AMx = matlabFunction(M * x * dxdym / EI, "Vars", x);
% Deformaci�n vertical por carga unitaria vertical
Axx = matlabFunction(x * x * dxdym / EI, "Vars", x);
% Deformaci�n horizontal por carga unitaria vertical
Axy = matlabFunction(x * y * dxdym / EI, "Vars", x);
% �ngulo de giro por carga unitaria vertical
Ax =  matlabFunction(x * dxdym / EI, "Vars", x);
% C�lculo de �reas mediante cuadraturas (integrales num�ricas)
AMx = area_cuadraturas(a, b, AMx); % Deformaci�n vertical positiva (Qy)
Axx = area_cuadraturas(a, b, Axx); % Deformaci�n vertical puntual negativa (Py)
Axy = area_cuadraturas(a, b, Axy); % Deformaci�n horizontal puntual negativa (Py)
Ax =  area_cuadraturas(a, b, Ax);  % �ngulo de giro puntual (Py)
%% Deformaci�n horizontal
% Deformaci�n horizontal por cargas externas
AMy = matlabFunction(M * y * dxdym / EI, "Vars", x);
% Deformaci�n horizontal por carga unitaria horizontal
Ayy = matlabFunction(y * y * dxdym / EI, "Vars", x);
% �ngulo de giro por carga unitaria horizontal
Ay = matlabFunction(y * dxdym / EI, "Vars", x);
% Deformaci�n vertical por carga unitaria horizontal
Ayx = matlabFunction(y * x * dxdym / EI, "Vars", x);
% C�lculo de �reas mediante cuadraturas para deformaciones horizontales
AMy = area_cuadraturas(a, b, AMy);  % Deformaci�n horizontal por carga (Qy)
Ayx = area_cuadraturas(a, b, Ayx);  % Deformaci�n vertical por carga horizontal (Px)
Ayy = area_cuadraturas(a, b, Ayy);  % Deformaci�n horizontal por carga horizontal (Px)
Ay =  area_cuadraturas(a, b, Ay);   % �ngulo de giro por carga horizontal (Px)
%% �ngulo de giro
% C�lculo del �ngulo de giro debido a las cargas externas
AM = matlabFunction(M * dxdym / EI, "Vars", x);
% �ngulo de giro por momento unitario
Adx = matlabFunction(1 * dxdym / EI, "Vars", x);
% Deformaci�n horizontal por momento unitario
Ayg = matlabFunction(y * dxdym / EI, "Vars", x);
% Deformaci�n vertical por momento unitario
Axg = matlabFunction(x * dxdym / EI, "Vars", x);
% C�lculo de �reas mediante cuadraturas para el �ngulo de giro
AM =  area_cuadraturas(a, b, AM);   % �ngulo de giro debido a Qy
Ayg = area_cuadraturas(a, b, Ayg);  % Deformaci�n horizontal debido a momento M
Axg = area_cuadraturas(a, b, Axg);  % Deformaci�n vertical debido a momento M
Adx = area_cuadraturas(a, b, Adx);  % �ngulo de giro puntual debido a momento M
%% Ecuaciones para resolver el sistema
% Las ecuaciones corresponden a las sumas de fuerzas y momentos en la viga
%ec1 = AxP1   * ryy + AxPy1 * rxx + AxM1 + xAP1 * mm == 0;  % Suma de fuerzas en Y
%ec2 = AxyP1  * ryy + AyP1  * rxx + AyM1 + yAP1 * mm == 0;  % Suma de fuerzas en X
%ec3 = G1yAP  * ryy +G1xAP1 * rxx +  G1  + GP1  * mm == 0;  % Suma de momentos
% Resoluci�n del sistema de ecuaciones
%[Rxx, Ryy, Mm] = solve(ec1, ec2, ec3, [rxx, ryy, mm]);
   %rxx   ryy    mm
ec=[Axx ,Axy ,Ax;  % Suma de fuerzas en Y
    Ayx ,Ayy ,Ay;  % Suma de fuerzas en X
    Axg ,Ayg ,Adx];% Suma de momentos
cein=[-AMx-DY;-AMy-DX;-AM-GM];
sol=ec\cein;
Ma=sol(3);
Ray=-sol(1);
Rax=sol(2);
Rby=qa*L - (L^2*(qa - qb))/(2*L)-Ray;
Rbx=-Rax;
Mb=-Rax*h-Ma-(-((qa - qb) * L^3) / (6 * L) + (qa * L^2) / 2)+Ray*L;
R=[Rax;Ray;Ma;Rbx;Rby;Mb];
if a<0
    Ma=-Ma;
    Rax=-Rax;
    Ray=-Ray;
    Mb=-(-Rax*h+Ma-(-((qa - qb) * L^3) / (6 * L) + (qa * L^2) / 2)-Ray*L);
    R=[-Rbx;-Rby;Mb;Rax;Ray;Ma];
end


