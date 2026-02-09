%% --- PARÁMETROS CORREGIDOS ---
E = 21525560.5; 
G = 0.4 * E; 
b_sec = 0.2; h_sec = 0.4;
A = b_sec * h_sec;
I = (b_sec * h_sec^3) / 12;
EI = E * I;
EA = E * A;
Ac = (5/6) * G * A; % Rigidez a cortante (sección rectangular)

hx =(x1^2*y2 - x2^2*y1 - x1^2*y3 + x3^2*y1 + x2^2*y3 - x3^2*y2)/(2*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
ky =(x1^4*y2^2 - 2*x1^4*y2*y3 + x1^4*y3^2 + 4*x1^3*x2*y2*y3 - 4*x1^3*x2*y3^2 - 4*x1^3*x3*y2^2 + 4*x1^3*x3*y2*y3 - 2*x1^2*x2^2*y1*y2 - 2*x1^2*x2^2*y1*y3 - 2*x1^2*x2^2*y2*y3 + 6*x1^2*x2^2*y3^2 + 4*x1^2*x2*x3*y1*y2 + 4*x1^2*x2*x3*y1*y3 - 8*x1^2*x2*x3*y2*y3 - 2*x1^2*x3^2*y1*y2 - 2*x1^2*x3^2*y1*y3 + 6*x1^2*x3^2*y2^2 - 2*x1^2*x3^2*y2*y3 + 4*x1*x2^3*y1*y3 - 4*x1*x2^3*y3^2 + 4*x1*x2^2*x3*y1*y2 - 8*x1*x2^2*x3*y1*y3 + 4*x1*x2^2*x3*y2*y3 - 8*x1*x2*x3^2*y1*y2 + 4*x1*x2*x3^2*y1*y3 + 4*x1*x2*x3^2*y2*y3 + 4*x1*x3^3*y1*y2 - 4*x1*x3^3*y2^2 + x2^4*y1^2 - 2*x2^4*y1*y3 + x2^4*y3^2 - 4*x2^3*x3*y1^2 + 4*x2^3*x3*y1*y3 + 6*x2^2*x3^2*y1^2 - 2*x2^2*x3^2*y1*y2 - 2*x2^2*x3^2*y1*y3 - 2*x2^2*x3^2*y2*y3 - 4*x2*x3^3*y1^2 + 4*x2*x3^3*y1*y2 + x3^4*y1^2 - 2*x3^4*y1*y2 + x3^4*y2^2)/(4*(x1 - x2)*(x1 - x3)*(x2 - x3)*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
p =((x1 - x2)*(x1 - x3)*(x2 - x3))/(4*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
 
yp=(x-hx)^2/(-4*p)+ky;
focoy=ky-p;
focox=hx;
dydx = diff(yp, x);
ds = sqrt(1 + dydx^2);

% Funciones de proyección trigonométrica local
cos_th = 1 / ds;       % dx/ds
sin_th = dydx / ds;    % dy/ds

%% --- 1. MATRIZ DE FLEXIBILIDAD (f_ij) ---
% Representa los desplazamientos en el apoyo liberado ante cargas unitarias

% f11: Desplazamiento vertical por Qy=1
f11 = (x^2/EI + sin_th^2/EA + cos_th^2/Ac)*ds;
% f22: Desplazamiento horizontal por Qx=1
f22 = (yp^2/EI + cos_th^2/EA + sin_th^2/Ac)*ds;
% f33: Giro por M=1
f33 = (1/EI)*ds;
% f12: Cruce Vertical-Horizontal (Simetría)
f12 = (x*yp/EI + sin_th*cos_th/EA - sin_th*cos_th/Ac)*ds;
% f13: Cruce Vertical-Giro
f13 = (x/EI)*ds;
% f23: Cruce Horizontal-Giro
f23 = (yp/EI)*ds;

% Integración numérica
AMxx = area_cuadraturas(a, b, matlabFunction(f11, 'Vars', x));
AMyy = area_cuadraturas(a, b, matlabFunction(f22, 'Vars', x));
AMgg = area_cuadraturas(a, b, matlabFunction(f33, 'Vars', x));
AMxy = area_cuadraturas(a, b, matlabFunction(f12, 'Vars', x));
AMxg = area_cuadraturas(a, b, matlabFunction(f13, 'Vars', x));
AMyg = area_cuadraturas(a, b, matlabFunction(f23, 'Vars', x));

% Construcción de la matriz (Simétrica)
ec = [AMxx, AMxy, AMxg;
      AMxy, AMyy, AMyg;
      AMxg, AMyg, AMgg];

%% --- 2. VECTOR DE CARGA EXTERNA (Deflexiones isostáticas) ---
% Proyectamos la carga externa vertical en componentes locales
Ne = V * sin_th; % Axial externo
Ve = V * cos_th; % Cortante externo

f_Mx = (M*x/EI + Ne*sin_th/EA + Ve*cos_th/Ac)*ds;
f_My = (M*yp/EI + Ne*cos_th/EA - Ve*sin_th/Ac)*ds;
f_Mg = (M/EI)*ds;

AMMx = area_cuadraturas(a, b, matlabFunction(f_Mx, 'Vars', x));
AMMy = area_cuadraturas(a, b, matlabFunction(f_My, 'Vars', x));
AMMg = area_cuadraturas(a, b, matlabFunction(f_Mg, 'Vars', x));

cein = [-AMMx; -AMMy; -AMMg];

%% --- 3. SOLUCIÓN ---
sol = ec \ cein; % [Ray; Rax; Ma]