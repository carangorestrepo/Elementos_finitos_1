syms x1 y1 x2 y2 x3 y3 x y
lcx = 0.35;   % Largo columna [m]
lcy = 0.65;   % Ancho columna [m]

cx0 = 1.1;   % Posición x columna [m]
cy0 = 1;     % Posición y columna [m]

Mx=200;
My=150;
P=-1000;
ex=Mx/P;
ey=My/P;
Z =(pm*(lcx^2*lcy^2+ 12*cx0*lcy^2*abs(ex)+ 12*cy0*lcx^2*abs(ey)+ 6*lcx*lcy^2*abs(ex)+ 6*lcx^2*lcy*abs(ey)- 12*lcy^2*x*abs(ex)- 12*lcx^2*y*abs(ey)))/(lcx^2*lcy^2);


x1=2.05;
y1=1.4;

x2=2.05;
y2=1.4;

x3=2.05;
y3=1.45;



% Rectas entre pares de puntos (pendiente-intersección)
L12 = (y2 - y1)/(x2 - x1)*(x - x1) + y1;
L13 = (y3 - y1)/(x3 - x1)*(x - x1) + y1;
L23 = (y3 - y2)/(x3 - x2)*(x - x2) + y2;

% Asumimos x1 < x2 < x3 para definir las regiones correctamente
A1 = int(int(1, y, min(L12, L13), max(L12, L13)), x, x1, x2);
A2 = int(int(1, y, min(L13, L23), max(L13, L23)), x, x2, x3);

% Área total (valor absoluto para asegurar área positiva)
Area = simplify(abs(A1 + A2), 'Steps', 100);
