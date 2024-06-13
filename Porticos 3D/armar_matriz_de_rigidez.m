clc
clear
syms Acy Acz EIz EIy L la lb EA GJ

%% EE
%{
AcyEIz12=(12*Acy*EIz)/(L*(Acy*L^2 + 12*EIz));
AcyEIz6= (6*Acy*EIz)/(Acy*L^2 + 12*EIz);
AcyEIz4= (4*EIz*(Acy*L^2 + 3*EIz))/(L*(Acy*L^2 + 12*EIz));
AcyEIz2= (2*EIz*(- Acy*L^2 + 6*EIz))/(L*(Acy*L^2 + 12*EIz));

AczEIy12=(12*Acz*EIy)/(L*(Acz*L^2 + 12*EIy));
AczEIy6= (6*Acz*EIy)/(Acz*L^2 + 12*EIy);
AczEIy4= (4*EIy*(Acz*L^2 + 3*EIy))/(L*(Acz*L^2 + 12*EIy));
AczEIy2= (2*EIy*(- Acz*L^2 + 6*EIy))/(L*(Acz*L^2 + 12*EIy));
%}      
%%NR  
AcyEIz12=(12*Acy*EIz)/((la - L + lb)*(Acy*L^2 - 2*Acy*L*la - 2*Acy*L*lb + Acy*la^2 + 2*Acy*la*lb + Acy*lb^2 + 12*EIz));
AcyEIz6=(6*Acy*EIz*(L + la - lb))/((la - L + lb)*(Acy*L^2 - 2*Acy*L*la - 2*Acy*L*lb + Acy*la^2 + 2*Acy*la*lb + Acy*lb^2 + 12*EIz));
AcyEIz4=(4*EIz*(Acy*L^2 + Acy*L*la - 2*Acy*L*lb + Acy*la^2 - Acy*la*lb + Acy*lb^2 + 3*EIz))/((la - L + lb)*(Acy*L^2 - 2*Acy*L*la - 2*Acy*L*lb + Acy*la^2 + 2*Acy*la*lb + Acy*lb^2 + 12*EIz));
AcyEIz2=(2*EIz*(Acy*L^2 + Acy*L*la + Acy*L*lb - 2*Acy*la^2 + 2*Acy*la*lb - 2*Acy*lb^2 - 6*EIz))/((la - L + lb)*(Acy*L^2 - 2*Acy*L*la - 2*Acy*L*lb + Acy*la^2 + 2*Acy*la*lb + Acy*lb^2 + 12*EIz));

AczEIy12=(12*Acz*EIy)/((la - L + lb)*(Acz*L^2 - 2*Acz*L*la - 2*Acz*L*lb + Acz*la^2 + 2*Acz*la*lb + Acz*lb^2 + 12*EIy));
AczEIy6= (6*Acz*EIy*(L + la - lb))/((la - L + lb)*(Acz*L^2 - 2*Acz*L*la - 2*Acz*L*lb + Acz*la^2 + 2*Acz*la*lb + Acz*lb^2 + 12*EIy));
AczEIy4= (4*EIy*(Acz*L^2 + Acz*L*la - 2*Acz*L*lb + Acz*la^2 - Acz*la*lb + Acz*lb^2 + 3*EIy))/((la - L + lb)*(Acz*L^2 - 2*Acz*L*la - 2*Acz*L*lb + Acz*la^2 + 2*Acz*la*lb + Acz*lb^2 + 12*EIy));
AczEIy2= (2*EIy*(Acz*L^2 + Acz*L*la + Acz*L*lb - 2*Acz*la^2 + 2*Acz*la*lb - 2*Acz*lb^2 - 6*EIy))/((la - L + lb)*(Acz*L^2 - 2*Acz*L*la - 2*Acz*L*lb + Acz*la^2 + 2*Acz*la*lb + Acz*lb^2 + 12*EIy));

Kloc=[[  EA/L,         0,         0,     0,       0,        0, -EA/L,        0,         0,     0,        0,       0];
      [     0,  AcyEIz12,         0,     0,       0,  AcyEIz6,     0,-AcyEIz12,         0,     0,        0, AcyEIz6];
      [     0,         0,  AczEIy12,     0,-AczEIy6,        0,     0,        0, -AczEIy12,     0, -AczEIy6,       0];
      [     0,         0,         0,  GJ/L,       0,        0,     0,        0,         0, -GJ/L,        0,       0];
      [     0,         0,  -AczEIy6,     0, AczEIy4,        0,     0,        0,   AczEIy6,     0, -AczEIy2,       0];
      [     0,   AcyEIz6,         0,     0,       0,  AcyEIz4,     0, -AcyEIz6,         0,     0,        0,-AcyEIz2];
      [ -EA/L,         0,         0,     0,       0,        0,  EA/L,        0,         0,     0,        0,       0];
      [     0, -AcyEIz12,         0,     0,       0, -AcyEIz6,     0, AcyEIz12,         0,     0,        0,-AcyEIz6];
      [     0,         0, -AczEIy12,     0, AczEIy6,        0,     0,        0,  AczEIy12,     0,  AczEIy6,       0];
      [     0,         0,         0, -GJ/L,       0,        0,     0,        0,         0,  GJ/L,        0,       0];
      [     0,         0,  -AczEIy6,     0,-AczEIy2,        0,     0,        0,   AczEIy6,     0,  AczEIy4,       0];
      [     0,   AcyEIz6,         0,     0,       0, -AcyEIz2,     0, -AcyEIz6,         0,     0,        0, AcyEIz4]];