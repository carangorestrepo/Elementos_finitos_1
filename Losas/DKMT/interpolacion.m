

syms  gtz4 gtz5 gtz6 xi eta

gxiz=  gtz4-(gtz4+(2)^(1/2)*gtz5-gtz6)*eta;
getaz=-gtz6+(gtz4+(2)^(1/2)*gtz5-gtz6)*xi;

gns=[gtz4 gtz5 gtz6].';
Ngama = equationsToMatrix([gxiz;getaz], gns)