function [Hy1,Hy2,Hy3,Hy4,Hy5,Hy6,Hy7,Hy8,Hy9]=Hy(xi)

    Hy1 = 1.5*(d6*N6 - d5*N5);
    Hy2 = -N1 + e5*N5 + e6*N6;
    Hy3 = -b5*N5 - b6*N6;
    Hy4 = 1.5*(d4*N4 - d6*N6);
    Hy5 = -N2 + e4*N4 + e6*N6;
    Hy6 = -b4*N4 - b6*N6;
    Hy7 = 1.5*(d5*N5 - d4*N4);
    Hy8 = -N3 + e4*N4 + e5*N5;
    Hy9 = -b4*N4 - b5*N5;