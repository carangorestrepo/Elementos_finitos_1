function [Hx1,Hx2,Hx3,Hx4,Hx5,Hx6,Hx7,Hx8,Hx9]=Hx(xi)

    Hx1 = 1.5*(a6*N6 - a5*N5);
    Hx2 = b5*N5 + b6*N6;
    Hx3 = N1 - c5*N5 - c6*N6;
    Hx4 = 1.5*(a4*N4 - a6*N6);
    Hx5 = b4*N4 + b6*N6;
    Hx6 = N2 - c4*N4 - c6*N6;
    Hx7 = 1.5*(a5*N5 - a4*N4);
    Hx8 = b4*N4 + b5*N5;
    Hx9 = N3 - c4*N4 - c5*N5;
