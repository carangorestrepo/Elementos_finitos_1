function [Nw,Ntx,Nty,dNw_dx,dNw_dy,dNtx_dx,dNtx_dy,dNty_dx,dNty_dy] = ...
    funciones_forma_Mindlin_TE(xi,eta,Lx,Ly,EI,Ks)
%==========================================================================
% FUNCIONES DE FORMA MINDLIN A PARTIR DE FUNCIONES EXACTAS DE TIMOSHENKO
% Elemento rectangular alineado con ejes x-y.
% GDL = [w1 tx1 ty1 w2 tx2 ty2 w3 tx3 ty3 w4 tx4 ty4]
%==========================================================================

rx = (1+xi)/2;
ry = (1+eta)/2;
sx = rx*Lx;
sy = ry*Ly;

[Nvx,Ntx1,dNvx,dNtx1] = funciones_Timoshenko_1D(sx,Lx,EI,Ks);
[Nvy,Nty1,dNvy,dNty1] = funciones_Timoshenko_1D(sy,Ly,EI,Ks);

T12 = zeros(4,12);
T12(1,1)=1; T12(2,2)=1; T12(3,4)=1; T12(4,5)=1;

T43 = zeros(4,12);
T43(1,10)=1; T43(2,11)=1; T43(3,7)=1; T43(4,8)=1;

T14 = zeros(4,12);
T14(1,1)=1; T14(2,3)=1; T14(3,10)=1; T14(4,12)=1;

T23 = zeros(4,12);
T23(1,4)=1; T23(2,6)=1; T23(3,7)=1; T23(4,9)=1;

Nw12 = Nvx*T12;
Nw43 = Nvx*T43;
Nw14 = Nvy*T14;
Nw23 = Nvy*T23;

Ntx12 = Ntx1*T12;
Ntx43 = Ntx1*T43;
Nty14 = Nty1*T14;
Nty23 = Nty1*T23;

dNw12_dx = dNvx*T12;
dNw43_dx = dNvx*T43;
dNw14_dy = dNvy*T14;
dNw23_dy = dNvy*T23;

dNtx12_dx = dNtx1*T12;
dNtx43_dx = dNtx1*T43;
dNty14_dy = dNty1*T14;
dNty23_dy = dNty1*T23;

N1 = (1-rx)*(1-ry);
N2 = rx*(1-ry);
N3 = rx*ry;
N4 = (1-rx)*ry;

Nbil = zeros(1,12);
Nbil(1)=N1; Nbil(4)=N2; Nbil(7)=N3; Nbil(10)=N4;

dNbil_dx = zeros(1,12);
dNbil_dx(1)=-(1-ry)/Lx;
dNbil_dx(4)= (1-ry)/Lx;
dNbil_dx(7)= ry/Lx;
dNbil_dx(10)=-ry/Lx;

dNbil_dy = zeros(1,12);
dNbil_dy(1)=-(1-rx)/Ly;
dNbil_dy(4)=-rx/Ly;
dNbil_dy(7)= rx/Ly;
dNbil_dy(10)=(1-rx)/Ly;

Nw = (1-ry)*Nw12 + ry*Nw43 + (1-rx)*Nw14 + rx*Nw23 - Nbil;

dNw_dx = (1-ry)*dNw12_dx + ry*dNw43_dx ...
         -Nw14/Lx + Nw23/Lx - dNbil_dx;

dNw_dy = -Nw12/Ly + Nw43/Ly ...
         +(1-rx)*dNw14_dy + rx*dNw23_dy - dNbil_dy;

Ntx = (1-ry)*Ntx12 + ry*Ntx43;
dNtx_dx = (1-ry)*dNtx12_dx + ry*dNtx43_dx;
dNtx_dy = -Ntx12/Ly + Ntx43/Ly;

Nty = (1-rx)*Nty14 + rx*Nty23;
dNty_dx = -Nty14/Lx + Nty23/Lx;
dNty_dy = (1-rx)*dNty14_dy + rx*dNty23_dy;
end
