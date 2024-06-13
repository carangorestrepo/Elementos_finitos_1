
syms xi eta

lamda=1-xi-eta

Nforma=[lamda;
        xi;
        eta]
dN_dxi=diff(Nforma,xi)
dN_deta=diff(Nforma,eta)
Pforma=[4*lamda*xi;
    4*xi*eta;
    4*lamda*eta];
dP_dxi=diff(Pforma,xi)
dP_deta=diff(Pforma,eta)
