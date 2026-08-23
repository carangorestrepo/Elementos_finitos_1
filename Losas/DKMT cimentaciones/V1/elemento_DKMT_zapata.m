function [Kb,Ks,Kw,fe,fe_comp,Bb_c,Bs_c,Ae] = ...
    elemento_DKMT_zapata(xe,E,nu,h,kappa,ks,cargas,en_huella,orden_gauss)
% Matrices del elemento triangular DKMT y carga nodal consistente.
% Orden: [w1 bx1 by1 w2 bx2 by2 w3 bx3 by3]

if size(xe,1)~=3 || size(xe,2)~=2
    error('xe debe ser una matriz 3x2.');
end

x=xe(:,1);
y=xe(:,2);

J=[x(2)-x(1),y(2)-y(1);
   x(3)-x(1),y(3)-y(1)];

detJ=det(J);

if detJ<=100*eps(max(1,norm(J,'fro')))
    error('Elemento degenerado o con orientacion no positiva.');
end

invJ=inv(J);
Ae=detJ/2;

xji=[x(2)-x(1),x(3)-x(2),x(1)-x(3)];
yji=[y(2)-y(1),y(3)-y(2),y(1)-y(3)];

Lk=hypot(xji,yji);
Ck=xji./Lk;
Sk=yji./Lk;

phi=2/(kappa*(1-nu))*(h./Lk).^2;

A_dbeta=diag((2/3)*Lk.*(1+phi));

Aw=[ 1,-xji(1)/2,-yji(1)/2,-1,-xji(1)/2,-yji(1)/2,0,0,0;
     0,0,0,1,-xji(2)/2,-yji(2)/2,-1,-xji(2)/2,-yji(2)/2;
    -1,-xji(3)/2,-yji(3)/2,0,0,0,1,-xji(3)/2,-yji(3)/2];

An=A_dbeta\Aw;

dN_nat=[-1,1,0;
        -1,0,1];

dNxy=invJ*dN_nat;
dNdx=dNxy(1,:);
dNdy=dNxy(2,:);

Bb_beta=zeros(3,9);

for i=1:3
    c=3*(i-1);
    Bb_beta(:,c+(1:3))=[0,dNdx(i),0;
                        0,0,dNdy(i);
                        0,dNdy(i),dNdx(i)];
end

C4=Ck(1); C5=Ck(2); C6=Ck(3);
S4=Sk(1); S5=Sk(2); S6=Sk(3);

A1=C4*S6-C6*S4;
A2=C5*S4-C4*S5;
A3=C6*S5-C5*S6;

if any(abs([A1,A2,A3])<1e-12)
    error('Geometria inadecuada para la matriz de cortante DKMT.');
end

M5=[ S6/A1,       0,-S4/A1;
    -C6/A1,       0, C4/A1;
    -S5/A2,   S4/A2,      0;
     C5/A2,  -C4/A2,      0;
          0, -S6/A3, S5/A3;
          0,  C6/A3,-C5/A3];

M8=diag(-(2/3)*phi);

Hb=E*h^3/(12*(1-nu^2))* ...
   [1,nu,0;
    nu,1,0;
    0,0,(1-nu)/2];

G=E/(2*(1+nu));
Hs=kappa*G*h*eye(2);

gp=cuadratura_triangular(orden_gauss);

Kb=zeros(9);
Ks=zeros(9);
Kw=zeros(9);
fe=zeros(9,1);

fe_comp.PM=zeros(9,1);
fe_comp.zapata=zeros(9,1);
fe_comp.suelo=zeros(9,1);
fe_comp.pedestal=zeros(9,1);

for p=1:size(gp,1)
    xi=gp(p,1);
    eta=gp(p,2);
    wg=gp(p,3);

    N=[1-xi-eta,xi,eta];
    Nw=[N(1),0,0,N(2),0,0,N(3),0,0];

    dP_dxi=[4-8*xi-4*eta;
            4*eta;
           -4*eta];

    dP_deta=[-4*xi;
              4*xi;
              4-4*xi-8*eta];

    dPxy=invJ*[dP_dxi.';dP_deta.'];
    dPdx=dPxy(1,:);
    dPdy=dPxy(2,:);

    Bb_dbeta=zeros(3,3);

    for k=1:3
        Bb_dbeta(:,k)=[dPdx(k)*Ck(k);
                       dPdy(k)*Sk(k);
                       dPdy(k)*Ck(k)+dPdx(k)*Sk(k)];
    end

    Bb=Bb_beta+Bb_dbeta*An;

    M6=[N(1),0,N(2),0,N(3),0;
        0,N(1),0,N(2),0,N(3)];

    Bs=M6*M5*M8*An;

    dA=detJ*wg;

    Kb=Kb+Bb.'*Hb*Bb*dA;
    Ks=Ks+Bs.'*Hs*Bs*dA;

    % Winkler solo actua sobre w, sin factor h
    Kw=Kw+Nw.'*ks*Nw*dA;

    xgp=N*xe(:,1);
    ygp=N*xe(:,2);

    [qtotal,qc]=evaluar_carga_zapata(cargas,xgp,ygp,en_huella);

    fe=fe+Nw.'*qtotal*dA;

    fe_comp.PM=fe_comp.PM+Nw.'*qc.qPM*dA;
    fe_comp.zapata=fe_comp.zapata+Nw.'*qc.q_zapata*dA;
    fe_comp.suelo=fe_comp.suelo+Nw.'*qc.q_suelo*dA;
    fe_comp.pedestal=fe_comp.pedestal+Nw.'*qc.q_pedestal*dA;
end

Kb=0.5*(Kb+Kb.');
Ks=0.5*(Ks+Ks.');
Kw=0.5*(Kw+Kw.');

% Operadores en el centroide para recuperacion
xi=1/3;
eta=1/3;
N=[1/3,1/3,1/3];

dP_dxi=[4-8*xi-4*eta;
        4*eta;
       -4*eta];

dP_deta=[-4*xi;
          4*xi;
          4-4*xi-8*eta];

dPxy=invJ*[dP_dxi.';dP_deta.'];
dPdx=dPxy(1,:);
dPdy=dPxy(2,:);

Bb_dbeta=zeros(3,3);
for k=1:3
    Bb_dbeta(:,k)=[dPdx(k)*Ck(k);
                   dPdy(k)*Sk(k);
                   dPdy(k)*Ck(k)+dPdx(k)*Sk(k)];
end

Bb_c=Bb_beta+Bb_dbeta*An;

M6=[N(1),0,N(2),0,N(3),0;
    0,N(1),0,N(2),0,N(3)];

Bs_c=M6*M5*M8*An;

end

function gp=cuadratura_triangular(orden)
% Pesos sobre triangulo natural de area 1/2.

if orden<=1
    gp=[1/3,1/3,1/2];

elseif orden<=2
    gp=[1/6,1/6,1/6;
        2/3,1/6,1/6;
        1/6,2/3,1/6];

else
    gp=[0.445948490915965,0.445948490915965,0.111690794839005;
        0.445948490915965,0.108103018168070,0.111690794839005;
        0.108103018168070,0.445948490915965,0.111690794839005;
        0.091576213509771,0.091576213509771,0.054975871827661;
        0.091576213509771,0.816847572980459,0.054975871827661;
        0.816847572980459,0.091576213509771,0.054975871827661];
end
end
