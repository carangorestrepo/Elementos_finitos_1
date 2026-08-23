function resultado = resolver_zapata_DKMT_modular( ...
    xnod,LaG,in_col,geo,mat,suelo,cargas,opciones)
% Ensambla y resuelve (Kplaca + Ksuelo)*a = f.

nno=size(xnod,1);
nef=size(LaG,1);
ngdl=3*nno;

gdl=[(1:3:ngdl).',(2:3:ngdl).',(3:3:ngdl).'];

Kplaca=sparse(ngdl,ngdl);
Ksuelo=sparse(ngdl,ngdl);

f_total=zeros(ngdl,1);
f_PM=zeros(ngdl,1);
f_zapata=zeros(ngdl,1);
f_suelo=zeros(ngdl,1);
f_pedestal=zeros(ngdl,1);

idx=cell(nef,1);
Bbc=zeros(3,9,nef);
Bsc=zeros(2,9,nef);
Ae=zeros(nef,1);

for e=1:nef
    nod=LaG(e,:);
    xe=xnod(nod,:);

    [Kbe,Kse,Kwe,fe,fe_comp,Bb_cent,Bs_cent,Aelem] = ...
        elemento_DKMT_zapata( ...
        xe,mat.E,mat.nu,geo.hz,mat.kappa,suelo.ks, ...
        cargas,in_col(e),opciones.orden_gauss);

    idx{e}=[gdl(nod(1),:),gdl(nod(2),:),gdl(nod(3),:)];

    ie=idx{e};

    Kplaca(ie,ie)=Kplaca(ie,ie)+Kbe+Kse;
    Ksuelo(ie,ie)=Ksuelo(ie,ie)+Kwe;

    f_total(ie)=f_total(ie)+fe;
    f_PM(ie)=f_PM(ie)+fe_comp.PM;
    f_zapata(ie)=f_zapata(ie)+fe_comp.zapata;
    f_suelo(ie)=f_suelo(ie)+fe_comp.suelo;
    f_pedestal(ie)=f_pedestal(ie)+fe_comp.pedestal;

    Bbc(:,:,e)=Bb_cent;
    Bsc(:,:,e)=Bs_cent;
    Ae(e)=Aelem;
end

K=Kplaca+Ksuelo;
K=0.5*(K+K.');

if isfield(opciones,'gdl_restringidos')
    c=opciones.gdl_restringidos(:);
else
    c=[];
end

d=setdiff((1:ngdl).',c);
a=zeros(ngdl,1);
q_reaccion=zeros(ngdl,1);

if isempty(c)
    a=K\f_total;
else
    ac=zeros(length(c),1);

    Kcc=K(c,c);
    Kcd=K(c,d);
    Kdc=K(d,c);
    Kdd=K(d,d);

    fc=f_total(d);
    fd=f_total(c);

    ad=Kdd\(fc-Kdc*ac);
    qd=Kcc*ac+Kcd*ad-fd;

    a(c)=ac;
    a(d)=ad;
    q_reaccion(c)=qd;
end

Rsuelo=-Ksuelo*a;

w=a(1:3:end);
betax=a(2:3:end);
betay=a(3:3:end);

Fw=f_total(1:3:end);
Rw=Rsuelo(1:3:end);

equilibrio.Fz=sum(Fw);
equilibrio.Rz=sum(Rw);

equilibrio.Mx_ext=sum(Fw.*(xnod(:,2)-cargas.yc));
equilibrio.My_ext=sum(Fw.*(xnod(:,1)-cargas.xc));

equilibrio.Mx_suelo=sum(Rw.*(xnod(:,2)-cargas.yc));
equilibrio.My_suelo=sum(Rw.*(xnod(:,1)-cargas.xc));

Hb=mat.E*geo.hz^3/(12*(1-mat.nu^2))* ...
   [1,mat.nu,0;
    mat.nu,1,0;
    0,0,(1-mat.nu)/2];

G=mat.E/(2*(1+mat.nu));
Hs=mat.kappa*G*geo.hz*eye(2);

momentos=zeros(nef,3);
cortantes=zeros(nef,2);

for e=1:nef
    ae=a(idx{e});
    momentos(e,:)=(Hb*Bbc(:,:,e)*ae).';
    cortantes(e,:)=(Hs*Bsc(:,:,e)*ae).';
end

[Mx_n,My_n,Mxy_n,Qx_n,Qy_n] = ...
    promedio_nodal_resultados(LaG,Ae,momentos,cortantes,nno);

resultado.K=K;
resultado.Kplaca=Kplaca;
resultado.Ksuelo=Ksuelo;

resultado.f=f_total;
resultado.f_PM=f_PM;
resultado.f_zapata=f_zapata;
resultado.f_suelo=f_suelo;
resultado.f_pedestal=f_pedestal;

resultado.a=a;
resultado.q_reaccion=q_reaccion;
resultado.Rsuelo=Rsuelo;

resultado.w=w;
resultado.betax=betax;
resultado.betay=betay;
resultado.presion_nodal=-suelo.ks*w;

resultado.equilibrio=equilibrio;
resultado.momentos_elemento=momentos;
resultado.cortantes_elemento=cortantes;

resultado.Mx=Mx_n;
resultado.My=My_n;
resultado.Mxy=Mxy_n;
resultado.Qx=Qx_n;
resultado.Qy=Qy_n;

resultado.idx=idx;
resultado.area_elemento=Ae;

resultado.resultantes.PM=sum(f_PM(1:3:end));
resultado.resultantes.zapata=sum(f_zapata(1:3:end));
resultado.resultantes.suelo=sum(f_suelo(1:3:end));
resultado.resultantes.pedestal=sum(f_pedestal(1:3:end));
resultado.resultantes.total=sum(f_total(1:3:end));

fprintf('\nRESULTANTES INTEGRADAS POR EL MEF\n');
fprintf('P,Mx,My: carga vertical = %12.6f kN\n',resultado.resultantes.PM);
fprintf('Peso zapata            = %12.6f kN\n',resultado.resultantes.zapata);
fprintf('Peso suelo             = %12.6f kN\n',resultado.resultantes.suelo);
fprintf('Peso pedestal          = %12.6f kN\n',resultado.resultantes.pedestal);
fprintf('Total externo          = %12.6f kN\n',resultado.resultantes.total);

fprintf('\nEQUILIBRIO GLOBAL\n');
fprintf('F externa + R suelo = %12.6e kN\n', ...
    equilibrio.Fz+equilibrio.Rz);
fprintf('Mx externo + Mx suelo = %12.6e kN*m\n', ...
    equilibrio.Mx_ext+equilibrio.Mx_suelo);
fprintf('My externo + My suelo = %12.6e kN*m\n', ...
    equilibrio.My_ext+equilibrio.My_suelo);

end

function [Mx_n,My_n,Mxy_n,Qx_n,Qy_n] = ...
    promedio_nodal_resultados(LaG,Ae,momentos,cortantes,nno)

Mx_n=zeros(nno,1);
My_n=zeros(nno,1);
Mxy_n=zeros(nno,1);
Qx_n=zeros(nno,1);
Qy_n=zeros(nno,1);
peso=zeros(nno,1);

for e=1:size(LaG,1)
    nod=LaG(e,:);
    pe=Ae(e);

    Mx_n(nod)=Mx_n(nod)+pe*momentos(e,1);
    My_n(nod)=My_n(nod)+pe*momentos(e,2);
    Mxy_n(nod)=Mxy_n(nod)+pe*momentos(e,3);
    Qx_n(nod)=Qx_n(nod)+pe*cortantes(e,1);
    Qy_n(nod)=Qy_n(nod)+pe*cortantes(e,2);
    peso(nod)=peso(nod)+pe;
end

Mx_n=Mx_n./peso;
My_n=My_n./peso;
Mxy_n=Mxy_n./peso;
Qx_n=Qx_n./peso;
Qy_n=Qy_n./peso;
end
