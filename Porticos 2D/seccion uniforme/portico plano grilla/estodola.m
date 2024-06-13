function [Phi,Omega,T]=estodola(k,m,iteraciones,modos)
sizek=size(k,1);
Bo=diag(ones(1,sizek));
x1n  = cell(modos,1);
xm1n  = cell(modos,1);
Phi=zeros(sizek,modos);
Omega=zeros(modos,1);
for n=1:modos
    Fb=k^-1*m*Bo;
    x1=ones(sizek,iteraciones-1);
    xm1=zeros(sizek,iteraciones-1);
    for i=2:iteraciones
        xm1(:,i-1)=Fb*x1(:,i-1);
        x1(:,i)=xm1(:,i-1)/max(abs(xm1(:,i-1)));
    end
    x1n{n}=x1;
    xm1n {n}=xm1;
    Bo=Bo-x1(:,i)*x1(:,i)'*m/(x1(:,i)'*m*x1(:,i));
    Phi(:,n)=x1(:,i-1);
    Omega(n)=sqrt(abs(1/max(abs(xm1(:,i-1)))));
end
[Omega,In] = sort(Omega);     % ordena las frecuencias
T = 2*pi./Omega;             % seg - periodo de vibracion 
Phi = Phi(:,In); % ordena los modos segun el orden de las frecuencias
Phi = Phi./repmat(sqrt(diag(Phi'*m*Phi)'),size(Phi,1),1); % normaliza los modos