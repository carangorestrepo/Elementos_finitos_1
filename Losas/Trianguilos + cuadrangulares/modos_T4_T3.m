%% Solución del problema de valores propios
OPT.maxit = 10000;
OPT.tol   = 1e-32;
OPT.issym = 1;  % La matriz Kdd es simetrica
%numero de modos
[Phi,Omega2,flag] = eigs(K(d,d),M(d,d),10,'SM',OPT);% autovalores y autovectores

Omega = sqrt(diag(Omega2));  % frecuencias angulares
[Omega,I] = sort(Omega);     % ordena las frecuencias
Te = 2*pi./Omega;             % seg - periodo de vibracion 
Phi = Phi(:,I); % ordena los modos segun el orden de las frecuencias
Phi = Phi./repmat(sqrt(diag(Phi'*M(d,d)*Phi)'),size(Phi,1),1); % normaliza los modos

modo=6;
%% Visualización de modos de vibración

def = zeros(ngdl,1);
def(d,1) = Phi(:,modo);
vect_mov = reshape(def,3,nno)';

figure;
hold on;
for e = 1:nef    
    if Numero_de_nodos_elem(e)==3 %% elementos triangulares
       fill3(xnod(LaG(e,[1 2 3 1]),X), ...
             xnod(LaG(e,[1 2 3 1]),Y), ...
             vect_mov(LaG(e,[1 2 3 1]),ww), ...
             vect_mov(LaG(e,[1 2 3 1]),ww));
   elseif Numero_de_nodos_elem(e)==4 %  elemetos cudrangulares 
      fill3(xnod(LaG(e,[1 2 3 4 1]),X), ...
            xnod(LaG(e,[1 2 3 4 1]),Y), ...
            vect_mov(LaG(e,[1 2 3 4 1]),ww), ...
           vect_mov(LaG(e,[1 2 3 4 1]),ww));
   end
end
nombre=double('modo');
nombre2=double('W (rad/sec)=');
modo1=num2str(modo);
val=num2str(Omega(modo));
nombre=[nombre,32,modo1,32,nombre2,val];
title(char(nombre));
colormap jet; colorbar;
view(3); axis tight;