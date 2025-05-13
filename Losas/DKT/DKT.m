%syms x23 x31 x12 y23 y31 y12 xi_gl eta_gl
clear
clc
x1=0;
x2=1;
x3=0;
y1=0;
y2=0;
y3=1;
x23 = x2-x3;
x31 = x3-x1;
x12 = x1-x2;
y23 = y2-y3;
y31 = y3-y1;
y12 = y1-y2;

Nforma = @(xi,eta) [2.0*(1.0-xi-eta)*(0.5-xi-eta)        % N1
                    xi*(2.0*xi-1.0)             % N2
                    eta*(2.0*eta-1.0)
                    4.0*xi*eta
                    4.0*eta*(1.0-xi-eta)
                    4.0*xi*(1.0-xi-eta)];      % N6
E=4*4880.0;
nu=1/3;
t=1;

D = E*t^3/(12*(1-nu^2)) * [
                           1.0  nu  0.0;
                           nu  1.0  0.0;
                           0.0 0.0  0.5*(1-nu)];
nef=1;   
x_gl = [1/2;1/2;0];
e_gl = [1/2;0;1/2];
w_gl =[1/3;1/3;1/3];
n_gl = size(x_gl,1);  %# Número de puntos de Gauss.
B  = cell(nef, n_gl, n_gl);
Kbe = zeros(9);%%

for e = 1:nef               % ciclo sobre todos los elementos finitos
    L23 = (x23^2 + y23^2);
    L31 = (x31^2 + y31^2);
    L12 = (x12^2 + y12^2);

    a4 = -x23/L23;
    a5 = -x31/L31;
    a6 = -x12/L12;
    b4 = 3*x23*y23/(4*L23);
    b5 = 3*x31*y31/(4*L31);
    b6 = 3*x12*y12/(4*L12);
    c4 = (x23^2 - 2*y23^2)/(4*L23);
    c5 = (x31^2 - 2*y31^2)/(4*L31);
    c6 = (x12^2 - 2*y12^2)/(4*L12);
    d4 = -y23/L23;
    d5 = -y31/L31;
    d6 = -y12/L12;
    e4 = (y23^2 - 2*x23^2)/(4*L23);
    e5 = (y31^2 - 2*x31^2)/(4*L31);
    e6 = (y12^2 - 2*x12^2)/(4*L12);

    p4 = -6*x23/L23;
    p5 = -6*x31/L31;
    p6 = -6*x12/L12;
    q4 =  3*x23*y23/L23;
    q5 =  3*x31*y31/L31;
    q6 =  3*x12*y12/L12;
    r4 =  3*y23^2/L23;
    r5 =  3*y31^2/L31;
    r6 =  3*y12^2/L12;
    t4 = -6*y23/L23;
    t5 = -6*y31/L31;
    t6 = -6*y12/L12;        
    for pp = 1:n_gl  
       for qq = 1:n_gl        
        xi_gl  = x_gl(pp);            
        eta_gl = e_gl(qq);
        [dHxdxi,dHxdeta,dHydxi,dHydeta]=dH(x12,x23,x31,y12,y23,y31,xi_gl,eta_gl);
                 B{e,pp,qq} =1.0/(x31*y12 - x12*y31) * [
                 y31*dHxdxi' + y12*dHxdeta';
                -x31*dHydeta' - x12*dHydeta';
                -x31*dHxdxi' - x12*dHxdeta' + y31*dHydxi' + y12*dHydeta'];

         Kbe = Kbe+B{e,pp,qq}'*D*B{e,pp}*w_gl(pp)*w_gl(qq);
       end
    end
a=1
end