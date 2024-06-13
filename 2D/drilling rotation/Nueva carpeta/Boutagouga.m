syms a1 a2 a3 a4 a5 a6 
syms b1 b2 b3 b4 b5 b6 
syms d1 d2 d3 d4 d5 d6

syms c12 c23 c31 
syms s12 s23 s31 
syms L12 L23 L31 

syms alpha_1 alpha_2 alpha_3 alpha_4 alpha_5 alpha_6 alpha_7 alpha_8 alpha_9

syms x y 

syms P4 P5 P6 
syms x1 x2 x3 
syms y1 y2 y3

%x1 = 0;
%y1 = 0;
%y2 = 0;
x4 = (x1+x2)/2;    x5 = (x2+x3)/2;    x6 = (x3+x1)/2;
y4 = (y1+y2)/2;    y5 = (y2+y3)/2;    y6 = (y3+y1)/2;

%# equations 14
eqs = [ a1 + a2*x1 + a3*y1 + a4*x1^2 + a5*x1*y1 + a6*y1^2 - 0;
        a1 + a2*x2 + a3*y2 + a4*x2^2 + a5*x2*y2 + a6*y2^2 - 0;
        a1 + a2*x3 + a3*y3 + a4*x3^2 + a5*x3*y3 + a6*y3^2 - 0;
        a1 + a2*x4 + a3*y4 + a4*x4^2 + a5*x4*y4 + a6*y4^2 - 1;
        a1 + a2*x5 + a3*y5 + a4*x5^2 + a5*x5*y5 + a6*y5^2 - 0;
        a1 + a2*x6 + a3*y6 + a4*x6^2 + a5*x6*y6 + a6*y6^2 - 0 ];
sol_a = solve(eqs, [a1, a2, a3, a4, a5, a6]);

%# ecuaciones 16,# print equations 17
eqs = [ b1 + b2*x1 + b3*y1 + b4*x1^2 + b5*x1*y1 + b6*y1^2 - 0;
        b1 + b2*x2 + b3*y2 + b4*x2^2 + b5*x2*y2 + b6*y2^2 - 0;
        b1 + b2*x3 + b3*y3 + b4*x3^2 + b5*x3*y3 + b6*y3^2 - 0;
        b1 + b2*x4 + b3*y4 + b4*x4^2 + b5*x4*y4 + b6*y4^2 - 0;
        b1 + b2*x5 + b3*y5 + b4*x5^2 + b5*x5*y5 + b6*y5^2 - 1;
        b1 + b2*x6 + b3*y6 + b4*x6^2 + b5*x6*y6 + b6*y6^2 - 0 ];
sol_b = solve(eqs, [b1, b2, b3, b4, b5, b6]);


%# ecuaciones 18
eqs = [ d1 + d2*x1 + d3*y1 + d4*x1^2 + d5*x1*y1 + d6*y1^2 - 0;
        d1 + d2*x2 + d3*y2 + d4*x2^2 + d5*x2*y2 + d6*y2^2 - 0;
        d1 + d2*x3 + d3*y3 + d4*x3^2 + d5*x3*y3 + d6*y3^2 - 0;
        d1 + d2*x4 + d3*y4 + d4*x4^2 + d5*x4*y4 + d6*y4^2 - 0;
        d1 + d2*x5 + d3*y5 + d4*x5^2 + d5*x5*y5 + d6*y5^2 - 0;
        d1 + d2*x6 + d3*y6 + d4*x6^2 + d5*x6*y6 + d6*y6^2 - 1 ];
sol_d = solve(eqs, [d1, d2, d3, d4, d5, d6]);

%# equations 13
dp4 = (a1 + a2*x + a3*y + a4*x^2 + a5*x*y + a6*y^2)*P4;
dp5 = (b1 + b2*x + b3*y + b4*x^2 + b5*x*y + b6*y^2)*P5;
dp6 = (d1 + d2*x + d3*y + d4*x^2 + d5*x*y + d6*y^2)*P6;

%# equations 20
up4 = dp4*c12;    up5 = dp5*c23;    up6 = dp6*c31;
vp4 = dp4*s12;    vp5 = dp5*s23;    vp6 = dp6*s31;

%# equations 10
ul = alpha_1 + alpha_2*x + alpha_3*y;
vl = alpha_4 + alpha_5*x + alpha_6*y;

%# equations 12
up = dp4*c12 + dp5*c23 + dp6*c31;
vp = dp4*s12 + dp5*s23 + dp6*s31;

%% # P4 = alpha_7,   P5 = alpha_8,   P6 = alpha_9
%# equations 9
u = subs(ul + up,{P4     ,P5     ,P6     ,a1      ,b1      ,b2      ,b3      ,b4      ,d1      ,d2      ,d4},...
                 {alpha_7,alpha_8,alpha_9,sol_a.a1,sol_b.b1,sol_b.b2,sol_b.b3,sol_b.b4,sol_d.d1,sol_d.d2,sol_d.d4});
v = subs(vl + vp,{P4     ,P5     ,P6     ,a1      ,b1      ,b2      ,b3      ,b4      ,d1      ,d2      ,d4},...
                 {alpha_7,alpha_8,alpha_9,sol_a.a1,sol_b.b1,sol_b.b2,sol_b.b3,sol_b.b4,sol_d.d1,sol_d.d2,sol_d.d4});
%# equation 7
phi = subs((diff(v, x) - diff(u, y))/2,{P4,P5,P6},{alpha_7,alpha_8,alpha_9});

%# equation 26
ex  = diff(u, x);
ey  = diff(v, y);
gxy = diff(u, y) + diff(v, x);
% # equation 27
alphas = [ alpha_1 alpha_2 alpha_3 alpha_4 alpha_5 alpha_6 alpha_7 alpha_8 alpha_9].';
B = equationsToMatrix([ ex; ey; gxy], alphas);

Axy=equationsToMatrix([ u; v; phi], alphas);



c12i =  (y2 - y1)/L12;    c23i =  (y3 - y2)/L23;    c31i =  (y1 - y3)/L31;
s12i = -(x2 - x1)/L12;    s23i = -(x3 - x2)/L23;    s31i = -(x1 - x3)/L31;

A=[subs(Axy  ,{x,y,c12,c23,c31,s12,s23,s31},{x1,y1,c12i,c23i,c31i,s12i,s23i,s31i});
   subs(Axy  ,{x,y,c12,c23,c31,s12,s23,s31},{x2,y2,c12i,c23i,c31i,s12i,s23i,s31i});
   subs(Axy  ,{x,y,c12,c23,c31,s12,s23,s31},{x3,y3,c12i,c23i,c31i,s12i,s23i,s31i})];
A=simplify(subs(A,{       a2,       a3,       a4,       a5,       a6,       b5,       b6,       d3,       d5,       d6},...
                  { sol_a.a2, sol_a.a3, sol_a.a4, sol_a.a5, sol_a.a6, sol_b.b5, sol_b.b6, sol_d.d3, sol_d.d5, sol_d.d6}));

% # equation 32
phi1 = A(3,:);
phi2 = A(6,:);
phi3 = A(9,:);

phiR=(phi1+phi2+phi3)/3;

% # equation 33
n1=phi2 - phi1;
n2=phi3 - phi2;
n3=phi1 - phi3;

%Estimate the rotation due to linear displacement:
phi0=simplify(subs(phi,{x               ,y               ,c12 ,c23 ,c31 ,s12 ,s23 ,s31,       a2,       a3,       a4,       a5,       a6,       b5,       b6,       d3,       d5,       d6},...
                       {(x1 + x2 + x3)/3,(y1 + y2 + y3)/3,c12i,c23i,c31i,s12i,s23i,s31i,sol_a.a2, sol_a.a3, sol_a.a4, sol_a.a5, sol_a.a6, sol_b.b5, sol_b.b6, sol_d.d3, sol_d.d5, sol_d.d6}));
%# equation 44
D = 1/L12;
E = 1/L23;
F = 1/L31;
syms k
%k=inf                  
A_amend = [[0, 0, 0, 0, 0, 0,   0,   0,   0];
           [0, 0, 0, 0, 0, 0,   0,   0,   0];
           [0, 0, 0, 0, 0, 0, k*D, k*E, k*F];
           [0, 0, 0, 0, 0, 0,   0,   0,   0];
           [0, 0, 0, 0, 0, 0,   0,   0,   0];
           [0, 0, 0, 0, 0, 0, k*D, k*E, k*F];
           [0, 0, 0, 0, 0, 0,   0,   0,   0];
           [0, 0, 0, 0, 0, 0,   0,   0,   0];
           [0, 0, 0, 0, 0, 0, k*D, k*E, k*F]];

%# equation 46
A2 = A + A_amend;
%# equation 47
Ainv = inv(A2);
%# equation 48
H = limit(Ainv,k,inf);


