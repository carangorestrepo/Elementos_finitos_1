clc
clear
syms C1 C2 C3 C4 C5 C6 EIz Acz GJx x L qZ qT
syms txi tyi wzi txf tyf wzf


Vz=int(qZ,x)+C1;
My=int(Vz,x)+C2;
ty=int(My/EIz,x)+C3;
wz=int(ty-Vz/Acz,x)+C4;% 

d2u=int(qT/GJx,x)+C5;
tx=int(d2u,x)+C6;

K_TE2 = sym(zeros(6));

[c1,c2,c3,c4,c5,c6]=solve(subs(tx,x,0)==txi,...
                          subs(ty,x,0)==tyi,...
                          subs(wz,x,0)==wzi,...
                          subs(tx,x,L)==txf,...
                          subs(ty,x,L)==tyf,...
                          subs(wz,x,L)==wzf,...
                          [C1,C2,C3,C4,C5,C6]);
V=simplify(subs(Vz,{C1,C2,C3,C4,C5,C6},{c1,c2,c3,c4,c5,c6}));
M=simplify(subs(My,{C1,C2,C3,C4,C5,C6},{c1,c2,c3,c4,c5,c6}));
w=simplify(subs(wz,{C1,C2,C3,C4,C5,C6},{c1,c2,c3,c4,c5,c6}));