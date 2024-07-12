clc
clear
syms C1 C2 C3 C4 C5 C6 EIz Acz GJx x L
qz=0;
q=0;
%L=1.2;
%EIz=30400.0000000000;
%GJx=15388.8171386719;
%Acz=791666.6666666670;

Vz=int(qz,x)+C1;
My=int(Vz,x)+C2;
ty=int(My/EIz,x)+C3;
wz=int(ty-Vz/Acz,x)+C4;% 

d2u=int(q/GJx,x)+C5;
tx=int(d2u,x)+C6;

K_TE2 = sym(zeros(6));
for i = 1:6
[c1,c2,c3,c4,c5,c6]=solve(subs(tx,x,0)==(i==1),...
                          subs(ty,x,0)==(i==2),...
                          subs(wz,x,0)==(i==3),...
                          subs(tx,x,L)==(i==4),...
                          subs(ty,x,L)==(i==5),...
                          subs(wz,x,L)==(i==6),...
                          [C1,C2,C3,C4,C5,C6]);                 
   % # se evaluan las reacciones horizontales y verticales y los momentos en los apoyos
            K_TE2(:,i)=[ -subs(d2u*GJx ,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,0}); % X1
                         -subs(My      ,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,0}); % Y2
                         -subs(Vz      ,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,0}); % M2
                          subs(d2u*GJx ,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L}); % X2 
                          subs(My      ,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L}); % Y2
                          subs(Vz      ,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L})] ;% M2
end
%K_TE2=double(K_TE2);                   
K_TE2(2,3)=(-K_TE2(2,3)); 
K_TE2(3,3)=(-K_TE2(3,3)); 
K_TE2(5,3)=(-K_TE2(5,3)); 
K_TE2(6,3)=(-K_TE2(6,3)); 

K_TE2(2,6)=(-K_TE2(2,6)); 
K_TE2(3,6)=(-K_TE2(3,6)); 

K_TE2(5,6)=(-K_TE2(5,6));
K_TE2(6,6)=(-K_TE2(6,6));

K_TE2=simplify(K_TE2);
