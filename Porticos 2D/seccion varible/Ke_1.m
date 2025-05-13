function ke=Ke_1(b1,b2,nb,h1,h2,nh,E,G,L,tw,tf,sec)

[k2]=Ax_Ix_As2(b1,b2,nb,h1,h2,nh,0,0,0,0,0,0,0,L,E,G,L,tw,tf,sec,0,1);
[k3]=Ax_Ix_As2(b1,b2,nb,h1,h2,nh,0,0,0,0,0,0,0,L,E,G,L,tw,tf,sec,-1,0);
[k5]=Ax_Ix_As2(b1,b2,nb,h1,h2,nh,0,0,0,0,0,0,-L,0,E,G,L,tw,tf,sec,0,-1);
[k6]=Ax_Ix_As2(b1,b2,nb,h1,h2,nh,0,0,0,0,0,0,-L,0,E,G,L,tw,tf,sec,-1,0);
k11=[k2(1);0;0;k2(4);0;0];
k22=[0;k2(2);k2(3);0;k2(5);k2(6)];
k33=[0;k3(2);k3(3);0;k3(5);k3(6)];
k44=[-k2(1);0;0;-k2(4);0;0];
k55=[0;-k5(5);k5(6);0;-k5(2);k5(3)];
k66=[0;-k6(5);k6(6);0;-k6(2);k6(3)];
ke=[k11,k22,k33,k44,k55,k66];