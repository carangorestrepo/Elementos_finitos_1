function ke=Ke(nq,sec,b1,b2,nb,h1,h2,nh,G,E,tf,tw,L)

[k2]=funcion(0,0,nq,sec,b1,b2,nb,h1,h2,nh,G,E,tf,tw,L,0,1,0,L);
[k3]=funcion(0,0,nq,sec,b1,b2,nb,h1,h2,nh,G,E,tf,tw,L,-1,0,0,L);
[k5]=funcion(0,0,nq,sec,b1,b2,nb,h1,h2,nh,G,E,tf,tw,L,0,-1,-L,0);
[k6]=funcion(0,0,nq,sec,b1,b2,nb,h1,h2,nh,G,E,tf,tw,L,-1,0,-L,0);
[k1]=funcion1(tf,tw,b1,b2,nb,h1,h2,nh,0,0,nq,L,E,G,1,0,L,sec);
k11=[k1(1);0;0;k1(2);0;0];
k22=[0;k2(1);k2(2);0;k2(3);k2(4)];
k33=[0;k3(1);k3(2);0;k3(3);k3(4)];
k44=-[k1(1);0;0;k1(2);0;0];
k55=[0;-k5(3);k5(4);0;-k5(1);k5(2)];
k66=[0;-k6(3);k6(4);0;-k6(1);k6(2)];
ke=[k11,k22,k33,k44,k55,k66];