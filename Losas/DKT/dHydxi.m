function [dd]=dHydxi(t4,t5,t6,r4,r5,r6,q4,q5,q6,xi,eta)

d1 =  t6*(1.0-2.0*xi) + eta*(t5-t6);
d2 =  1.0 + r6*(1.0-2.0*xi) - eta*(r5+r6);
d3 = -q6*(1.0-2.0*xi) + eta*(q5+q6);
d4 = -t6*(1.0-2.0*xi) + eta*(t4+t6);
d5 = -1.0 + r6*(1.0-2.0*xi) + eta*(r4-r6);
d6 = -q6*(1.0-2*xi) - eta*(q4-q6);
d7 = -eta*(t4+t5);
d8 =  eta*(r4-r5);
d9 = -eta*(q4-q5);
dd=[d1,d2,d3,d4,d5,d6,d7,d8,d9];