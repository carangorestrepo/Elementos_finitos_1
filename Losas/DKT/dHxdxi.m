function [db]=dHxdxi(p4,p5,p6,r4,r5,r6,q4,q5,q6,xi,eta)

    d1 =  p6*(1.0-2.0*xi) + (p5-p6)*eta;
    d2 =  q6*(1.0-2.0*xi) - (q5+q6)*eta;
    d3 = -4.0 + 6.0*(xi+eta) + r6*(1.0-2.0*xi) - eta*(r5+r6);
    d4 = -p6*(1.0-2.0*xi) + eta*(p4+p6);
    d5 =  q6*(1.0-2.0*xi) - eta*(q6-q4);
    d6 = -2.0 + 6.0*xi + r6*(1.0-2.0*xi) + eta*(r4-r6);
    d7 = -eta*(p5+p4);
    d8 =  eta*(q4-q5);
    d9 = -eta*(r5-r4);
db=[d1,d2,d3,d4,d5,d6,d7,d8,d9];