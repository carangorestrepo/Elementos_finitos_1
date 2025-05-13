function [da]=dHxdeta(p4,p5,p6,r4,r5,r6,q4,q5,q6,xi,eta)
d1 = -p5*(1.0-2.0*eta) - xi*(p6-p5);
d2 =  q5*(1.0-2.0*eta) - xi*(q5+q6);
d3 = -4.0 + 6.0*(xi+eta) + r5*(1.0-2.0*eta) - xi*(r5+r6);
d4 =  xi*(p4+p6);
d5 =  xi*(q4-q6);
d6 = -xi*(r6-r4);
d7 =  p5*(1.0-2.0*eta) - xi*(p4+p5);
d8 =  q5*(1.0-2.0*eta) + xi*(q4-q5);
d9 = -2.0 + 6.0*eta + r5*(1.0-2.0*eta+xi*(r4-r5));

da=[d1,d2,d3,d4,d5,d6,d7,d8,d9];