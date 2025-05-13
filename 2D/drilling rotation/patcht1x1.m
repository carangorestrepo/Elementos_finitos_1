function [coor,connex,climit,charg,E,nu,th]=patcht1x1

E=100000;
nu=0.25;
th=0.1;

coor=[0.	0.;
      0.24	0.;
      0.24	0.12;
      0.	0.12];

connex=[1 2 4;4 2 3];

L=123456;
climit=[1 0 0 L;4 0 L L];

charg=[2 1 0 0;
3 1 0 0];

end
