function [coor,connex,climit,charg,E,nu,th]=pureshear

E=1000000;
nu=0.25;
th=0.001;

coor=[0.	0.;
0.2	0.;
0.2	0.2;
0.	0.2;
0.1	0.1];

connex=[1 2 5;5 2 3;4 5 3;1 5 4];

L=123456;
climit=[1 0 0 0];

charg=[2 1 -1 0;3 -1 -1 0;4 -1 1 0];

end