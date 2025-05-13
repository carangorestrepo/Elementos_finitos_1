function [coor,connex,climit,charg,E,nu,th]=patcht4

E=100000;
nu=0.25;
th=0.1;

coor=[0.	0.;
0.24	0.;
0.24	0.12;
0.	0.12;
0.16	0.08];

connex=[1 2 5;5 2 3;4 5 3;1 5 4];

L=123456;
%climit=[1 0 0 L;2 2.4 1.2 L;3 3 2.4 L;4 0.6 1.2 L;5 2 1.6 L];
%charg=[];

climit=[1 0 0 L;4 0 L L];
charg=[2 1 0 0;3 1 0 0];
end


