function [coor,connex,climit,charg,E,nu,th]=torquetest

E=1000000;
nu=0.25;
th=0.01;

coor=[0.	0.;
1.	0.;
0.	1.;
-1.	0.;
0.	-1.];

connex=[1 2 3;4 1 3;4 5 1;5 2 1];

L=123456;

%climit=[1 0 0 0.1];
%climit=[1 0 0 L;2 L 0.1 L];
%charg=[];

%climit=[1 0 0 0];
%charg=[2 0 5 0;
%4 0 -5 0];

climit=[1 0 0 L;4 L 0 L];
charg=[1 0 0 5];


end