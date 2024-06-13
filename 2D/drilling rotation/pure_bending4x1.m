function [coor,connex,climit,charg,E,nu,th]=pure_bending4x1

E=1500;
nu=0.;
th=1;

coor=[0.	-0.5;
2.5	-0.5;
5.	-0.5;
7.5	-0.5;
10.	-0.5;
0.	0.5;
2.5	0.5;
5.	0.5;
7.8	0.5;
10.	0.5];
coor(:,2)=coor(:,2)*2;

connex=[1 2 6;6 2 7;2 3 7;7 3 8;3 4 8;8 4 9;4 5 9;9 5 10];

L=123456;
climit=[1 0 0 L;6 0 0 L];

charg=[5 10 0 0;
10 -10 0 0];
%charg=[5 0 0 10;
%10 0 0 10];

end