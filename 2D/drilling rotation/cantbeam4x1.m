function [coor,connex,climit,charg,E,nu,th]=cantbeam4x1

E=30000;
nu=0.25;
th=1;

coor=[0.	-6.;
12.	-6.;
24.	-6.;
36.	-6.;
48.	-6.;
0.	6.;
12.	6.;
24.	6.;
36.	6.;
48.	6.];

connex=[1 2 6;6 2 7;2 3 7;7 3 8;3 4 8;8 4 9;4 5 9;9 5 10];

L=123456;
climit=[1 0 0 0;6 0 0 0];

charg=[5 0 20 0;
10 0 20 0];

end