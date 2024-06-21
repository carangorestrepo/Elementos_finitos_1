clc
clear
syms w1 w2 w3 w4
sym S4 S5 S6 
syms C4 C5 C6
syms Bx1 Bx2 Bx3 Bx4
syms By1 By2 By3 By4
syms L4 L5 L6  
syms alfa4 alfa5 alfa6
b=1/5;

P44=(1/3+b); c11=(C4*C4+S4*S4)*alfa4; 
P45=b;       c12=(C4*C5+S4*S5)*alfa5; 
P46=b;       c13=(C4*C6+S4*S6)*alfa6; 

P54=b;       c21=(C5*C4+S5*S4)*alfa4; 
P55=(1/3+b); c22=(C5*C5+S5*S5)*alfa5;
P56=b;       c23=(C5*C6+S5*S6)*alfa6;

P64=b;       c31=(C6*C4+S6*S4)*alfa4; 
P65=b;       c32=(C6*C5+S6*S5)*alfa5; 
P66=(1/3+b); c33=(C6*C6+S6*S6)*alfa6;

gamas4=(w2-w1)/L4+1/2*C4*(Bx1+Bx2)+1/2*S4*(By1+By2);      
gamas5=(w3-w2)/L4+1/2*C5*(Bx2+Bx3)+1/2*S5*(By1+By2);   
gamas6=(w1-w3)/L4+1/2*C6*(Bx3+Bx1)+1/2*S6*(By1+By2);   
       
gamas4a=P44*c11+P45*c12+P46*c13;      
gamas5a=P54*c21+P55*c22+P56*c23;   
gamas6a=P64*c31+P65*c32+P66*c33;   

gamas4b=P44*c11+P45*c12+P46*c13;      
gamas5b=P54*c21+P55*c22+P56*c23;   
gamas6b=P64*c31+P65*c32+P66*c33;  


Un = [ w1 Bx1 By1 w2 Bx2 By2 w3 Bx3 By3].';
alfak = [ alfa4 alfa5 alfa6].';


An1 = simplify(equationsToMatrix([ gamas4a; gamas5a; gamas6a ], alfak));

A0=[[ (8*C4^2)/15 + (8*S4^2)/15,     (S4*S5)/5 + (C4*C5)/5,     (S4*S6)/5 + (C4*C6)/5];
    [     (S4*S5)/5 + (C4*C5)/5, (8*C5^2)/15 + (8*S5^2)/15,     (S5*S6)/5 + (C5*C6)/5];
    [     (S4*S6)/5 + (C4*C6)/5,     (S5*S6)/5 + (C5*C6)/5, (8*C6^2)/15 + (8*S6^2)/15]];







