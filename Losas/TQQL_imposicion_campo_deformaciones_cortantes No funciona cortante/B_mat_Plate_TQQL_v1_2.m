function [bmat_b,bmat_s,Ngp,area] = B_mat_Plate_TQQL_v1_2(x,y,xgs,ygs)

%% B_Mat Computes the strain-displacement matrix for bending moments
%        and shear forces
%
%  Parameters:
%
%    Input, x   : X Coordinates of the element
%           y   : Y Coordinates of the element
%           xgs : Local X coordinate of the Gauss point
%           ygs : Local Y coordinate of the Gauss point
%   
%    Output, bmat_b the strain-displacement matrix for bending moments
%            bmat_s the strain-displacement matrix for shear forces
%               Ngp the local shape functions
%              area the element area

%==================

  L1=1.0-xgs-ygs;
  L2=xgs;
  L3=ygs;

  Ngp(1) = (2*L1-1)*L1;    % OK
  Ngp(2) = 4*L1*L2;        % OK
  Ngp(3) = (2*L2-1)*L2;    % OK
  Ngp(4) = 4*L2*L3;        % OK
  Ngp(5) = (2*L3-1)*L3;    % OK
  Ngp(6) = 4*L1*L3;        % OK

  dxNloc(1)=1.0-4.0*L1;  % OK
  dxNloc(2)=4.0*(L1-L2); % OK
  dxNloc(3)=4.0*L2-1.0;  % OK
  dxNloc(4)=4.0*L3;      % OK
  dxNloc(5)=0.0;         % OK
  dxNloc(6)=-4.0*L3;     % OK

  dyNloc(1)=1.0-4.0*L1;  % OK
  dyNloc(2)=-4.0*L2;     % OK
  dyNloc(3)=0.0;         % OK
  dyNloc(4)=4.0*L2;      % OK
  dyNloc(5)=4.0*L3-1.0;  % OK
  dyNloc(6)=4.0*(L1-L3); % OK

  xjacm(1,1) = x*dxNloc';
  xjacm(1,2) = y*dxNloc';
  xjacm(2,1) = x*dyNloc';
  xjacm(2,2) = y*dyNloc';
  
  xjaci = inv(xjacm);

  area2 = abs(xjacm(1,1)*xjacm(2,2) - xjacm(2,1)*xjacm(1,2));
  area  = area2/2;
  
  dxN(1) = xjaci(1,1)*dxNloc(1)+xjaci(1,2)*dyNloc(1);
  dxN(2) = xjaci(1,1)*dxNloc(2)+xjaci(1,2)*dyNloc(2);
  dxN(3) = xjaci(1,1)*dxNloc(3)+xjaci(1,2)*dyNloc(3);
  dxN(4) = xjaci(1,1)*dxNloc(4)+xjaci(1,2)*dyNloc(4);
  dxN(5) = xjaci(1,1)*dxNloc(5)+xjaci(1,2)*dyNloc(5);
  dxN(6) = xjaci(1,1)*dxNloc(6)+xjaci(1,2)*dyNloc(6);
 
  dyN(1) = xjaci(2,1)*dxNloc(1)+xjaci(2,2)*dyNloc(1);
  dyN(2) = xjaci(2,1)*dxNloc(2)+xjaci(2,2)*dyNloc(2);
  dyN(3) = xjaci(2,1)*dxNloc(3)+xjaci(2,2)*dyNloc(3);
  dyN(4) = xjaci(2,1)*dxNloc(4)+xjaci(2,2)*dyNloc(4);
  dyN(5) = xjaci(2,1)*dxNloc(5)+xjaci(2,2)*dyNloc(5);
  dyN(6) = xjaci(2,1)*dxNloc(6)+xjaci(2,2)*dyNloc(6);

%==================

  bmat_b1  = [ 0,-dxN(1),     0  ;
               0,      0,-dyN(1) ;
               0,-dyN(1),-dxN(1)];
               
  bmat_b2  = [ 0,-dxN(2),     0  ;
               0,      0,-dyN(2) ;
               0,-dyN(2),-dxN(2)];
               
  bmat_b3  = [ 0,-dxN(3),     0  ;
               0,      0,-dyN(3) ;
               0,-dyN(3),-dxN(3)];
 
  bmat_b4  = [ 0,-dxN(4),     0  ;
               0,      0,-dyN(4) ;
               0,-dyN(4),-dxN(4)];

  bmat_b5  = [ 0,-dxN(5),     0  ;
               0,      0,-dyN(5) ;
               0,-dyN(5),-dxN(5)];

  bmat_b6  = [ 0,-dxN(6),     0  ;
               0,      0,-dyN(6) ;
               0,-dyN(6),-dxN(6)];

  bmat_b = [bmat_b1,bmat_b2,bmat_b3,bmat_b4,bmat_b5,bmat_b6];
  
%==================

 %== Puntos de colocacion:
 
   r3 = 1/sqrt(3);
   ap = 0.5 - r3;
   
   cx = [ 0.5-r3 , 0.5+r3,     ap,   1-ap ,     0 ,     0  ];
   cy = [      0 ,     0 , 1.0-ap,     ap ,  0.5+r3, 0.5-r3];
  
   c     = zeros(12,12);
   b_bar = [];
   for i = 1 : 6
      
     L1=1.0-cx(i)-cy(i);
     L2=cx(i);
     L3=cy(i);
 
     N(1) = (2*L1-1)*L1;    % OK
     N(2) = 4*L1*L2;        % OK
     N(3) = (2*L2-1)*L2;    % OK
     N(4) = 4*L2*L3;        % OK
     N(5) = (2*L3-1)*L3;    % OK
     N(6) = 4*L1*L3;        % OK
 
     dxNloc(1)=1.0-4.0*L1;  % OK
     dxNloc(2)=4.0*(L1-L2); % OK
     dxNloc(3)=4.0*L2-1.0;  % OK
     dxNloc(4)=4.0*L3;      % OK
     dxNloc(5)=0.0;         % OK
     dxNloc(6)=-4.0*L3;     % OK
 
     dyNloc(1)=1.0-4.0*L1;  % OK
     dyNloc(2)=-4.0*L2;     % OK
     dyNloc(3)=0.0;         % OK
     dyNloc(4)=4.0*L2;      % OK
     dyNloc(5)=4.0*L3-1.0;  % OK
     dyNloc(6)=4.0*(L1-L3); % OK
 
     xjacm(1,1) = x*dxNloc';
     xjacm(1,2) = y*dxNloc';
     xjacm(2,1) = x*dyNloc';
     xjacm(2,2) = y*dyNloc';
   
     xjacip = inv(xjacm);

     dxN(1) = xjacip(1,1)*dxNloc(1)+xjacip(1,2)*dyNloc(1);
     dxN(2) = xjacip(1,1)*dxNloc(2)+xjacip(1,2)*dyNloc(2);
     dxN(3) = xjacip(1,1)*dxNloc(3)+xjacip(1,2)*dyNloc(3);
     dxN(4) = xjacip(1,1)*dxNloc(4)+xjacip(1,2)*dyNloc(4);
     dxN(5) = xjacip(1,1)*dxNloc(5)+xjacip(1,2)*dyNloc(5);
     dxN(6) = xjacip(1,1)*dxNloc(6)+xjacip(1,2)*dyNloc(6);
 
     dyN(1) = xjacip(2,1)*dxNloc(1)+xjacip(2,2)*dyNloc(1);
     dyN(2) = xjacip(2,1)*dxNloc(2)+xjacip(2,2)*dyNloc(2);
     dyN(3) = xjacip(2,1)*dxNloc(3)+xjacip(2,2)*dyNloc(3);
     dyN(4) = xjacip(2,1)*dxNloc(4)+xjacip(2,2)*dyNloc(4);
     dyN(5) = xjacip(2,1)*dxNloc(5)+xjacip(2,2)*dyNloc(5);
     dyN(6) = xjacip(2,1)*dxNloc(6)+xjacip(2,2)*dyNloc(6);

     jpos = [ i*2-1 , i*2 ];
     c(jpos,jpos) = xjacm;
   
     bmat_s1  = [ dxN(1), -N(1),    0  ;
                  dyN(1),     0, -N(1)];
              
     bmat_s2  = [ dxN(2), -N(2),    0  ;
                  dyN(2),     0, -N(2)];
                
     bmat_s3  = [ dxN(3), -N(3),    0  ;
                  dyN(3),     0, -N(3)];
      
     bmat_s4  = [ dxN(4), -N(4),    0  ;
                  dyN(4),     0, -N(4)];
              
     bmat_s5  = [ dxN(5), -N(5),    0  ;
                  dyN(5),     0, -N(5)];
 
     bmat_s6  = [ dxN(6), -N(6),    0  ;
                  dyN(6),     0, -N(6)];
 
     bmat_s = [bmat_s1,bmat_s2,bmat_s3,bmat_s4,bmat_s5,bmat_s6];
   
     b_bar = [ b_bar ;
               bmat_s];
   end
  
  a = sqrt(2)/2;
  
  T_mat = [  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ;
             0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 ;
             0,  0,  0,  0, -a,  a,  0,  0,  0,  0,  0,  0 ;
             0,  0,  0,  0,  0,  0, -a,  a,  0,  0,  0,  0 ;
             0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0 ;
             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1 ];
 
  P_mat = [ 1,   cx(1),   cy(1),  0,       0,      0  ;
            1,   cx(2),   cy(2),  0,       0,      0  ;
           -a,-a*cx(3),-a*cy(3),  a, a*cx(3), a*cy(3) ;
           -a,-a*cx(4),-a*cy(4),  a, a*cx(4), a*cy(4) ;
            0,       0,       0,  1,   cx(5),   cy(5) ;
            0,       0,       0,  1,   cx(6),   cy(6) ];
           
       
 A_mat = [ 1 , xgs, ygs, 0 ,  0 ,   0 ;
           0 ,   0,   0, 1 , xgs, ygs ]; 
       
 bmat_s = xjaci * A_mat * inv(P_mat) * T_mat * c * b_bar;

