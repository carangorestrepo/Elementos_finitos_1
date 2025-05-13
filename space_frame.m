function [Pg,Pl,Hu]  =  space_frame
%             SPACE ELASTIC-FRAME SYSTEM ANALYSIS WITH LINER STATIC
%             METHOD.         
%                                      H.C.E.  Ali ?ZG?L (C) (R) 08-09-2007
% Description:
%
%             In this matlap application had descriptioned, space elastic
%             frame system of structural analysis as linear static methods.
%            
% Base-Algorithm:
% Example: [K] solution
%                         _______________________________
%                        | Local axis system per element | [K]_local
%                        |________stiffness matrix_______|
%                                        |
%                                        |                 [T]
%         ______________________________\|/____________________________
%        | Tranformation matrix for global to local axis system and    |
%        |____ Conversion matrix for local to global axis system and___|
%                                        | 
%              _________________________\|/_______________________
%              |   Global axis system per element stiffness matrix |
%              |___________________________________________________|
%                                        |                 [K]global
%              _________________________\|/_______________________
%              | System stiffness matrix with superposition       |
%              |__________________________________________________|
%                                        |                 [K]sis
%              _________________________\|/_______________________
%              | System loads vector  with superposition          |
%              |__________________________________________________|
%                                        |                [Q]sis & [P]sis
%              _________________________\|/_______________________
%              | System displacements [K]{D}+{f}={P}              |
%              |__________________________________________________|
%                                        |              [K]sis*{D}={P-Q}sis
%                                        |
%          _____________________________\|/____________________________
%          | Global axis system node reactions  {P}global              |
%          | Local axis system node reactions   {P}local =[T].{P}global|
%          |___________________________________________________________|
%
%
% Example   :     [Pg,Pl,Hu]  =  space_frame       <--|
% Pglobal =
%
%   [ 0.0000]     [ 0.0000 ]     [ 0.0000 ]     [ 0.0000 ]     fx
%   [ 6.7015]     [-6.7015 ]     [ 8.4030 ]     [-6.5970 ]     fy
%   [ 0.0000]     [-0.0000 ]     [ 0.0000 ]     [-0.0000 ]     fz
%   [16.7537]     [-16.7537]     [-0.0000 ]     [-0.0000 ]     Mxx
%   [-0.0000]     [-0.0000 ]     [-0.0000 ]     [ 0.0000 ]     Myy
%   [ 1.2092]1(i) [-1.2092 ]2(i) [ 2.4184 ]3(i) [-11.1934]4(i) Mzz
%
%   [-0.0000]     [-0.0000 ]     [-0.0000 ]     [-0.0000 ]   
%   [-6.7015]     [ 6.7015 ]     [1.5970  ]     [ 6.5970 ]
%   [-0.0000]     [ 0.0000 ]     [-0.0000 ]     [ 0.0000 ]
%   [16.7537]     [-16.7537]     [ 0.0000 ]     [-0.0000 ]
%   [-0.0000]     [-0.0000 ]     [-0.0000 ]     [ 0.0000 ]
%   [-1.2092]1(j) [ 1.2092 ]2(j) [ 11.1934]4(j) [-15.1947]4(j)
%   
% References:
%{
[1]  Robert D. Cook, David S. Malkus, Michael E. Plesha, Concepts and 
     Applications of Finite Element Analysis, John-Wiley Sons,1989,
     no:88-27929.
 
[2]  T.Y. Yang, Finite Element Structural Analysis, Prentice-Hall 
     International Series in Civil Engineering and Engineering Mechanics,
     1986, ISBN:0-13-317116-7.
 
[3]  Kasimzade A.A., Finite Element Method : Foundation and Application to
     Earthquake Engineering (is included education and finite element 
     analysis programs CD), Istanbul, Beta Publication,(First edition 1997)
     Second edition, 2004, p.827.(ISBN 975-511-379-7).
 
[4]  Kasimzade A.A., Structural Dynamics : Theory and Application to 
     Earthquake Engineering (is included education and dynamic analysis 
     programs CD) , Istanbul, Beta Publication, (First edition 1998) Second
     edition, 2004, p.527. (ISBN 975-511-381-9).
 
[5]  J.N.Reddy, D.K. Gartling, The Finite Element Method in Heat Transfer
     and Fluid Dynamics, Second edition,CRC Press,no:DE-AC04-76DP00789.
 
[6]  Chuen-Yuan Chia, Nonlinear Analysis of Plates, McGraww-Hill Press,1980,
     ISBN: 0-07-010746-7.
 
[7]  Harry Kraus, Thin Elastic Shells, John Wiley &Sons inc.,1967,Library
     of Congress Catalog Card Number: 67-23328.
 
[8]  William Weaver, Paul R. Johnston, Finite Elementsd for Structural 
     Analysis Prentice-Hall Inc.,1984, ISBN: 0-13-317099-3.
 
 
[9]  A. Laulusa,O.A. Bauchau ,J-Y. Choi,V.B.C. Tan c, L. Li, Evaluation of 
     some shear deformable shell elements, Elsevier  Science, 18 October 
     2005,pp:43-(2006)-5033-5054
 
[10] Kasimzade A.A., Theory of Elasticity and Structural Analysis of space 
     frame, membrane, plate, shell, asolid, solid systems (is included 
     education and finite element analysis programs-2 diskettes ) , 
     Istanbul ,"Beta" Publication, 2000, p.401 . (ISBN 975-486-866-2).
%}
global node_coordinates
global node_position
global freedom
clc
%____________________________
% INPUT-VERIABLES [1,2,3,4,5]
%_________________|
%
%[1]____________[ GLOBAL COORDINATES ]
%                    {x_g}    {y_g}  {z_g}
node_coordinates = [ 0.00    0.00   5.00
                     0.00    0.00   0.00
                     0.00    0.00  -5.00
                     4.00    0.00   0.00
                     8.00    3.00   0.00];
%___________________________________|
Node=size(node_coordinates,1);
%[2]______________[ POSITION MATRIX ]
%          Node(i)--->--Node(j)
node_position = [1    2 
                 2    3
                 2    4
                 4    5];
                 
%___________________________________|
No=size(node_position,1);
%[3]_____________[ MATERIAL PROPERTIES ]
% E_x ,E_z, G        (KN/m**2)                (SI)
% I_x ,Iz ,I_xz ,Jp    (m**4)
% A  = cross sectional area   (m**2)
%      E_x    E_z    G    I_x   I_z     I_xz    A            Jp
m_p = [ 2E+6  2E+6 1E+5 2.6E-3 0.651E-3 0.00  0.125   ];   
%  frame_beta (B1) load_beta (B2)
p_p =  [  0.00      0.00 ];
m_p = m_p';
p_p = p_p';
 
%if all element properties is equal than. 
for i=1:No
 p_p(:,i)=p_p(:,1);
 m_p(:,i)=m_p(:,1);
 
 %Polar inertia moment
% J_p = I_x + I_z - I_xz
  J_p = m_p(4,i) + m_p(5,i) - m_p(6,i);
  m_p(8,i)=J_p;   
end
%__________________________________|
for i=1:Node;    Re(i,:)=[1 1 1 1 1 1]; end
Nom=size(Re,2); %One node's total degree of freedom.
%[4]______________[ SYSTEM SUPPORT ]
%Re(Node number,:)=[ u(x) v(y) w(z) Qxx Qyy Qzz]
 Re(1,:)=[1 0 0 0 0 0];
 Re(3,:)=[1 0 0 0 0 0];
 Re(5,:)=[0 0 0 0 0 0];
%_________________________________|
% if not is empty P and load_matrix than hide this line else
% active this line for ;
[freedom ,R,Re]  =  frame_element_topology (No,Node,Nom,Re);
P(freedom) = 0;
%load_matrix = [];
%[5]_______________[ SYSTEM LOAD ]
%Re(Node number, freedoom number)[fx fy fz Mxx Myy Mzz]
 P(Re(2,2)) =  -5.00;
 P(Re(4,2)) =  -5.00;
%                                                                     |---|
%             Element_no  qo       P  a   b   L  logical , load_type ,load_case )
 load_matrix=[   3      , 2.00 , 0.00 ,0 , 0 , 5 , true   , 1    ,  2   ];
%                2      , 2.00 , 0.00 ,0 , 0 , 5 , true   , 1    ,  1   ];
%                4      , 4.00 , 0.00 ,0 , 0 , 5 , true   , 12   ,  1   ];
%________________________________|
%Analysis Space-Frame system
%[freedom ,R]            =  frame_element_topology (No,Node,Nom,Re);
[global_loads,  Qsis ]  =  frame_system_loads (R,load_matrix,p_p ,freedom ,No);
[T,K_g,Hu]              =  frame_system_stiffness (m_p,p_p,No,P,R,Qsis) ;
[Pg,Pl]                 =  frame_node_reactions(No,K_g,T,Hu,global_loads);
%==========================
% PROGRAM SUB-FUNCTIONS
%==========================
%____________________________________________________________|||||||||||||
function [freedom,R,Re] = frame_element_topology (No,Node,Nom,Re)
% 
% Description:
%               In this sub-function calculated per element reology matrix
%
% Syntax:       
%               Node = System total node value
%                Nom = one node total freedom value.
    
    
global node_position
%Accumulate method
freedom=0;
 for i=1:Node;
     for j=1:Nom;
         if Re(i,j)==1 ;
             freedom = freedom +1;
             Re(i,j) = freedom;
         end
     end
 end
  
%Topology method.
for i=1:No
    R(i,:)=[Re(node_position(i,1),:) Re(node_position(i,2),:) ];
end
%________________________________|
%________________________________________________________________||||||||||
function [global_loads, Qsis ] = ...
          frame_system_loads (R,load_matrix,p_p ,freedom ,No)
%
% Description:
%
%             In sub-function produced system internal-load vector with
%             global and local node reactions.
%
% Syntax :
%             load_matrix = Uniform loads description matrix
%             p_p = beta angle matrix
%         freedom = total system freedom
%              No = total element no
%                                                                [ERROR ]
if isempty(load_matrix) == true
global_loads= 0;
Qsis = 0;
return
end
    
    
%Open dimension for cumulative variables.
Qsis(freedom)      = 0 ;
global_loads(No,:) = zeros(1,12);
local_loads(No,:)  = zeros(1,12);
for i=1:size(load_matrix,1)
    if load_matrix(i,9)==1     %look system_load_1.pdf
    load_vector = uniform_load_1( load_matrix(i,2), ...    
                                  load_matrix(i,3), ...    
                                  load_matrix(i,4), ...    
                                  load_matrix(i,5), ...    
                                  load_matrix(i,6), ...    
                                  load_matrix(i,7), ...    
                                  load_matrix(i,8));
                              
    else                        %look system_load_2.pdf
    load_vector = uniform_load_2( load_matrix(i,2), ...    
                                  load_matrix(i,3), ...    
                                  load_matrix(i,4), ...    
                                  load_matrix(i,5), ...    
                                  load_matrix(i,6), ...    
                                  load_matrix(i,7), ...    
                                  load_matrix(i,8));    
    end
%degree    
betan = p_p(2,i);
s     = load_matrix(i,1);
%[loads_transformation_matrix] = load_transformation(s,betan);
T = load_transformation(s,betan);
%global load_vector                        
 global_load_vector = T'*load_vector;
 
%This matrix need global and local system node reactions
global_loads(load_matrix(i,1),:) = global_load_vector;
 local_loads(load_matrix(i,1),:) = load_vector;
 
%Global axis per element uniform out load superposition
%Qsis(freedom)=0;
      for n=1:No;
           for sat=1:12;
               if (R(n,sat)~=0)
                   Qsis(R(n,sat))=Qsis(R(n,sat)) + global_loads(n,sat);
              end
          end
      end
end % if-than
%____________________________________|
%________________________________________________________________||||||||||
function [T,K_g,Hu] = frame_system_stiffness (m_p,p_p,No,P,R,Qsis)
%
% Description:
%               In this sub-function calculated frame-element local and
%               global system stiffness matrix, system stiffness matrix,
%               global to local axis system transformation matrix and
%               global nodes model displacement matrix [Hu].
%
% Syntax:
%             m_p = material properties matrix.
%             p_p = axial rotation matrix angle.
%           betan = axial rotation angle
%               s = elemet no.
%
%            K_l = Local axis system stiffness matrix
%              T = Tranformation matrix
%            K_g = Global axis system per element global stiffness matrix
%           Ksis = System  stiffness matrix
%             Hu = Per nodes global system modal displacement matrix
global freedom
for s=1:No
%Material base properties    
    E_x  = m_p(1,s) ;  E_z   = m_p(2,s)  ; G   = m_p(3,s)    ; 
    I_x  = m_p(4,s) ;  I_z   = m_p(5,s)  ; J_p = m_p(8,s)    ;
    A    = m_p(7,s) ;  betan = p_p(1,s)  ;
%Global--->Local transformation matrix [T]
    [T(:,:,s) L]  = frame_transformation( betan , s) ;
%Local axis system stiffness matrix
    K_l(:,:,s) = local_stiffness(A,L,E_x,E_z,G,I_x,I_z,J_p);
%Global axis system stiffness matrix
    K_g(:,:,s)=T(:,:,s)'*K_l(:,:,s)*T(:,:,s);
end %--s
%System stiffness matrix superposing with topology method.
Ksis(freedom,freedom)=0;
     for n=1:No;
         for sat=1:12;
             for sut=1:12;
               if (R(n,sat)~=0)
                  if (R(n,sut)~=0);
          Ksis(R(n,sut),R(n,sat))=Ksis(R(n,sut),R(n,sat)) + K_g(sat,sut,n);
                  end    
               end
            end
         end
     end
%System stiffness matrix singularity control                    [ERROR - 2]
if size(Ksis,2) == rank(Ksis)
D = Ksis^-1*(P-Qsis)';
else
disp('Warning: System stiffness matrix solution error')
disp('Control the system boundary conditions')
disp(R);
return
end
%Moving per node freedom-displacement global system
for  v = 1 : No;
      for m = 1 :12;
          u = R(v, m);
          if u ~=0 
             Hu(v, m,:) = D(u) ;
          else    
             Hu(v,m)=0;
          end
      end
end
%_____________________________|
%_______________________________________________________________|||||||||||
function [Pg,Pl] = frame_node_reactions(No,K_g,T,Hu,global_loads)
%
% Description:
%               In this sub-function calculated frame element has
%               global and local node reactions
% Syntax:
%             No = Total frame-element no
%             Kg = Global axis system stiffness matrix
%              T = Transformation matrix global--->local
%             Hu = Global system modal displacement matrix 
%   global_loads = Uniform load edge node reactions on frame element
for s=1 :No;
% [K]{D} + {f} = {P}    
%Global reactions
if global_loads==0;
     Pg(:,s) = K_g(:,:,s)*Hu(s,:)';
else
     Pg(:,s) = K_g(:,:,s)*Hu(s,:)'+global_loads(s,:)';
end
%Local axis reactions.
%local_loads(:,v)=T(:,:,v)*global_loads(v,:)
     Pl(:,s) = T(:,:,s)*Pg(:,s);
end
%__________________________________|
    
%________________________________________________________________||||||||||
function [node_coor_x , node_coor_y ,node_coor_z]=frame_node_move(s)
% Description:
%            In this matlab sub-function moved global system nodes to
%            local frame-element (i) and (j) nodes.
%
% Syntax:
%            s = Frame element no:
%            Here other variables are global 
global node_coordinates
global node_position
      
%Moving global-axis-system coordinates frame local axis
node_coor_x = [node_coordinates(node_position(s,1),1) ...
               node_coordinates(node_position(s,2),1)];
        
node_coor_y = [node_coordinates(node_position(s,1),2) ...
               node_coordinates(node_position(s,2),2)];
         
node_coor_z = [node_coordinates(node_position(s,1),3) ...
               node_coordinates(node_position(s,2),3)];
%_______________________________|                 
%________________________________________________________________||||||||||
function [loads_transformation_matrix] = load_transformation(s,betan)
% Description:
%               In this sub-function is produce tranformation matrix for
%               local axis system of uniform load reactions.
%
% Syntax:   
%             s = Uniform load on frame element no
%         betan = Rotation matrix for uniform loads x-local axis 
%global node_coordinates
%global node_position
[load_coor_x  load_coor_y load_coor_z]= frame_node_move(s);
%Cartesian coordinates
x1 = load_coor_x(1);   x2 = load_coor_x(2);
y1 = load_coor_y(1);   y2 = load_coor_y(2);
z1 = load_coor_z(1);   z2 = load_coor_z(2);
%Direction cosinesses
L = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2) ;
ly     = (x2-x1) / L ;
my     = (y2-y1) / L ;
ny     = (z2-z1) / L ;
Q      = sqrt(1 - ny ^ 2)  ;
%radii<-----degree
betan=betan*pi/180;
%"x" axis direction direction cosinesses for load position
%base_transform = [ly  my   ny
%                 -my  ly   0
%                 -ny  0   ly ];
switch abs(ny)
    case{+1}
base_transform = [ 0.00  0.00  1.00      
                   0.00  1.00  0.00
                   1.00  0.00  0.00 ];
    otherwise
        
base_transform = [ ly          my           ny
                   my / Q      ly / Q      0.00 
                  -ly*ny / Q   -my*ny / Q   Q    ]';            
end
%x local axis rotation matrix
load_rotation = [ 1.00   0.00  0.00
                  0.00  cos(betan) -sin(betan)
                  0.00  sin(betan) cos(betan)];
                  
%Multi transformation
Tn=load_rotation*base_transform;              
%Frame element [T]
loads_transformation_matrix= [  Tn        zeros(3,3)  zeros(3,3) zeros(3,3)
                              zeros(3,3)   Tn         zeros(3,3) zeros(3,3)
                              zeros(3,3)  zeros(3,3)  Tn         zeros(3,3)
                              zeros(3,3)  zeros(3,3)  zeros(3,3)   Tn  ];
                          
%____________________________________|                          
%_______________________________________________________________|||||||||||
%
% Description:
%         In this sub-function had produced of elastic-space-frame element
%         conversion matrix local to global axis system ( [C] matrix )
%         transform matrix  global to local axis system ( [T] matrix )
%         for local axis "y" direction ny==-+1 of special station and
%         ny<>-+1 general frame of space-position.
% 
% Syntax: node_coordinates: System all node's global cartesian coordinates
%         node_position   : System all frame-element (i) and (j) node's
%                           position matrix.
%
%         betan            : Frame-element has surface rotation angle as
%                           right-hand-rule.
%
function [frame_transformation_matrix  L] = frame_transformation(betan,s)
%global node_coordinates
%global node_position
[node_coordinates_x ,...
 node_coordinates_y ,...
 node_coordinates_z]= frame_node_move(s);
%Cartesian coordinates
x1 = node_coordinates_x(1);   x2 = node_coordinates_x(2);
y1 = node_coordinates_y(1);   y2 = node_coordinates_y(2);
z1 = node_coordinates_z(1);   z2 = node_coordinates_z(2);
%Direction cosinesses
L = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2) ;
ly     = (x2-x1) / L ;
my     = (y2-y1) / L ;
ny     = (z2-z1) / L ;
Q      = sqrt(1 - ny ^ 2)  ;
%radii<-----degree
betan=betan*pi/180;
                            
                            
                                                     % |||||||
                                                        % |_______ny==+1
if ny==1 ; %In special conversion for ny=+1 and Q=0 than
% Local axis system->>[C]-->> Special axis system             
% x_s   =   x_o*cos(betan) - z_o*sin(betan)
% y_s   =   y_o
% z_s   =   x_o*sin(betan) + z_o*cos(betan)
%Special station of frame_transformation matrix 
% {X}         [0.00     0.00    1.00] {x_s}
% { }         [                     ] {   } 
% {Y}      =  [1.00     0.00    0.00]*{y_s}   
% { }         [                     ] {  }
% {Z}global   [0.00     1.00    0.00] {z_s}local
% {f}global = [C].{f}local
%(Conversion matrix local------>global)
   frame_rotation_matrix = [cos(betan)   0.00    -sin(betan)
                              0.00       1.00     0.00
                            sin(betan)   0.00     cos(betan)];            
                                    
                                     
         frame_transform = [  1.00     0.00    0.00
                              0.00     0.00   -1.00
                             0.00     1.00    0.00 ];
%Conversion matrix local to global axis system
frame_conversion_matrix = frame_transform*frame_rotation_matrix;
%Transformation matrix global to local axis system
frame_trans = frame_conversion_matrix^-1;
end
                                                      % |||||||
                                                         % |_______ny==-1
if ny==-1 %In special conversion for ny=-1 and Q=0 than
             
%Conversion matrix for axial (y-y) rotation as betan (right-hand-rule)
        frame_rotation_matrix =[ cos(betan)  0.00   sin(betan)
                                    0.00     1.00    0.00 
                                -sin(betan)  0.00   cos(betan) ];
                                       
            frame_transform  = [ 1.00   0.00   0.00
                                 0.00   0.00   1.00
                                 0.00  -1.00   0.00];
%{x'}               {xo}
%{y'} = [C]Rotation*{yo}  =>  [XYZ]'=[C]r*[XYZ]o
%{z'}special        {zo}local
%
%{xo}                {X}
%{yo} = [C]conversion{Y}  =>  [XYZ]o=[C]c*[XYZ]g =>
%{zo}special         {Z}global
%
%[XYZ]special = [C]r*[C]c*[XYZ]global  (Conversion method)
                                      
                                      
%Conversion matrix local to global axis system
frame_conversion_matrix = frame_transform*frame_rotation_matrix;
%Transformation matrix global to local axis system
frame_trans = frame_conversion_matrix^-1;
end
                                             % |||||||
                                                % |_______ny<>-1 and ny<>+1
if abs(ny)~= +1  %General conversion for ny<>-1 and ny<>+1,Q<>0 than
%Conversion matrix for axial (y-y) rotation as betan (right-hand-rule)
% my:"y" axis direction direction cosines
            frame_rotation_matrix =[ cos(betan)    0.00    -my*sin(betan)
                                      0.00       1.00        0.00 
                                     my*sin(betan)  0.00      cos(betan) ];
%Conversion matrix [C]=[T]'=[T]^-1
frame_transform  = [  my / Q      -ly / Q       0.00 
                       ly           my           ny 
                     -ly*ny / Q   -my*ny / Q    Q    ]';
%Conversion matrix local-->global      [C]
frame_conversion_matrix = frame_transform*frame_rotation_matrix;
%Transformation matrix global--->local [T]
frame_trans = frame_conversion_matrix^-1;
end
zero = zeros(3,3);   %          f[x,y,z]i  M[x,y,z]i  f[x,y,z]j  M[x,y,z]j  
frame_transformation_matrix = [ frame_trans  zero         zero        zero 
                                 zero        frame_trans  zero        zero  
                                 zero        zero         frame_trans zero
                                 zero        zero         zero        frame_trans];
%________________________________|
%_______________________________________________________________||||||||||
function frame_local_stiffness = local_stiffness(A,L,E_x,E_z,G,I_x,I_z,J_p)
%
% Description:  
%               In this sub-function is produce elastic space-frame
%               element has local axis stiffness matrix.
%
% Syntax:      
%              A = Cross section area  (m**2)        (SI).
%            E_x = "x" axis elasticity module        (KN/m**2).
%            E_z = "z" axis elasticity module        (KN/m**2).
%             G  = Shear elasticity module           (KN/m**2).
%            I_x = "x" axis inertia moment           (m**4).
%            I_z = "y" axis inertia moment           (m**4).
%            J_p = I_x+I_z-I_xz Polar inertia moment (m**4).
%             L  = Element's length.                 (m).
%         u_x(i)               v_y(i)        w_z(i)           Qx_x(i)         Qy_y(i)      Qz_z(i)      
K_i_i = [ 12*E_z*I_z/(L^3)     0.00          0.00              0.00             0.00    -6*E_z*I_z/(L^2)   %u_x (i)
          0.00                A*E_x/L        0.00              0.00             0.00       0.00            %v_y (i) 
          0.00                 0.00        12*E_x*I_x/(L^3)  6*E_x*I_x/(L^2)    0.00       0.00            %w_z (i)
          0.00                 0.00        6*E_x*I_x/(L^2)   4*E_x*I_x/L        0.00       0.00            %Qx_x(i) 
          0.00                 0.00          0.00              0.00            G*J_p/L     0.00            %Qy_y(i)
         -6*E_z*I_z/(L^2)      0.00          0.00              0.00             0.00      4*E_z*I_z/L   ]; %Qz_z(i)
%         u_x(j)               v_y(j)        w_z(j)           Qx_x(j)         Qy_y(j)      Qz_z(j)           
K_i_j = [ -12*E_z*I_z/(L^3)     0.00       0.00              0.00            0.00      -6*E_z*I_z/(L^2)   %u_x (i)
           0.00               -A*E_x/L     0.00              0.00            0.00         0.00            %v_y (i)
           0.00                 0.00    -12*E_x*I_x/(L^3)  6*E_x*I_x/(L^2)   0.00         0.00            %w_z (i)
           0.00                 0.00    -6*E_x*I_x/(L^2)   2*E_x*I_x/L       0.00         0.00            %Qx-x(i)
           0.00                 0.00       0.00              0.00           -G*J_p/L      0.00            %Qy_y(i)
          6*E_z*I_z/(L^2)       0.00       0.00              0.00            0.00       2*E_z*I_z/L  ];   %Qz_z(i)
 
 
%         u_x(i)               v_y(i)        w_z(i)           Qx_x(i)         Qy_y(i)      Qz_z(i)                 
K_j_i = [ -12*E_z*I_z/(L^3)     0.00         0.00              0.00            0.00    6*E_z*I_z/(L^2)    %u_x (j)
           0.00              -A*E_x/L        0.00              0.00            0.00      0.00             %v_y (j)
           0.00                 0.00    -12*E_x*I_x/(L^3)  -6*E_x*I_x/(L^2)    0.00      0.00             %w_z (j)
           0.00                 0.00      6*E_x*I_x/(L^2)   2*E_x*I_x/L        0.00      0.00             %Qx-x(j)
           0.00                 0.00         0.00              0.00          -G*J_p/L    0.00             %Qy-y(j)
          -6*E_z*I_z/(L^2)      0.00         0.00              0.00            0.00      2*E_z*I_z/L ];   %Qz-z(j)
%         u_x(j)               v_y(j)        w_z(j)           Qx_x(j)         Qy_y(j)      Qz_z(j)                 
K_j_j = [  12*E_z*I_z/(L^3)    0.00        0.00               0.00            0.00    6*E_z*I_z/(L^2)    %u_x (j)
           0.00                A*E_x/L     0.00               0.00            0.00      0.00             %v_y (j)
           0.00                0.00      12*E_x*I_x/(L^3)  -6*E_x*I_x/(L^2)   0.00      0.00             %w_z (j)
           0.00                0.00      -6*E_x*I_x/(L^2)   4*E_x*I_x/L       0.00      0.00             %Qx_x(j)
           0.00                0.00        0.00               0.00           G*J_p/L    0.00             %Qy_y(j)
         6*E_z*I_z/(L^2)       0.00        0.00               0.00            0.00      4*E_z*I_z/L  ];  %Qz_z(j)
     
     
%Space frame-element local axis system stiffness matrix     
frame_local_stiffness = [ K_i_i   K_i_j
                          K_j_i   K_j_j ];
                      
%__________________________________|
%_________________________________________________________________|||||||||
function uniform_load_vector = uniform_load_1(qo,P,a,b,L,logical,load_type)
%
% Description:
%
%              In this function had calculated internal force effects under
%              the uniform or single load on elastic frame system. 
%
% Syntax:  
%              qo = Uniform or not-uniform distribute load (KN/m)
%              P  = Single loads
%          a,b,L  = Frame element lenght components
%        logical  = Frame system load position
%                   if you see (i)-->---(j) logical = true  or =1
%                   if you see (i)--<---(j) logical = false or =0
%
%    load_type = Please look loads combinations in "load_types1.pdf"
%
switch load_type
    case {1}
            fx_i = 0.00      ;
            fy_i = 1/2*qo*L  ;
            Mz_i = 1/8*qo*L^2;
            fx_j = 0.00      ;
            fy_j = 1/2*qo*L  ;
            Mz_j = 0.00      ;
    case{2}
            fx_i = 0.00       ;
            fy_i = qo*a-(qo*a^2)/(2*L)  ;
            Mz_i = 1/8*qo*a^2*(2-a/L)^2 ;
            fx_j = 0.00        ;
            fy_j = 1/2*qo*a^2/L;
            Mz_j = 0.00        ;
    case{3}
      c=L;L=a+b;
            fx_i = 0.00     ; 
            fy_i = 1/2*qo*c ;
            Mz_i = 1/(16*L)*qo*c*(3*L^2-c^2) ;
            fx_j = 0.00     ;
            fy_j = 1/2*qo*c ;
            Mz_j = 0.00     ;
    case{4}
      c=L;L=a+b;
            fx_i = 0.00;
            fy_i = qo*b*c/L;
            Mz_i = 1/8*qo*b*c/L^2*(4*(L^2-b^2)-c^2);
            fx_j = 0.00;
            fy_j = qo*b*c/L;
            Mz_j = 0.00;
    case{5}
            fx_i = 0.00 ;
            fy_i = qo*a ;
            Mz_i = 1/(4*L)*qo*a^2*(3*L-2*a);
            fx_j = 0.00 ;
            fy_j = qo*a ;
            Mz_j = 0.00 ;
    case{6}
            fx_i = 0.00;
            fy_i = 1/4*qo*L;
            Mz_i = 5/64*qo*L^2;
            fx_j = 0.00;
            fy_j = 1/4*qo*L;
            Mz_j = 0.00;
    case{7}
            fx_i = 0.00;
            fy_i = qo*a/2-qo*a^2/(3*L)+qo*b^2/(3*L);
            Mz_i = qo*L/120*(L+b)*(7-3*b^2/L^2);
            fx_j = 0.00;
            fy_j = qo*b/2 + qo*a^2/(3*L)-qo*b^2/(3*L);
            Mz_j = 0.00;
    case{8}
            fx_i = 0.00;
            fy_i = 1/3*qo*L;
            Mz_i = 1/15*qo*L^2;
            fx_j = 0.00;
            fy_j = 1/6*qo*L;
            Mz_j = 0.00;
    case{9}
            fx_i = 0.00;
            fy_i = 1/2*qo*a;
            Mz_i = 1/(8*L)*qo*a^2*(2*L-a);
            fx_j = 0.00;
            fy_j = 1/2*qo*a;
            Mz_j = 0.00;
    case{10}
            fx_i = 0.00;
            fy_i = qo*L;
            Mz_i = 1/8*qo*(L^2-a^2*(2-a/L));
            fx_j = 0.00;
            fy_j = qo*L;
            Mz_j = 0.00;
    case{11}
            fx_i = 0.00;
            fy_i = qo(1)*L/3 + qo(2)*L/6;
            Mz_i = L^2/120*(8*qo(1)+7*qo(2));
            fx_j = 0.00;
            fy_j = qo(1)*L/6 + qo(2)*L/3;
            Mz_j = 0.00;
    case{12}
            fx_i = 0.00;
            fy_i = qo*L/pi;
            Mz_i = 1/10*qo*L^2;
            fx_j = 0.00;
            fy_j = qo*L/pi;
            Mz_j = 0.00;
    case{13}
            fx_i = 0.00;
            fy_i = 1/(L^3)*(P*b+1/2*P*a*b*(b+L));
            Mz_i = (1/2*P*a*b*(b+L))/(L^2);
            fx_j = 0.00;
            fy_j = 1/(L^3)*(P*a-1/2*P*a*b*(b+L));
            Mz_j = 0.00;
    case{14}
            %qo=n
            %
            fx_i = 0.00;
            fy_i = 1/qo*sum(P*(1:qo-1));
            Mz_i = 1/8*P*L*(1-1/qo);
            fx_j = 0.00;
            fy_j = 1/qo*sum(P*(1:qo-1));
            Mz_j = 0.00;
end
        
uniform_load_vector = load_logical(fx_i,fy_i,Mz_i,fx_j,fy_j,Mz_j,logical);
%______________________________|
%________________________________________________________________|||||||||
function uniform_load_vector = uniform_load_2(qo,P,a,b,L,logical,load_type)
%
% Description:
%
%              In this function is calculate internal force effects under
%              the uniform or single load on elastic frame system. 
%
% Syntax:  
%              qo = Uniform or not-uniform distribute load (KN/m)
%              P  = Single loads
%          a,b,L  = Frame element lenght components
%        logical  = Frame system load position
%                   if you see (i)-->---(j) logical = true  or =1
%                   if you see (i)--<---(j) logical = false or =0
%
%    load_type = Please look loads combinations in "load_types2.pdf"
%
switch load_type
    case {1}
            fx_i = 0.00;
            fy_i = 1/2*qo*L;
            Mz_i = 1/12*qo*L^2;
            fx_j = 0.00;
            fy_j = 1/2*qo*L;
            Mz_j =-1/12*qo*L^2;
    case{2}
            fx_i = 0.00 ;
            fy_i = qo*a-(qo*a^2)/(2*L);
            Mz_i = qo*a^2/4*(2-a*L*(8/3-a/L));
            fx_j = 0.00;
            fy_j = qo*a/(2*L);
            Mz_j =-qo/a^2/(12*L^2)*(4*L-3*a);
    case{3}
        c=L/qo;
            fx_i = 0.00;
            fy_i = 1/2*qo*c;
            Mz_i = 1/(24*L)*qo*c*(3*L^2-c^2);
            fx_j = 0.00;
            fy_j = 1/2*qo*c;
            Mz_j =-1/(24*L)*qo*c*(3*L^2-c^2);
    case{4}
        c=L;L=a+b;
            fx_i = 0.00;
            fy_i = qo*c/L*(b-c/2);
            Mz_i = qo*c/(12*L^2)*((4*L^2-c^2)*(2*b-a)-4*(2*b^3-a^3));
            fx_j = 0.00;
            fy_j = qo*c/L(a-c*2);
            Mz_j =-qo*c/(12*L^2)*((4*L^2-c^2)*(2*b-a)-4*(2*a^3-b^3));
    case{5}
            fx_i = 0.00;
            fy_i = qo*a;
            Mz_i = 1/(6*L)*qo*a^2*(3*L-2*a);
            fx_j = 0.00;
            fy_j = qo*a;
            Mz_j =-1/(6*L)*qo*a^2*(3*L-2*a);
    case{6}
            fx_i = 0.00;
            fy_i = 1/4*qo*L;
            Mz_i = 5/64*qo*L;
            fx_j = 0.00;
            fy_j = 1/4*qo*L;
            Mz_j =-5/64*qo*L^2;
    case{7}
            fx_i = 0.00;
            fy_i = 1/3*qo*(b^2-a^2)+qo*a/2;
            Mz_i = qo/(180*L)*(7*L^3-7*L^2*(2*b-a)+3*L*(2*a^2-b^2)-3*(2*a^3-b^3));
            fx_j = 0.00;
            fy_j = 1/3*qo*(b^2-a^2)+qo*a/2;
            Mz_j =-qo/(180*L)*(7*L^3+7*L^2*(2*b-a)-3*L*(2*a^2-b^2)-3*(2*a^3-b^3));
    case{8}
            fx_i = 0.00;
            fy_i = 2/3*qo*L;
            Mz_i = 1/20*qo*L^2;
            fx_j = 0.00;
            fy_j = 1/3*qo*L;
            Mz_j =-1/30*qo*L^2;
    case{9}
            fx_i = 0.00;
            fy_i = 1/2*qo*a;
            Mz_i = 1/(12*L)*qo*a^2*(2*L-a);
            fx_j = 0.00;
            fy_j = 1/2*qo*a;
            Mz_j =-1/(12*L)*qo*a^2*(2*L-a);
    case{10}
            fx_i = 0.00;
            fy_i = qo*a/2-qo*L/2;
            Mz_i = 1/12*qo*(L^2-a^2*(2-a/L));
            fx_j = 0.00;
            fy_j = qo*L;
            Mz_j =-1/(12)*qo*(L^2-a^2*(2-a/L));
    case{11}
            fx_i = 0.00;
            fy_i = qo(1)*L/6+qo(2)*L/3;
            Mz_i = L^2/60*(3*qo(1)+qo(2));
            fx_j = 0.00;
            fy_j = qo(1)*L/3 + qo(2)*L/6;
            Mz_j =-L^2/60*(2*qo(1)+3*qo(2));
    case{12}
            fx_i = 0.00;
            fy_i = qo*L/pi;
            Mz_i = 1/15*qo*L^2;
            fx_j = 0.00;
            fy_j = qo*L/pi;
            Mz_j =-1/15*qo*L^2;
    case{13}
            fx_i = 0.00;
            fy_i = P*b/L;
            Mz_i = P*a*b^2/L^2;
            fx_j = 0.00;
            fy_j = P*a/L;
            Mz_j =-P*b*a^2/L^2;
    case{14}
            fx_i = 0.00;
            fy_i = 1/qo*sum(P*(1:qo-1));
            Mz_i = 1/12*P*L*(qo+1/(2*qo));
            fx_j = 0.00;
            fy_j = 1/qo*sum(P*(1:qo-1));
            Mz_j =-1/12*P*L*(qo+1/(2*qo));
end
uniform_load_vector = load_logical(fx_i,fy_i,Mz_i,fx_j,fy_j,Mz_j,logical);
%________________________________|
%_________________________________________________________________|||||||||
function uniform_load_vector_n = load_logical(fx_i , fy_i , Mz_i ,...
                                              fx_j , fy_j , Mz_j ,...
                                              logical) 
                                          
% Desctiption:
%            In this sub-function is selecting frame-element's uniform load
%            station. Here is possible two different station for loads.
%            logical = true 
%                          (i)-----<---||||support-(j) node edges
%            logical = false
%                          (i)||||support----->--(j) node edges
%
% Syntax:  
%           fx,fy...Mz = Local axis system uniform load reactions
%           logical    = Local axis system uniform load station.
                                          
switch logical
    case {false}
        % (i)-----<---||||support-(j) node edges
             uniform_load_vector_n = [fx_j fy_j 0.00 0.00 0.00 -Mz_j ...
                                      fx_i fy_i 0.00 0.00 0.00 -Mz_i ]';
    case {true}             
        % (i)||||support----->--(j) node edges
             uniform_load_vector_n = [fx_i fy_i 0.00 0.00 0.00 Mz_i ...
                                      fx_j fy_j 0.00 0.00 0.00 Mz_j]';
             
end 
%_________________________________|