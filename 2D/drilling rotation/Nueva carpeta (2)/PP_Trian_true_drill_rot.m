% =========================================================================
%
%   Main FEM Program for 2D Elasticity : T3_true_DR  
%
% =========================================================================
clear all;
close all;
clc;

global E nu th nne ddln
format long

% read Problem DATA
[coor,connex,climit,charg,E,nu,th]= cantbeam32x8;

%[coor,connex,climit,charg,E,nu,th]=pure_bending32x8
%[coor,connex,climit,charg,E,nu,th]=pure_bending25x4

% Problem Initiation .....
NTN= size(coor,1) ; % nombre de noeuds
NTE = size(connex,1) ; % nombre d'éléments
ddln = 3 ; % nombre de degres de liberte par noeud
nne = 3 ; % nombre de noeuds par element
NTDDL = NTN*ddln ; % nombre total de degres de liberte du systeme
FG = zeros(NTDDL,1) ;  % vecteur forces nodales
KG = zeros(NTDDL,NTDDL) ;   % matrice de rigidité golobale
DG = zeros(NTDDL,1);   % vecteur déplacement nodaux

% Construction de la matrice de raideur globale
for iel = 1:NTE
% matrice de rotation
[rog,xe] = rotaT3(coor,connex,iel);
% matrices élémentaires
kel = kelem_Trian_true_drill_rot(xe);    

%AAAA=eig(kel);

keg=rog'*kel*rog;
% assemblage
KG = assemblage(KG,connex,keg,nne,ddln,iel);
end

% construction du vecteur forces
[FG]= charge (FG,charg,ddln);

KGA=KG;
% Application des CL:   numnoued  ddl1   ddl2]
% L: libre            cl=[1         L    0.03] 
[KG,FG]=cond_lim(KG,FG,climit,NTDDL,ddln);

% Résolution du systeme
DG = KG\FG ;

% calcules les reactions
React=KGA*DG;

% somme de réactions et équilibre
SFx=0;
SFy=0;
for i=1:2:NTDDL-1
SFx=SFx+React(i);
SFy=SFy+React(i+1);
end

% effort internes (contraintes généralisées)
for iel = 1:NTE

% matrice de rotation
[rog,xe] = rotaT3(coor,connex,iel);
[depl] = local_disp(DG,connex,rog,iel);
% calcul des contraintes
[forcm] = efforts_elem(xe,depl); % sigmax, sigmay, ..... sigmaxy

% ********  contraintes dans le repere global   ********
ro=rog(1:2,1:2);
[forcmg] = stressg(forcm,ro);% sigmaX, sigmaY, ..... sigmaXY

n1=(iel-1)*nne+1;
n2=iel*nne;
efforts(n1:n2,:)=forcmg;

end
%----------------------------------------------------
%----------------------------------------------------
%1 sigmaX,   2 sigmaY,   3 sigmaXY
plot_stress_T3(coor,connex,NTE,efforts(:,1));


aaa=0;