% -------------------------------------------------------------------------
%                      Program Description
% -------------------------------------------------------------------------
%   
% Purpose:
%     - Main program for one country with fixed foreign demand
%     - Saving Glut Project
%  
% Author:
%     - Xin Tang @ International Monetary Fund
%  
% Record of Revisions:
%         Date:                 Description of Changes
%     ============        =================================
%      11/11/2019                 Original Version
% =========================================================================
clc;
clear all;

% Parameters
theta = 0.2;
zbar = 1.0;
zmin = 0.86;
zmax = 2.0*zbar - zmin ;
a_min = 1e-9;

wbar = (1-theta)*zbar^theta ;
wgt1 = 0.796 ;
surv_rate = 0.975;
rbar = 0.03 ;
bbeta = surv_rate/(1+rbar);
nt = 10 ;
dforeign = 0 ;

% initialization
nb = 20;
maxgrid = 200;
nz = 21;

bmin = 0.0;
bmax = 0.8;
bVec = linspace(bmin,bmax,nb);
DebChoiceGrid = linspace(bmin,bmax,maxgrid);

zvec = linspace(zmin,zmax,nz);
pzvec = (1.0/nz)*ones(nz,1) ;
pzvec = pzvec/sum(pzvec);

% variables
v1wMx = zeros(nb,1);
v1eMx = zeros(nb,1);
val1Mx = zeros(nb,1);
Pri1Mx = zeros(nb,1);
R1Mx = zeros(nb,1);
phi1Mx = zeros(nb,1);
DebPol1 = zeros(nb,1);

v1w_seqa = zeros(nt+1,nb);
v1e_seqa = zeros(nt+1,nb);
val1_seqa = zeros(nt+1,nb);
Pri1_seqa = zeros(nt+1,nb);
R1_seqa = zeros(nt+1,nb);
phi1_seqa = zeros(nt+1,nb);
DebPol1_seqa = zeros(nt,nb);

v1w_seqa(nt+1,:) = v1wMx ;
v1e_seqa(nt+1,:) = v1eMx ;
val1_seqa(nt+1,:) = val1Mx ;
Pri1_seqa(nt+1,:) = Pri1Mx ;
R1_seqa(nt+1,:) = R1Mx ;
phi1_seqa(nt+1,:) = phi1Mx ;

% Terminal period
indt = nt ;
DebPol1 = 0 ;
DebPol1_seqa(nt,:) = DebPol1 ;

for indb = 1:1:nb
   b1 = bVec(indb) ;
   b1pr = DebPol1(indb) ;
   % solve_system
   be1 = wgt1*b1 - dforeign ;
   be1pr = wgt1*b1pr ;
   
   alfa_t = (1-bbeta^(nt-indt+1))/(1-bbeta);
   eta_t = (bbeta-bbeta^(nt-indt+1))/(1-bbeta^(nt-indt+1)) ;
   
   p1pr = 0 ;
   afun = theta*zvec/(zbar^(1-theta))+(zvec-zbar)*p1pr/zbar;
   phi1 = ((afun+p1pr)./(afun+p1pr+be1pr))*pzvec ;
   afunbar = theta*zbar/(zbar^(1-theta));
   p1 = eta_t*phi1*(afunbar+be1)/(1-eta_t*phi1);
   R1 = 1.0 ;
   ahat = afun + p1 + be1 ;
   
   
end

































