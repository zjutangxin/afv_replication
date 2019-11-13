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
% clc;
clear all;
tic;
% Parameters
theta = 0.2;
zbar = 1.0;
zmin = 0.86;
zmax = 2.0*zbar - zmin ;
a_min = 1e-9;

dforeign = 0.10 ;

wbar = (1-theta)*zbar^theta ;
wgt1 = 0.796 ;
surv_rate = 0.975;
rbar = 0.03 ;
bbeta = surv_rate/(1+rbar);
nt = 20 ;

% initialization
nb = 40;
maxgrid = 500;
nz = 21;

% bmin = 0.0;
bmin = dforeign/wgt1 + 1e-9 ;
bmax = bmin+0.8;
bVec = linspace(bmin,bmax,nb);
DebChoiceGrid = linspace(bmin,bmax,maxgrid);

zvec = linspace(zmin,zmax,nz);
pzvec = (1.0/nz)*ones(nz,1) ;
pzvec = pzvec/sum(pzvec);

% variables
% initialized for nt+1
v1wMx = zeros(nb,1);v1eMx = zeros(nb,1);
val1Mx = zeros(nb,1);
Pri1Mx = zeros(nb,1);
R1Mx = zeros(nb,1);
phi1Mx = zeros(nb,1);
DebPol1 = zeros(nb,1);

vv1wMx = zeros(nb,1);
vv1eMx = zeros(nb,1);
vval1Mx = zeros(nb,1);
PPri1Mx = zeros(nb,1);
RR1Mx = zeros(nb,1);
pphi1Mx = zeros(nb,1);
DDebPol1 = zeros(nb,1);

v1w_seqa = zeros(nt+1,nb);
v1e_seqa = zeros(nt+1,nb);
val1_seqa = zeros(nt+1,nb);
Pri1_seqa = zeros(nt+1,nb);
R1_seqa = zeros(nt+1,nb);
phi1_seqa = zeros(nt+1,nb);
DebPol1_seqa = zeros(nt,nb);

% assign value for nt+1
v1w_seqa(nt+1,:) = v1wMx ;
v1e_seqa(nt+1,:) = v1eMx ;
val1_seqa(nt+1,:) = val1Mx ;
Pri1_seqa(nt+1,:) = Pri1Mx ;
R1_seqa(nt+1,:) = R1Mx ;
phi1_seqa(nt+1,:) = phi1Mx ;

% Terminal period
indt = nt ;
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
   ahat = max(ahat,a_min) ;
   EU1_e = log(1-eta_t) ;
   EU1_e = EU1_e + log(ahat)*pzvec ;
   c1_w = wbar + wgt1*(b1pr/R1-b1);
   U1_w = log(max(c1_w,a_min));
%    v1wpr = interp1(bVec,v1wMx,b1pr);
%    v1epr = interp1(bVec,v1eMx,b1pr);
   v1wpr = 0;
   v1epr = 0;
   v1w = U1_w + bbeta*v1wpr;
   v1e = EU1_e + bbeta*v1epr;
   
   % update value for current period
   vv1wMx(indb,1) = v1w ;
   vv1eMx(indb,1) = v1e ;
   PPri1Mx(indb,1) = p1 ;
   RR1Mx(indb,1) = R1 ;
   pphi1Mx(indb,1) = phi1 ;
   vval1Mx(indb,1) = wgt1*v1w+(1-wgt1)*v1e;
end

v1wMx = vv1wMx ;
v1eMx = vv1eMx ;
val1Mx = vval1Mx ;
Pri1Mx = PPri1Mx ;
R1Mx = RR1Mx ;
phi1Mx = pphi1Mx ;

v1w_seqa(indt,:) = v1wMx ;
v1e_seqa(indt,:) = v1eMx ;
val1_seqa(indt,:) = val1Mx ;
Pri1_seqa(indt,:) = Pri1Mx ;
R1_seqa(indt,:) = R1Mx ;
phi1_seqa(indt,:) = phi1Mx ;

% disp('terminal period')

% Period t
for indt = nt-1:-1:1
%     if indt == 11
%         pause;
%     end
   % compute the optimal policy
   for indb = 1:1:nb
%        v1wtemp = zeros(maxgrid,1);
%        v1etemp = zeros(maxgrid,1);
%        r1temp = zeros(maxgrid,1);
%        c1wtemp = zeros(maxgrid,1);
%        v1wprtemp = zeros(maxgrid,1);
%        p1temp = zeros(maxgrid,1);
%        phi1temp = zeros(maxgrid,1);

       val1temp = zeros(maxgrid,1);
       b1 = bVec(indb) ;
       for indbp = 1:1:maxgrid
          b1pr = DebChoiceGrid(indbp) ; 
          be1 = wgt1*b1 - dforeign ;
          be1pr = wgt1*b1pr - dforeign ;
          
          alfa_t = (1-bbeta^(nt-indt+1))/(1-bbeta);
          eta_t = (bbeta-bbeta^(nt-indt+1))/(1-bbeta^(nt-indt+1)) ;
          
          p1pr = interp1(bVec,Pri1Mx,b1pr);
          afun = theta*zvec/(zbar^(1-theta))+(zvec-zbar)*p1pr/zbar;
          phi1 = ((afun+p1pr)./(afun+p1pr+be1pr))*pzvec ;
          afunbar = theta*zbar/(zbar^(1-theta));
          p1 = eta_t*phi1*(afunbar+be1)/(1-eta_t*phi1);
          
          R1 = (1-eta_t*phi1)*be1pr/(eta_t*(1-phi1)*(afunbar+be1));
          ahat = afun + p1 + be1 ;
          ahat = max(ahat,a_min) ;
          EU1_e = log(1-eta_t)+(alfa_t-1)*log(eta_t*phi1/p1);
          EU1_e = EU1_e + alfa_t*log(ahat)*pzvec ;
          
          c1_w = wbar + wgt1*(b1pr/R1-b1);
          U1_w = log(max(c1_w,a_min));
          
          v1wpr = interp1(bVec,v1wMx,b1pr);
          v1epr = interp1(bVec,v1eMx,b1pr);
          
          v1w = U1_w + bbeta*v1wpr;
          v1e = EU1_e + bbeta*v1epr;
          
%           r1temp(indbp,1) = R1;
%           p1temp(indbp,1) = p1;
%           c1wtemp(indbp,1) = c1_w;
%           v1wprtemp(indbp,1) = v1wpr;
%           v1wtemp(indbp,1) = v1w;
%           v1etemp(indbp,1) = v1e;
%           phi1temp(indbp,1) = phi1;

          val1temp(indbp,1) = wgt1*v1w + (1-wgt1)*v1e ;
       end
       % find optimal policy
       [xtemp xoptim] = max(val1temp);
       DDebPol1(indb) = DebChoiceGrid(xoptim) ;
   end
   
   % check convergence of policy
   disp(['policy iter = ',num2str(nt-indt),', error = ', ...
       num2str(sum(abs(DebPol1-DDebPol1))/nb)])
   
   DebPol1 = DDebPol1 ;
   DebPol1_seqa(indt,:) = DebPol1 ;
   
   % Compute the equilibrium at the optimal policy
   for indb = 1:1:nb
       b1 = bVec(indb) ;
       b1pr = DebPol1(indb) ;
       
       be1 = wgt1*b1 - dforeign ;
       be1pr = wgt1*b1pr - dforeign ;
          
       alfa_t = (1-bbeta^(nt-indt+1))/(1-bbeta);
       eta_t = (bbeta-bbeta^(nt-indt+1))/(1-bbeta^(nt-indt+1)) ;
       
       p1pr = interp1(bVec,Pri1Mx,b1pr);
       afun = theta*zvec/(zbar^(1-theta))+(zvec-zbar)*p1pr/zbar;
       phi1 = ((afun+p1pr)./(afun+p1pr+be1pr))*pzvec ;
       afunbar = theta*zbar/(zbar^(1-theta));
       p1 = eta_t*phi1*(afunbar+be1)/(1-eta_t*phi1);
          
       R1 = (1-eta_t*phi1)*be1pr/(eta_t*(1-phi1)*(afunbar+be1));
       ahat = afun + p1 + be1 ;
       ahat = max(ahat,a_min) ;
       EU1_e = log(1-eta_t)+(alfa_t-1)*log(eta_t*phi1/p1);
       EU1_e = EU1_e + alfa_t*log(ahat)*pzvec ;
          
       c1_w = wbar + wgt1*(b1pr/R1-b1);
       U1_w = log(max(c1_w,a_min));
          
       v1wpr = interp1(bVec,v1wMx,b1pr);
       v1epr = interp1(bVec,v1eMx,b1pr);
          
       v1w = U1_w + bbeta*v1wpr;
       v1e = EU1_e + bbeta*v1epr;      
       
       % update value for current period
       vv1wMx(indb,1) = v1w ;
       vv1eMx(indb,1) = v1e ;
       PPri1Mx(indb,1) = p1 ;
       RR1Mx(indb,1) = R1 ;
       pphi1Mx(indb,1) = phi1 ;
       vval1Mx(indb,1) = wgt1*v1w+(1-wgt1)*v1e;       
   end
   
   v1wMx = vv1wMx ;
   v1eMx = vv1eMx ;
   val1Mx = vval1Mx ;
   Pri1Mx = PPri1Mx ;
   R1Mx = RR1Mx ;
   phi1Mx = pphi1Mx ;

   v1w_seqa(indt,:) = v1wMx ;
   v1e_seqa(indt,:) = v1eMx ;
   val1_seqa(indt,:) = val1Mx ;
   Pri1_seqa(indt,:) = Pri1Mx ;
   R1_seqa(indt,:) = R1Mx ;
   phi1_seqa(indt,:) = phi1Mx ;   
end

disp(['time = ', num2str(toc)]);

% save('./results/autarky.mat');
% save('./results/df_1e2.mat');






























