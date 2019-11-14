% debug sequence
indb = 20; % autarky

figure(1)
plot(DebChoiceGrid,v1wtemp);
title('Aut: v1w')
saveas(gcf,'./debug/aut_v1w.emf','emf');

figure(2)
plot(DebChoiceGrid,v1etemp);
title('Aut: v1e')
saveas(gcf,'./debug/aut_v1e.emf','emf');

figure(3)
plot(DebChoiceGrid,val1temp);
title('Aut: val1')
saveas(gcf,'./debug/aut_val1.emf','emf');

figure(4)
plot(DebChoiceGrid,r1temp);
title('Aut: r1')
saveas(gcf,'./debug/aut_r1.emf','emf');

figure(5)
plot(DebChoiceGrid,p1temp);
title('Aut: p1')
saveas(gcf,'./debug/aut_p1.emf','emf');

figure(6)
plot(DebChoiceGrid,phi1temp);
title('Aut: phi1')
saveas(gcf,'./debug/aut_phi1.emf','emf');

figure(7)
plot(DebChoiceGrid,v1wprtemp);
title('Aut: v1wpr')
saveas(gcf,'./debug/aut_v1wpr.emf','emf');

% =========================================================================
indb = 14; % fin

figure(1)
plot(DebChoiceGrid,v1wtemp);
title('Fin: v1w')
saveas(gcf,'./debug/fin_v1w_cmpr_dpos.emf','emf');

figure(2)
plot(DebChoiceGrid,v1etemp);
title('Fin: v1e')
saveas(gcf,'./debug/fin_v1e_cmpr_dpos.emf','emf');

figure(3)
plot(DebChoiceGrid,val1temp);
title('Fin: val1')
saveas(gcf,'./debug/fin_val1_cmpr_dpos.emf','emf');

figure(4)
plot(DebChoiceGrid,r1temp);
title('Fin: r1')
saveas(gcf,'./debug/fin_r1_cmpr_dpos.emf','emf');

figure(5)
plot(DebChoiceGrid,p1temp);
title('Fin: p1')
saveas(gcf,'./debug/fin_p1_cmpr_dpos.emf','emf');

figure(6)
plot(DebChoiceGrid,phi1temp);
title('Fin: phi1')
saveas(gcf,'./debug/fin_phi1_cmpr_dpos.emf','emf');

figure(7)
plot(DebChoiceGrid,v1wprtemp);
title('Fin: v1wpr')
saveas(gcf,'./debug/fin_v1wpr_cmpr_dpos.emf','emf');

figure(8)
plot(DebChoiceGrid,log(c1wtemp));
title('Fin: c1w')
saveas(gcf,'./debug/fin_c1w_cmpr_dpos.emf','emf');

% =========================================================================
pt = floor(nt/2);
pt = 5;
plot(bVec,val1_seqa(pt,:));
plot(bVec,v1w_seqa(pt,:));
plot(bVec,v1e_seqa(pt,:));
plot(bVec,Pri1_seqa(pt,:));
plot(bVec,R1_seqa(pt,:));
plot(bVec,phi1_seqa(pt,:));
plot(bVec,DebPol1_seqa(pt,:));


% =========================================================================
load('./debug/dzero.mat')
clearvars -except DebChoiceGrid c1wtemp r1temp ;
c1wtemp_dz = c1wtemp;
r1temp_dz = r1temp;
load('./debug/dpos.mat')
clearvars -except DebChoiceGrid c1wtemp r1temp c1wtemp_dz r1temp_dz;

wbar = 0.8;
wgt1 = 0.796;
b1 = 0.5154;
c1zero = wbar + wgt1*DebChoiceGrid./r1temp_dz' - wgt1*b1;
c1pos_cons = wbar + wgt1*DebChoiceGrid./(r1temp_dz'-0.22) - wgt1*b1;
c1pos = wbar + wgt1*DebChoiceGrid./(r1temp') - wgt1*b1;

figure(1)
hold on;
plot(DebChoiceGrid,c1wtemp_dz);
plot(DebChoiceGrid,c1wtemp,'r');
legend('D=0','D>0');
title('consumption')
hold off;

figure(2)
hold on;
plot(DebChoiceGrid,r1temp_dz./r1temp);
title('R(D=0)/R(D>0)')
hold off;

figure(3)
hold on;
plot(DebChoiceGrid,r1temp_dz-r1temp);
title('R(D=0)-R(D>0)')
hold off;

figure(4)
hold on;
plot(DebChoiceGrid,r1temp_dz);
plot(DebChoiceGrid,r1temp,'r');
legend('D=0','D>0');
title('interest')
hold off;

figure(5)
hold on;
plot(DebChoiceGrid,c1zero);
plot(DebChoiceGrid,c1pos_cons,'r');
plot(DebChoiceGrid,c1pos,'g');
legend('D=0','D>0, shift','D>0, original')
hold off;

figure(6)
hold on;
plot(DebChoiceGrid,r1temp_dz-0.2245);
plot(DebChoiceGrid,r1temp,'r');
legend('D=0, shift','D>0');
title('interest')
hold off;

c1pos_cons1 = wbar + wgt1*DebChoiceGrid./(r1temp_dz'-0.20) - wgt1*b1;
c1pos_cons2 = wbar + wgt1*DebChoiceGrid./(r1temp_dz'-0.15) - wgt1*b1;
c1pos_cons3 = wbar + wgt1*DebChoiceGrid./(r1temp_dz'-0.10) - wgt1*b1;
c1pos_cons4 = wbar + wgt1*DebChoiceGrid./(r1temp_dz'-0.05) - wgt1*b1;

figure(7)
hold on;
plot(DebChoiceGrid,c1pos_cons);
plot(DebChoiceGrid,c1pos_cons1);
plot(DebChoiceGrid,c1pos_cons2);
plot(DebChoiceGrid,c1pos_cons3);
plot(DebChoiceGrid,c1pos_cons4);
plot(DebChoiceGrid,c1zero);
legend('0.2245','0.20','0.15','0.10','0.05','dzero');
title('interest')
hold off;

% smooth policy function
DebPol_coef = polyfit(bVec',DebPol1,2);
DebPol_smooth = DebPol_coef(1)*bVec.^2 + DebPol_coef(2)*bVec ...
    + DebPol_coef(3);
plot(bVec,DebPol_smooth)

% steady state debt
figure(1)
hold on;
plot(bVec,bVec);
plot(bVec,DebPol1_seqa(25,:));
hold off;

% forward simulation
figure(2)
plot(bsim)

