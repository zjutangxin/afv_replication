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
saveas(gcf,'./debug/fin_v1w.emf','emf');

figure(2)
plot(DebChoiceGrid,v1etemp);
title('Fin: v1e')
saveas(gcf,'./debug/fin_v1e.emf','emf');

figure(3)
plot(DebChoiceGrid,val1temp);
title('Fin: val1')
saveas(gcf,'./debug/fin_val1.emf','emf');

figure(4)
plot(DebChoiceGrid,r1temp);
title('Fin: r1')
saveas(gcf,'./debug/fin_r1.emf','emf');

figure(5)
plot(DebChoiceGrid,p1temp);
title('Fin: p1')
saveas(gcf,'./debug/fin_p1.emf','emf');

figure(6)
plot(DebChoiceGrid,phi1temp);
title('Fin: phi1')
saveas(gcf,'./debug/fin_phi1.emf','emf');

figure(7)
plot(DebChoiceGrid,v1wprtemp);
title('Fin: v1wpr')
saveas(gcf,'./debug/fin_v1wpr.emf','emf');

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
