%
% ..compute Jflux from Sneyd-Keener
%
p=0.002;
checkJfluxSneydKeenerFig_5_9
Jflux_p002=Jflux;
p=2;
checkJfluxSneydKeenerFig_5_9
Jflux_p2=Jflux;
p=1;
checkJfluxSneydKeenerFig_5_9
Jflux_p1=Jflux;
p=0.2;
checkJfluxSneydKeenerFig_5_9
Jflux_p02 =Jflux;

semilogx(c,Jflux_p002, c,Jflux_p02, c,Jflux_p1), legend('p=0.002', 'p=0.2','p=1')
grid
xlabel('[ Ca^{2+} ] concentration  (\muM)')
ylabel('[ Ca^{2+} ] efflux (\muM/s)')
title('Ca^{2+} release flux from IP_3 dependent receptors (from Keener & Sneyd)')

