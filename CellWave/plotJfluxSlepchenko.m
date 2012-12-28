%
% ..compute Jflux from Slepchenko
%
slepchenkoParams, sparams.J0=75; checkJfluxSlepchenko
Jflux_J10=Jflux;
Jpump_J10=Jpump;
slepchenkoParams, sparams.J0=50; checkJfluxSlepchenko
Jflux_J1000=Jflux;
Jpump_J1000=Jpump;
slepchenkoParams, sparams.J0=100.; checkJfluxSlepchenko
Jflux_J100 =Jflux;
Jpump_J100 =Jpump;

semilogx(c,Jflux_J10, c,Jflux_J100, c,Jflux_J1000, ...
	  c,Jpump_J10 );
legend('Release; J_0=75', 'Release; J_0=100', 'Release; J_0=50', 'SERCA')
grid
xlabel('[ Ca^{2+} ] concentration  (\muM)')
ylabel('[ Ca^{2+} ] efflux (\muM/s)')
title('Ca^{2+} release flux from IP_3 dependent receptors (from Slepchenko et al.)')

