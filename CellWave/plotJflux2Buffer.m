%
% ..compute Jflux from the 2Buffer model
%
ibuffer=1; % use buffering
params=dbg_2buffer_fluxbc6; params.b1tot=200; params.b2tot=10;
p=0.002;
checkJflux2Buffer
Jflux_p0002=Jflux;

p=0.02;
checkJflux2Buffer
Jflux_p002=Jflux;

p=1.;
checkJflux2Buffer
Jflux_p1=Jflux;

loglog(c,Jflux_p0002, c,Jflux_p002, c,Jflux_p1, ...
	  c,Jpump );
legend('Release; p=0.002 \muM', 'Release; p=0.02 \muM', 'Release; p=1 \muM', 'SERCA')
grid
xlabel('[ Ca^{2+} ] concentration  (\muM)')
ylabel('[ Ca^{2+} ] efflux (\muM/s)')
title('Ca^{2+} release flux from IP_3 dependent receptors (from Li-Rinzel w/ 2 Buffers)');

