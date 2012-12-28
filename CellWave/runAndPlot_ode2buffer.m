%% VIZ code for looking at results from 
outfile='out_ref2';
params = ref_hiro_1;

close all
nullClines2Buffer(params, params.ip3_0)
figure(2)
[tcomp,q,qAll]=odeLoad2Buffer( 'OUT',outfile);

plot2Buffer(tcomp, q);

