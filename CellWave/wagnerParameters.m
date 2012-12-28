function params = wagnerParameters
%
%  set  parameters as in Wagner
%

% ..dimensional parameters, from Table 1

params.nu_L       = 5e-4;  % nondimensional
params.d_I        = 0.025; % micro M
params.k_p        = 0.4;   % micro M
params.nu_P       = 0.1;   % micro M
params.I_s        = 0.12;  % micro M
params.C_er       = 10.;   % micro M
params.tau_0      = 4;     % sec
%params.tau_0      = 25;     % sec
params.d_act      = 1.2; %1.2;   % micro M
params.diffCalcium= 300.;  % micro m^2/sec
%diffCalcium =0.; %TRY THIS
params.lambda     = 112.5; % 1/sec
params.d_inh      = 1.5;   % micro M
params.beta       = 0.053; % nondimensional

params.c_0 = 0.115;
params.p_0 = 0.96; % IP3 concentration
params.h_0 = 0.93;
