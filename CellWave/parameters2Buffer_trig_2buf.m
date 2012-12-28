function params = parameters2Buffer
%
%  set  parameters for Li-Rinzel 2 Buffer model
%

% ..dimensional parameters, LiRinzel/2Buffer notes, table 1

params.nu_L       = 0.01;  % nondimensional
params.nu_c       = 500.;      % nondimensional
params.d_act      = 0.7;   % micro M
params.d_inh      = 0.6;   % micro M
params.d_I        = 0.025; % micro M
params.V_m        = 1.;   % micro M
params.k_p        = 0.25;   % micro M
params.C_er       = 10.;   % micro M
params.k_on       = 8.0;   % 1/ (micro M * sec)

params.gamma      = 0.;    % 1/sec
params.b1tot      = 200.;  % micro M
params.b2tot      = 1.;    % micro M
params.K1         = 10.;   % micro M
params.K2         = 0.24;  % micro M
params.diffCalcium= 300.;  % micro m^2/sec
params.H          = 2.;    % micro M/sec

%diffCalcium =0.; %TRY THIS

params.c_0 = 0.115;
params.p_0 = 0.96; % IP3 concentration
params.h_0 = 0.93;

