%
% Slepchenko params, table II
%

%iopt =1; % high excitability Table II
%iopt =2; % low  excitability Table II
iopt =3; % X. Laevis params, Table III

%ibuffer =0;  % no buffers
ibuffer =1;   % use buffers

if (iopt==1),
  sparams.J0      =  100;
  sparams.d_act   =  0.7;
  sparams.Vm      =  1.0;
  sparams.Kp      =  0.25;
  sparams.leak    =  0.0151;
  sparams.k_on    =  2.0; % 1/(uM s)
  sparams.d_inh   =  0.6; % uM

  sparams.b1tot   =  200;    % low affinity buffer,  K=  10    uM
  sparams.K1      =  10;

  sparams.b2tot   =  10;     % high affinity buffer, K=   0.25 uM
  sparams.K2      =  0.24;      

elseif (iopt==2),
  sparams.J0      =  110;
  sparams.d_act   =  0.7;
  sparams.Vm      =  1.0;
  sparams.Kp      =  0.1;
  sparams.k_on    =  2.0; % 1/(uM s)
  sparams.d_inh   =  0.6; % uM
  sparams.leak    =  0.1755;

  sparams.b1tot   =  200;    % low affinity buffer,  K=  10    uM
  sparams.K1      =  10;

  sparams.b2tot   =  10;     % high affinity buffer, K=   0.25 uM
  sparams.K2      =  0.24;      

elseif (iopt==3),
  sparams.J0      =  1000;
  sparams.d_act   =  0.7;
  sparams.Vm      =  1.0;
  sparams.Kp      =  0.25;
%  sparams.leak    =  0.0151;
  sparams.leak    =  0.11;

  sparams.k_on    =  2.0; % 1/(uM s)
  sparams.d_inh   =  0.6; % uM

  sparams.b1tot   =  200;    % low affinity buffer,  K=  10    uM
  sparams.K1      =  10;

  sparams.b2tot   =  10;     % high affinity buffer, K=   0.25 uM
  sparams.K2      =  0.24;      

end
