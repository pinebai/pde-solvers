function dydt = careaction1( t, y )
global pglobal parglobal 

params = parglobal;
nu_L = params.nu_L;
d_I  = params.d_I ;
k_p  = params.k_p;
nu_P = params.nu_P  ;
I_s  = params.I_s  ;
C_er = params.C_er  ;
tau_0= params.tau_0 ;
d_act= params.d_act ;
d_inh= params.d_inh;
beta = params.beta;
lambda= params.lambda;

p = pglobal;
q= ( p./( p+ d_I ));

%----------------------evaluate rhs
c=y(1,:);
h=y(2,:);
J_flux = ( nu_L +q.*(c./(c+d_act)).^3.*h.^3).*(C_er-c);
J_pump =  nu_P*(c.^2./(c.^2 + k_p.^2));
dct = beta*lambda*(J_flux -J_pump);
%dht = (d_inh - (d_inh + c).*h)/(d_inh*tau_0);
dht = (d_inh./(d_inh + c) -h)/tau_0;

dct = dct(:);
dht = dht(:);
%%size(dct)

%qip3=0.925774*(1-t/350);
%%qip3=1.5*(1-t/350);
%if(t>150),
%  qip3=0.8;
%end

%.. Wagner et al version
%%dct = beta*lambda*(J_flux -J_pump);
%%dht = (d_inh./(d_inh + c) - h)/tau_0;
dydt = [dct, dht]';
