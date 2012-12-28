function [caFlux, J_flux, J_pump]=ip3rOpenProb( params, p, c )
%function [caFlux, J_flux, J_pump]=ip3rOpenProb( params, p, c )
%   compute IP3 receptor 'open probability' as a fcn of [Ca2+]
%   .. assumes RBA, which is valid for near-equil. Ca2+ buffers
%

nu_L   =params.nu_L;
nu_c   =params.nu_c ;
d_act  =params.d_act ;
d_inh  =params.d_inh;
d_I    =params.d_I  ;
V_m    =params.V_m ;
k_p    =params.k_p  ;
C_er   =params.C_er ;
k_on   =params.k_on ;

b1tot =params.b1tot ; 
b2tot =params.b2tot ; 
K1    =params.K1    ; 
K2    =params.K2    ; 

q= nu_c*( p./( p+ d_I ));

%----------------------evaluate rhs
h= d_inh./(d_inh+c);

theta = getBufferScaling( params, c );

J_flux = ( ( nu_L +q.*(c./(c+d_act)).^3.*h.^3).*(C_er-c)) ./(1+theta);
%%J_flux = ( ( nu_L +q.*(c./(c+d_act)).^3.*h.^3).*C_er) ./(1+theta);
J_pump =  V_m*(c.^2./(c.^2 + k_p.^2))./ (1 +theta);

caFlux = J_flux-J_pump;
