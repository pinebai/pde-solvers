%
% computes the ca2+ RHS with 2 buffers
% ASSUMES:
%    params contains the valid 2Buffer parameters
%    c,h are set         (=vectors)
%    p   has been chosen (=parameter)
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
d_inh =params.d_inh ;

q= nu_c*( p./( p+ d_I )).^3;

%----------------------evaluate rhs

h = d_inh * ( d_inh + c).^(-1);

Ju_flux = ( nu_L +q.*(c./(c+d_act)).^3.*h.^3).*(C_er-c);
Ju_pump =  V_m*(c.^2./(c.^2 + k_p.^2));
theta= 1+ b1tot*c./(c+K1)  + b2tot*c./(c+K2);

openProb_2buffer = q.*(c./(c+d_act)).^3.*h.^3/nu_c;

J0 = q.*(C_er - c);
J0Min= min(J0);
J0Max= max(J0);

%.. Ju = unbuffered fluxes. These are reduced by buffering (&less available Ca)

if( ibuffer == 1),
    usingBuffers='';
    theta1= b1tot*K1*( c + K1  ).^(-2);
    theta2= b2tot*K2*( c + K2  ).^(-2);
else
    usingBuffers='not';
    theta1=1.;
    theta2=1.;
end

disp(['Computing Li-Rinzel/2Buffer  Ca2+ flux,  ',usingBuffers, ...
       ' using buffers, ',num2str(J0Min),' < J0 < ',num2str(J0Max) ]);

Jflux = Ju_flux./(1+theta1+theta2);
Jpump = Ju_pump./(1+theta1+theta2);


