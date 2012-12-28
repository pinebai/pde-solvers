%
% assumes 'sparams' are defined as
%    sparams.J0     = IP3 effects
%    sparams.d_act  = ER channel activation
%    sparams.Vm     = SERCA pump rate
%    sparams.Kp     = SERCA pump dissociation constant 
%    sparams.leak   = Ca leak across cell membrane
%    sparams.K1     = dissociation const for buffer 1
%    sparams.K2     = dissociation const for buffer 2   
%    sparams.b1tot  = total buffer amount for buf 1
%    sparams.b2tot  = total buffer amount for buf 2
%
% assumes vectors c,h are defined
%
% returns J_flux, J_pump from Slepchenko, Schaff & Choi
%

sp= sparams;

h = sparams.d_inh * (sparams.d_inh + c).^(-1);

Ju_flux = (sp.leak + sp.J0*h.^3.*( c./(c+sp.d_act) ).^3);
Ju_pump = sp.Vm* ( c.^2 ./( c.^2 + sp.Kp.^2));

%.. Ju = unbuffered fluxes. These are reduced by buffering (&less available Ca)

if( ibuffer == 1),
    usingBuffers='';
    theta1= sp.b1tot*sp.K1*( c + sp.K1  ).^(-2);
    theta2= sp.b2tot*sp.K2*( c + sp.K2  ).^(-2);
else
    usingBuffers='not';
    theta1=1.;
    theta2=1.;
end
disp(['Computing Slepchenko Ca2+ flux,  ',usingBuffers, ...
       ' using buffers, J0 = ',num2str(sp.J0)]);

Jflux = Ju_flux./(1+theta1+theta2);
Jpump = Ju_pump./(1+theta1+theta2);


