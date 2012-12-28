function dydt = ca2BufferReaction( t, y )

global pglobal parglobal 
params = parglobal;

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

p = pglobal;
q= nu_c*( p./( p+ d_I )).^3;

%----------------------evaluate rhs
c=y(1,:);
h=y(2,:);
J_flux = ( nu_L +q.*(c./(c+d_act)).^3.*h.^3).*(C_er-c);
J_pump =  V_m*(c.^2./(c.^2 + k_p.^2));
%theta= 1+ b1tot*c./(c+K1)  + b2tot*c./(c+K2);
theta= 1+ b1tot*K1*(c+K1).^(-2)  + b2tot*K2*(c+K2).^(-2);

dct = (J_flux -J_pump)./theta;
dht = k_on*(d_inh - (d_inh + c).*h);

dct = dct(:);
dht = dht(:);
%%size(dct)

%qip3=0.925774*(1-t/350);
%%qip3=1.5*(1-t/350);
%if(t>150),
%  qip3=0.8;
%end

dydt = [dct, dht]';
