%
% checking J_flux from Sneyd-Keener's book, Fig. 5.9 (params from Table 12.2, p.352)
%
sk.kf      =  1.4;    % uM/s
%sk.kf      = 1;    % uM/s
sk.mu0     = 0.567;   
sk.mu1     = 0.433;  
sk.kmu     = 4.;     % uM
sk.b       = 0.11;
sk.k1      = 0.7;    % uM
sk.k2      = 0.7;    % uM
sk.kgamma  = 0.1;    % uM
sk.beta    = 0.02;   % uM/s
sk.gamma   = 2;

h= sk.k2^2*(sk.k2^2 + c.^2 ).^(-2);

Jflux = sk.beta+ sk.kf*(sk.mu0+ sk.mu1*p/(p+sk.kmu)).*(sk.b + (1-sk.b)*c./(sk.k1 + c)).*h;
Jpump = sk.gamma*c./( sk.kgamma + c );

%.. model from p. 351 in Keener & Sneyd
%mu_p  = p^3/(sk.kmu^3 + p^3);
%Jflux = sk.kf*(sk.mu0+ sk.mu1*p/(p+sk.kmu)).*(sk.b + (1-sk.b)*c./(sk.k1 + c)).*h;

