function [J0, leak] = computeSlepchenkoParams( params, p, c)
%function [J0, leak] = computeSlepchenkoParams( params, p, c)
%  convert 2Buffer model parameters to Slepchenko et al's
%  lumped J0 & leak variables

q=params;

leak = q.nu_L*( q.C_er - c );
J0   = q.nu_c*( p/(p+q.d_I)).^3 * (q.C_er - c);

