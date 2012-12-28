function [b1, b2] = computeBuffers( params, c )
%
%  Compute equilibrium values for the Ca2+ buffers
%

b1= params.b1tot*c./(params.K1 + c);
b2= params.b2tot*c./(params.K2 + c);


