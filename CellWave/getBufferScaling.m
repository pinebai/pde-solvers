function theta = getBufferScaling( params, c )
%function theta = getBufferScaling( params, c )
%   computes scaling factor for rapid buffering approx

q=params;

theta= q.b1tot*q.K1*(q.K1 + c).^(-2) ...
       +q.b2tot*q.K2*(q.K2 + c).^(-2);
