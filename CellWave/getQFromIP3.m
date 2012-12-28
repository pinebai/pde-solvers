function q = getQForIP3(params, p)
%
%  compute (nondim) q from a dimensional [IP3] concentration p
%

q= (p./( p + params.d_I )).^3;

