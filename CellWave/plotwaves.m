function plotwaves( ndim, ic, istep,q )
%function plotwaves( ndim, ic, istep )
%  ndim  = number of points on r-axis
%  ic    = component to plot, we assume ic=0 is r, and there are 5 components
%  istep = which timestep to plot
%  q     = contents of radialSolverSlepchenko2Buffer output file
%
%  ic: 0=r, 1=c, 2=h, 3=p, 4=b1, 5=b2

ncomp=6;
ii=3:2+ndim;
rline = 1+ncomp*istep;
icline= rline+ic;

r=q(rline,ii);
cc=q(icline,ii);

name=['r','c','h','p','b','b'];
cname=name(ic+1);
if(ic>=4),
  cname=[cname,num2str(ic-3)];
end

fprintf(1,'component %d (%s)', [q(icline,1),cname]);
fprintf(1,' at t=%f\n',q(icline,2));
plot( r, cc );
