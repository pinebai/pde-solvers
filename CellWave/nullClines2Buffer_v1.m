function caNullClines_v1( params, p )
%
% Compute null clines for the Li-Rinzel 2 buffer system
% -- dimensional version... use parameters2Buffer_v1
%      --> shows nice 
global parglobal  pglobal
 
nu_L   =params.nu_L;
nu_c   =params.nu_c ;
d_act  =params.d_act ;
d_inh  =params.d_inh;
d_I    =params.d_I  ;
nu_m   =params.nu_m ;
k_p    =params.k_p  ;
C_er   =params.C_er ;
k_on   =params.k_on ;

q= nu_c*getQFromIP3( params, p);
J_0 = q
pglobal = p;
parglobal=params;
tcomp=0.;

ca=0.02:0.1:5;
h = 0.02:0.1:1; 

plotbox = [ 0, 1.4, 0.2, 1.0]; % [ cmin, cmax, hmin, hmax]


%1)  Null cline from dh/dt=0 --> 0=g
%ca1=0.02:0.1:5;
powers=log10(0.01):0.05:1;
ca1=10.^(powers);
hcline1 = d_inh./(d_inh+ca1);

dz1 = ca2BufferReaction(tcomp, [ca1; hcline1]);
dc1 = dz1(1,:);
dh1 = dz1(2,:);

%2)  Null cline from dc/dt=0 --> 0=f
powers=log10(0.01):0.003:.9;
ca2=10.^(powers);
zz=(ca2.^2 + k_p.^2).*(C_er - ca2);
%%zz=(ca2.^2 + k_p.^2);
hcubed  = (1./q).*( (ca2+d_act)./ca2 ).^3.* ...
           (-nu_L + nu_m*ca2.^2./zz);
hcline2 = sign(hcubed).*abs(hcubed).^(1/3);
%dh2= (d_inh - (d_inh+ca2))/tau_0;
%dc2= beta*lambda*(( nu_L +q*(ca2./(ca2+d_act)).^3.*hcline2.^3).*(C_er-ca2) ...
%- nu_P.*(ca2.^2./(ca2.^2 + k_p.^2)));
dz2 = ca2BufferReaction(tcomp, [ca2; hcline2]);
dc2 = dz2(1,:);
dh2 = dz2(2,:);


%figure(4)
%semilogx(ca1,hcline1,'bo-', ca2, hcline2,'ro-'), grid, legend('h_t=0','c_t=0');
%plot(hcline1,ca1,'bo-', hcline2,ca2,'ro-'), grid, legend('h_t=0','c_t=0'); 
xlabel('h'), ylabel('Ca^{2+}')  
plot(ca1,hcline1,'bo-', ca2, hcline2,'r-'), grid, legend('h_t=0','c_t=0');
%%hold on
%%nmax=size(ca2,2);
%%jj=1:4:nmax;
%quiver(ca1,hcline1,dc1,dh1,0,'k');
%quiver(ca2(jj),hcline2(jj),dc2(jj),dh2(jj),0,'g');
%%hold off
axis(plotbox);
xlabel('Ca^{2+}'), ylabel('h');
title(['Null clines, [IP_3]=',num2str(p)])
