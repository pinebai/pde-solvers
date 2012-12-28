function caNullClines( params, p )
%
% Compute null clines for the Wagner-Li-Keizer system
% -- dimensional version
%
global parglobal  pglobal
nu_L = params.nu_L;
d_I  = params.d_I ;
k_p  = params.k_p;
nu_P = params.nu_P  ;
I_s  = params.I_s  ;
C_er = params.C_er  ;
tau_0= params.tau_0 ;
d_act= params.d_act ;
d_inh= params.d_inh;
beta = params.beta;
lambda= params.lambda;
 
q= qGet( params, p);
pglobal = p;
parglobal=params;
tcomp=0.;

ca=0.02:0.1:5;
h = 0.02:0.1:1; 
%plotbox = [ 0,2, 0, 1.2]; % [ cmin, cmax, hmin, hmax]
plotbox = [ 0,3, 0.1, 1.1]; % [ cmin, cmax, hmin, hmax]
%plotbox = [ 0.1, 1.1, 0, 7]; % [ hmin, hmax, cmin, cmax]

%1)  Null cline from dh/dt=0 --> 0=g
%ca1=0.02:0.1:5;
powers=log10(0.01):0.05:1;
ca1=10.^(powers);
hcline1 = d_inh./(d_inh+ca1);

dz1 = careaction1(tcomp, [ca1; hcline1]);
dc1 = dz1(1,:);
dh1 = dz1(2,:);

%2)  Null cline from dc/dt=0 --> 0=f
powers=log10(0.01):0.003:.9;
ca2=10.^(powers);
zz=(ca2.^2 + k_p.^2).*(C_er - ca2);
hcubed  = (1./q).*( (ca2+d_act)./ca2 ).^3.* ...
           (-nu_L + nu_P*ca2.^2./zz);
hcline2 = sign(hcubed).*abs(hcubed).^(1/3);
%dh2= (d_inh - (d_inh+ca2))/tau_0;
%dc2= beta*lambda*(( nu_L +q*(ca2./(ca2+d_act)).^3.*hcline2.^3).*(C_er-ca2) ...
%- nu_P.*(ca2.^2./(ca2.^2 + k_p.^2)));
dz2 = careaction1(tcomp, [ca2; hcline2]);
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
title(['Null clines ,[IP_3]=',num2str(p)])
