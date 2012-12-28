%
%  Ca flux plot from Sneyd & Keener's "Heuristic model", chap. 5
%

%%--- parameters, Table 12.1 in S&K
mu0 =  0.567 ; % mu M
mu1 =  0.433 ; % mu M
kmu =  4.0   ; % mu M
k1  =  0.7   ; % mu M
k2  =  0.7   ; % mu M
b   =  0.11  ; % [nondim]
kf  =  8.1   ; % mu M/s

beta  = 0.02 ; % mu M/s  Ca2+ leak into the cell cytoplasm
gamma = 2.0  ; % mu M/s  pump large [Ca2+] strength
kgamma= 0.1  ; % mu M    pump K factor

ip3list = [0.25 0.5 1.0 2.0 ];
c       = 10.^[-2:0.1:1];  %in mu M, compare w. Fig. 5.9 which is in [M]

ifig=[4 5 6];
plotSubFigures=1;

if( ~plotSubFigures ),
  figure(ifig(1)), clf
  figure(ifig(2)), clf
  figure(ifig(3)), clf
else
  clf
end
for p =  ip3list
  
  p1= mu0 +   mu1*p./(p + kmu );
  p2= b   + (1-b)*c./(c + k1  );
  p3= k2.^2 ./ ( k2.^2 + c.^2 );

  openProb = p1.*p2.*p3;
  Jflux    = kf* (p1.*p2.*p3);
  Jpump    = gamma*c./(kgamma + c);
  
  totalFlux = Jflux - Jpump + beta;
  
  %%-------------------FIGURE 1-------------------------
  figNumber=1;
  if(plotSubFigures),
    eval(['subplot(31',int2str(figNumber),');']);
  else
    figure(figNumber);
  end
  semilogx( c, openProb); 
  title(['(HEURISTIC MODEL) Open probability of [Ca2+] channels, [IP3]=',num2str(p),' \mu M']);
  if( ~plotSubFigures),
    xlabel('[Ca2+] concentration (\mu M)');
  end
  ylabel('Open probability');
  if(plotSubFigures),
    grid on
  end
  hold on

  %%-------------------FIGURE 2-------------------------
  figNumber=2;
  if(plotSubFigures),
    eval(['subplot(31',int2str(figNumber),');']);
  else
    figure(figNumber);
  end
 
  semilogx( c, 100*Jflux/max(Jflux),'b',c,100*totalFlux/max(totalFlux),'r-.'); 
  title(['HEURISTIC MODEL) Percentage of full [Ca2+] efflux, [IP3]=',num2str(p),' \mu M']);
  if( ~plotSubFigures),
    xlabel('[Ca2+] concentration (\mu M)');
  end
  ylabel('[Ca2+] efflux (%)');
  axis([1e-2,1e1, -150, 100])
  if(plotSubFigures),
    grid on
  end
  hold on

  %%-------------------FIGURE 3-------------------------
  figNumber=3;
  if(plotSubFigures),
    eval(['subplot(31',int2str(figNumber),');']);
  else
    figure(figNumber);
  end

  semilogx( c, Jflux, 'b', c, totalFlux, 'r-.' );
  title(['(HEURISTIC MODEL) [Ca2+] flux vs. calcium concentration, [IP3]=',num2str(p),' \mu M']);
  xlabel('[Ca2+] concentration (\mu M)')
  ylabel('Ca^{2+} efflux (\mu M/sec)');
  if(plotSubFigures),
    grid on
  end
  hold on

  disp(['Showing [IP3]=',num2str(p),' \mu M;  Press enter for next case.']);
  pause
end

if( ~plotSubFigures),
  figure(ifig1), grid on,hold off
  figure(ifig2), grid on,hold off
  figure(ifig3), grid on,hold off
end
