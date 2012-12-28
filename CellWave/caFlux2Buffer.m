%
%  Ca flux plot ala Sneyd & Keener, Fig. 5.9
%   but for my Li-Rinzel 2 Buffer model
%

%%--- parameters -- MUST DEFINE params BEFORE CALLING THIS
%% params = stuff 

params = ref_hiro_1;

ip3list = [0.25 0.5 1.0 2.0 ];
c       = 10.^[-2:0.1:1];  %in mu M, compare w. Fig. 5.9 which is in [M]

ifig=[1 2 3];
plotSubFigures=1;

if( ~plotSubFigures ),
  figure(ifig(1)), clf
  figure(ifig(2)), clf
  figure(ifig(3)), clf
else
  clf
end

for p =  ip3list

  ibuffer=0;
  checkJflux2Buffer;
  
  %openProb = p1.*p2.*p3;
  %Jflux    = kf* (p1.*p2.*p3);

  %%-------------------FIGURE 1-------------------------
  figNumber=1;
  if(plotSubFigures),
    eval(['subplot(31',int2str(figNumber),');']);
  else
    figure(figNumber);
  end
  semilogx( c, openProb); 
  title(['(LR 2BUFF) Open probability of [Ca2+] channels, [IP3]=',num2str(p),' \mu M']);
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
 
  semilogx( c, 100*Jflux/max(Jflux),'b', c,100*(Jflux - Jpump)/max(Jflux-Jpump),'r-.'); 
  title(['(LR 2BUFF) Percentage of full [Ca2+] efflux, [IP3]=',num2str(p),' \mu M']);
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
  semilogx( c, Jflux, 'b', c, Jflux - Jpump, 'r-.'   );
  title(['(LR 2BUFF) [Ca2+] flux vs. calcium concentration, [IP3]=',num2str(p),' \mu M']);
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
