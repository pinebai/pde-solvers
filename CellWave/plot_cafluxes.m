%% PLOTS the [Ca2+] efflux data for comparison w/ Sneyd & Keener
%         Fig 5.9
%  MAKE sure plotSubFigures is =1

figure(7), caFlux2Buffer

orient portrait, print -depsc FIGS_CAFLUX/caflux-2buffermodel.eps
orient portrait, print -dpsc FIGS_CAFLUX/caflux-2buffermodel.ps  

figure(8), caFluxHeuristic

orient portrait, print -depsc FIGS_CAFLUX/caflux-heuristicmodel.eps
orient portrait, print -dpsc FIGS_CAFLUX/caflux-heuristicmodel.ps  

