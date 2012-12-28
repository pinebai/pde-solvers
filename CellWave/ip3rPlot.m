function ip3rPlot( c, jflux, jpump )

semilogx( c, jflux, c,jpump,'--')
xlabel('log(Ca2+)');
ylabel('Ca2+ efflux \mu M/sec');
grid
