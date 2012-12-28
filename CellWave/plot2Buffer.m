function plot2Buffer( tcomp, q )

tmax=max(tcomp);
tmin=min(tcomp);

qmax=max( [max( q.c ), max( q.h ), max( q.b1), max( q.b2)] );

plot( tcomp, q.c, tcomp, q.h, tcomp, q.b1, tcomp, q.b2 );
axis([tmin tmax -0.1 qmax+0.1 ]),  grid
xlabel('time'), ylabel('concentration')
legend( 'c', 'h', 'b1', 'b2')
