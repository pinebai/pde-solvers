load OUT/out_slep2buf.dat
t=out_slep2buf(:,1);
c=out_slep2buf(:,2);
h=out_slep2buf(:,3);
b1=out_slep2buf(:,4);
b2=out_slep2buf(:,5);

plot(t,c,'bo-',t,h,'k',t,b1,'rx',t,b2,'g+')
legend('c', 'h', 'b1','b2')
title('t vs. c & h, for out_slep2buf.dat')
xlabel('time')
ylabel('concentration \mu M')

