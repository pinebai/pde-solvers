load OUT/out_lirinzel.dat
t=out_lirinzel(:,1);
c=out_lirinzel(:,2);
h=out_lirinzel(:,3);
plot(t,c,t,h), legend('c', 'h')
title('t vs. c & h, for out_lirinzel.dat')
xlabel('time')
ylabel('concentration \mu M')

