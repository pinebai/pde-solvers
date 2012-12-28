load OUT/out_wagner_ode.dat
t=out_wagner_ode(:,1);
c=out_wagner_ode(:,2);
h=out_wagner_ode(:,3);
plot(t,c,t,h), legend('c', 'h')
title('t vs. c & h, for out_wagner_ode.dat')
xlabel('time')
ylabel('concentration \mu M')

