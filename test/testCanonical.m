% test script for canonical system

clc; close all; clear classes; 

%% Test discrete canonical system
h = 0.01; ax = 2; tau = 1; nst = 100; pat = 'd';
can = Canonical(h,ax,tau,nst,pat);
x_tr = can.evolve();
plot(x_tr);
legend('Canonical system - discrete');
title('Phase rollout');
xlabel('time (number of steps)');
ylabel('phase x');

%% Test rhythmic canonical system
h = 0.01; ax = 2; tau = 1; nst = 100; pat = 'r';
can = Canonical(h,ax,tau,nst,pat);
x_tr = can.evolve();
figure;
plot(x_tr);
legend('Canonical system - rhythmic');
title('Phase rollout');
xlabel('time (number of steps)');
ylabel('phase x');