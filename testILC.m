%% Test script for ILC (Arimoto-style and Model-based ILC tests)

%# store breakpoints
tmp = dbstatus;
save('tmp.mat','tmp')

%# clear all
close all
clear classes %# clears even more than clear all
clc

%# reload breakpoints
load('tmp.mat')
dbstop(tmp)

%# clean up
clear tmp
delete('tmp.mat')

%% ILC Example 
% taken from the ILC survey by Bristow et al.

dimx = 2;
dimu = 1;

% simulation variables
t0 = 0;
tf = 50;
h = 1;
t = t0:h:tf;
N = length(t)-1;
numIt = 100;

% system and weighting matrices
Q = 1;
R = eye(dimu);
% discrete time matrices
A = [0 1; 1.8 -0.81];
B = [0;1];
C = [0 1];

% Create the 2-norm cost function
cost.fnc = @(y,s) diag((y-s)'*Q*(y-s));

% track smoothed step
s = 1./(1+exp(-(2*t-20)/10));

% initialize states and inputs
x0 = [0;s(1)];
x = zeros(dimx,N+1);
x(:,1) = x0;
y(:,1) = C*x0;
u0 = zeros(dimu,1);
u = zeros(dimu,N);
u(:,1) = u0;

% evolve model with zero input
for j = 1:N
    x(:,j+1) = A*x(:,j) + B*u(:,j);
    y(:,j+1) = C*x(:,j+1);
end

% form a trajectory
trj = Trajectory(t,s,u,[]);
trj.addPerformance(u,y,cost,'zeros');

% create the simpler ilc
ilc = bILC(trj);

% Perform ILC updates
for i = 1:numIt
    % Update the controls
    u = ilc.feedforward(trj,y);
    % Evolve system with both ILC inputs
    for j = 1:N
        x(:,j+1) = A*x(:,j) + B*u(:,j);
        y(:,j+1) = C*x(:,j+1);
    end
    % Add performances
    trj.addPerformance(u,y,cost,ilc);

end

figure(1);
plot(1:numIt,ilc.error);
title('Squared-2-Norm of ILC error');
figure(2);
plot(t,y,'.-',t,s,'-');
title('Last iteration result');
legend('ILC trajectory','Reference');

%% More complicated example
% Example taken from http://www.egr.msu.edu/classes/me851/jchoi/lecture/Lect_14.pdf

close all; clear; clc;
dimx = 3;
dimu = 1;
dimy = 2;

% simulation variables
t0 = 0;
tf = 1;
h = 0.02;
t = t0:h:tf;
N = length(t)-1;

% system and weighting matrices
Q = 100*eye(dimy);
R = 0.01*eye(dimu);
% continuous time matrices
A = [0 1 0; 0 0 1; -0.4 -4.2 -2.1];
B = [0;0;1];
C = [1 0 0; 0 1 0];

% create the structures
SIM.discrete = false;
SIM.dimx = dimx;
SIM.dimu = dimu;
SIM.dimy = dimy;
SIM.h = h;
SIM.eps = 0;
SIM.int = 'Euler';
PAR.A = A;
PAR.B = B;
PAR.C = C;
CON = [];
COST.Q = Q;
COST.R = R;
lin = Linear(PAR,CON,COST,SIM);

% track the sin trajectory
%s = sin(pi/6*t);
s = t.^2;
%s = [s; 2*t];
s = [s; 2*t];

% create yin with zero velocity
x0 = zeros(dimx,1);
y0 = C * x0;

% create trajectory and execute LQR
traj = lin.trajectory(t,y0,s);
[y,us] = lin.observeWithFeedback(traj,x0);
traj.addPerformance(us,y,lin.COST,'LQR');

% or instead create zero input and observe outcome
% us = zeros(dimu,N); 
% y = lin.observe(t,x0,us);
% traj = Trajectory(t,s,us,[]);
% traj.addPerformance(us,y,lin.COST,'zeros');

lin.plot_inputs(traj);
lin.plot_outputs(traj);

% Create an ilc controller
%ilc = bILC(traj);
ilc = mILC(lin,traj);
ilc.u_last = us;
num_trials = 1;

for i = 1:num_trials
    
    us = ilc.feedforward(traj,y);     
    % observe the actual states
    xact = lin.evolve(t,x0,us);
    % get the measurements
    y = lin.observe(t,x0,us);
    traj.addPerformance(us,y,lin.COST,ilc);

end

lin.plot_inputs(traj);
lin.plot_outputs(traj);

figure;
plot(1:num_trials,ilc.error);
title('Squared-2-Norm of ILC error');
legend(ilc.name);