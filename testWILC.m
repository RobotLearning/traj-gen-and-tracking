%% Example to test DMP-learning wILC algorithm on linear systems with drift
% Example taken from http://www.egr.msu.edu/classes/me851/jchoi/lecture/Lect_14.pdf

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

%% 
dimx = 3;
dimu = 1;
dimy = 3;

% simulation variables
t0 = 0;
tf = 1;
h = 0.02;
t = t0:h:tf;
N = length(t)-1;

% system and weighting matrices
Q = 100*eye(dimy);
%Q = diag([100,0,0]);
R = 1*eye(dimu);
% continuous time matrices
A = [0 1 0; 0 0 1; -0.4 -4.2 -2.1];
B = [0;0;1];
% partial observation model
C = eye(dimy);
%C = [1 0 0; 0 1 0];
%C = [1 0 0];

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
ref = t.^2;
ref = [ref; 2*t];
ref = [ref; 0, 2*ones(1,length(t)-1)];

% create yin with zero velocity
x0 = zeros(dimx,1);
y0 = C * x0;
% y0 = y0(1:2);

% create trajectory and execute LQR
traj = lin.trajectory(t,y0,ref);
[y,us] = lin.observeWithFeedback(traj,x0);
traj.addPerformance(us,y,lin.COST,'LQR');
% us = zeros(dimu,N); 
% y = lin.observe(t,x0,us);
% traj = Trajectory(t,ref,us,[]);
% traj.addPerformance(us,y,lin.COST,'zeros');

lin.plot_inputs(traj);
lin.plot_outputs(traj);

% Create an ilc controller
%ilc = bILC(traj);
ilc = mILC(lin,traj);
ilc.u_last = us;
num_trials = 1;

% estimate d
a = [0.1;0.1;0.1];
dhat = zeros(dimx,N+1);
dhat(:,1) = 0;
for i = 1:N
    dhat(:,i+1) = lin.Ad * dhat(:,i) + a;
end

for i = 1:num_trials
    
    us = ilc.feedforward(traj,y);     
    % observe the actual states
    xact = lin.evolve(t,x0,us);
    xnom = lin.simulate(t,x0,us,@lin.nominal);
    d = xact - xnom;
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