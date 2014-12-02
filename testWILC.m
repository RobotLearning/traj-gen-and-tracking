%% Example to test DMP-learning wILC algorithm on linear systems with drift

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

% Example taken from http://www.egr.msu.edu/classes/me851/jchoi/lecture/Lect_14.pdf
dimx = 3;
dimu = 1;
dimy = 1;

% simulation variables
t0 = 0;
tf = 1;
h = 0.02;
t = t0:h:tf;
N = length(t)-1;

% system and weighting matrices
Q = 100*eye(dimy);
%Q = diag([100,0,0]);
R = 10*eye(dimu);
% continuous time matrices
A = [0 1 0; 0 0 1; -0.4 -4.2 -2.1];
B = [0;0;1];
% partial observation model
%C = eye(dimy);
%C = [1 0 0; 0 1 0];
C = [1 0 0];

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
%ref = [ref; 2*t];
%ref = [ref; 0, 2*ones(1,length(t)-1)];
x0 = zeros(dimx,1);
y0 = C * x0;
% create yin with zero velocity
yin = [y0;0];      

% create trajectory and execute LQR
[traj,dmp] = lin.generateDMP(t,yin,ref);
[y,us] = lin.observeWithDMPFeedback(dmp,traj,x0);
traj.addPerformance(us,y,lin.COST,'LQR'); 
% us = zeros(dimu,N); 
% y = lin.observe(t,x0,us);
% traj = Trajectory(t,ref,us,[]);
% traj.addPerformance(us,y,lin.COST,'zeros');

lin.plot_inputs(traj);
lin.plot_outputs(traj);

% Create an ilc controller
ilc = wILC(lin,traj);
num_trials = 50;

for i = 1:num_trials
    % update the weights of the dmp
    ilc.feedforward(dmp,traj,y);     
    % get the measurements
    [y,us] = lin.observeWithDMPFeedback(dmp,traj,x0);
    traj.addPerformance(us,y,lin.COST,ilc);

end

lin.plot_inputs(traj);
lin.plot_outputs(traj);

figure;
plot(1:num_trials,ilc.error);
title('Squared-2-Norm of ILC error');
legend(ilc.name);