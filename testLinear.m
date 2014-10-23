%% Simulate trajectories for the 2D linear dynamics model

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

%% Set system parameters and constraints

% parameter values of the experimental setup

% Simulation Values 
% for simulation purposes notify if model already discretized
SIM.discrete = true;
% dimension of the x vector
SIM.dimx = 2;
% dimension of the output vector
SIM.dimy = 1;
% dimension of the control input
SIM.dimu = 1;
% time step h 
SIM.h = 1;
% measurement noise 
SIM.eps = 0;
% integration method if continous model is given
SIM.int = 'Euler';

% discrete time matrices
PAR.Ad = [0 1; 1.8 -0.81];
PAR.Bd = [0;1];
PAR.C = [0 1];

% constraints on the system dynamics
% TODO
CON = [];

% cost structure
% only penalize positions
COST.Q = 1;
COST.R = 1;

% initialize model
lin = Linear(PAR,CON,COST,SIM);

%% Create desired trajectory 

dimu = SIM.dimu;
dimx = SIM.dimx;
h = SIM.h;
tfin = 50;
t = 0:h:tfin;
N = length(t);
Nu = N - 1;

% track smoothed step
s = 1./(1+exp(-(2*t-20)/10));

%% Evolve system dynamics with a nominal input

% initialize states and inputs
x0 = [0;s(1)];
us = zeros(dimu,Nu);
% observe output
y = lin.observe(t,x0,us);
traj = Trajectory(t,s,us,[]);
traj.addPerformance(us,y,lin.COST,'zeros');

%% Iterative Learning Control

num_trials = 1;
ilc = mILC(lin,traj);

for i = 1:num_trials
    
    us = ilc.feedforward(traj,y);     
    % observe the actual states
    xact = lin.evolve(t,x0,us);
    % get the measurements
    y = lin.observe(t,x0,us);
    traj.addPerformance(us,y,lin.COST,ilc);
    
end

% Plot the controls
lin.plot_inputs(traj);
lin.plot_outputs(traj);

figure;
plot(1:num_trials,ilc.error);
title('Squared-2-Norm of ILC error');
legend('Monotonic ILC');

%% More complicated example

clc; close all; clear classes; 
dimx = 3;
dimu = 1;
dimy = 1;

% simulation variables
t0 = 0;
tf = 50;
h = 1;
t = t0:h:tf;
N = length(t)-1;

% system and weighting matrices
Q = 1;
R = 1;
% continuous time matrices
A = [0 1 0; 0 0 1; -0.4 -4.2 -2.1];
B = [0;0;1];
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

% track sin for the measured variable y = x1
s = sin(2*pi*t);
x0 = zeros(dimx,1);
y0 = C * x0;

% create trajectory and execute LQR
lin.trajectory(t,y0,s);

[y,us] = lin.observeWithFeedback(trj,x0,K);
trj.addPerformance(us,y,lin.COST,'LQR');

lin.plot_inputs(trj);
lin.plot_outputs(trj);

% Create an ilc controller
ilc = bILC(trj);
num_trials = 5;

for i = 1:num_trials
    
    us = ilc.feedforward(traj,y);     
    % observe the actual states
    xact = lin.evolve(t,x0,us);
    % get the measurements
    y = lin.observe(t,x0,us);
    traj.addPerformance(us,y,lin.COST,ilc);

    lin.plot_inputs(trj);
    lin.plot_outputs(trj);
end