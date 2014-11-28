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
ref = 1./(1+exp(-(2*t-20)/10));

%% Evolve system dynamics with a nominal input

% initialize states and inputs
x0 = [0;ref(1)];
us = zeros(dimu,Nu);
% observe output
y = lin.observe(t,x0,us);
traj = Trajectory(t,ref,us,[]);
traj.addPerformance(us,y,lin.COST,'zeros');

%% Iterative Learning Control

num_trials = 1;
%ilc = bILC(traj);
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
% Example taken from http://www.egr.msu.edu/classes/me851/jchoi/lecture/Lect_14.pdf

close all;
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
Q = diag([100,0,0]);
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