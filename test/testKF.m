%% Test Kalman Filter and EKF

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

%% Simple Test for Kalman Filter

clc; clear; close all;
A = 0.9;
C = 1;
eps = 0.01;
N = 50;
x0 = 10 * rand;
x(1) = x0;
for i = 1:N-1
    x(i+1) = A*x(i);
    y(i) = C*x(i) + sqrt(eps) * randn;
end
y(N) = C*x(end) + sqrt(eps) * randn;

% initialize KF
mats.A = A;
mats.B = 0;
mats.C = C;
mats.D = 0;
% process noise var
mats.O = 0;
% observation noise var
mats.M = eps;
filter = KF(1,mats);
filter.initState(y(1),eps);
yKF(1) = C * filter.x;
for i = 1:N-1
    filter.predict(0);
    filter.update(y(i),0);
    yKF(i+1) = C * filter.x;
end

w_cut = 45/180;
yButter = filterButter2nd(y,w_cut);

plot(1:N,x,1:N,y,1:N,yKF,1:N,yButter);
legend('traj','noisy traj','Kalman Filter', 'Butterworth');
SSE = yKF - y

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
SIM.eps_m = 0;
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
dimy = SIM.dimy;
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

%% Find the control inputs

num_trials = 2;
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

%% See performance of Kalman Filter

close all;
eps = 1e-4;
mats.A = lin.Ad;
mats.B = lin.Bd;
mats.C = lin.C;
mats.D = zeros(dimy,dimu);
% process noise covar
mats.O = 0 * eye(dimx);
% observation noise covar
mats.M = eps * eye(SIM.dimy);
filter = KF(dimx,mats);

lin.SIM.eps_m = eps;
yNoisy = lin.observe(t,x0,us);

filter.initState(x0,eps*eye(dimx));
yKF(1) = lin.C * filter.x;
for i = 1:Nu
    filter.predict(us(i));
    filter.update(yNoisy(i),us(i));
    yKF(i+1) = lin.C * filter.x;
end

plot(t,y,t,yNoisy,t,yKF);
legend('traj','noisy traj','estimated traj');