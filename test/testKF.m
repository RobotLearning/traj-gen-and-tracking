%% Test Kalman Filter

clc; clear; close all;

%% Simple Test for Kalman Filter

seed = 5;
rng(seed);
A = 0.5 * rand(2);
C = [0,1];
eps = 0.1;
N = 20;
x0 = [10;1];
x(:,1) = x0;
y(1) = C*x0 + sqrt(eps) * randn;
for i = 2:N
    x(:,i) = A*x(:,i-1);
    y(i) = C*x(:,i) + sqrt(eps) * randn;
end

% initialize KF
mats.A = A;
mats.B = 0;
mats.C = C;
mats.D = 0;
% process noise var
mats.O = 0;
% observation noise var
mats.M = eps;
filter = KF(2,mats);
filter.initState(x0,eps);
for i = 1:N-1
    filter.update(y(i),0);
    yKF(i) = C * filter.x;
    filter.predict(0);
end
filter.update(y(i),0);
yKF(N) = C * filter.x;

w_cut = 45/180;
yButter = filterButter2nd(y,w_cut);

plot(1:N, x(2,:), 'ks-', 1:N, y, 'b-', 1:N, yKF, 'rx:', 1:N, yButter, 'g*');
%legend('traj','noisy traj','Kalman Filter', 'Butterworth');
SSE(1) = (yKF - C*x)*(yKF - C*x)';
SSE(2) = (yButter - C*x)*(yButter - C*x)';

%% Set system parameters and constraints

clc; clear; close all;
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
%lin.plot_outputs(traj);

%% See performance of Kalman Filter

close all;
eps = 1e-1;
mats.A = lin.Ad;
mats.B = lin.Bd;
mats.C = lin.C;
mats.D = zeros(dimy,dimu);
% process noise covar
mats.O = 0 * eye(dimx);
% observation noise covar
mats.M = eps * eye(dimy);
filter = KF(dimx,mats);

lin.SIM.eps_m = eps;
yNoisy = lin.observe(t,x0,us);

filter.initState(x0,eps*eye(dimx));
for i = 1:Nu
    filter.update(yNoisy(i),us(i));
    yKF(i) = lin.C * filter.x;
    filter.predict(us(i));
end
filter.update(yNoisy(N),0);
yKF(N) = lin.C * filter.x;

plot(t,y,'ks-',t,yNoisy,'b-',t,yKF,'rx:');
%legend('traj','noisy traj','estimated traj');
SSE(1) = (yNoisy - y)*(yNoisy - y)';
SSE(2) = (yKF - y)*(yKF - y)';