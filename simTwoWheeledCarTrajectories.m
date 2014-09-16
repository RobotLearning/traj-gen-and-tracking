%% Simulate trajectories for the two-wheeled robot kinematical model

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
% radius of the robot
PAR.wheel1.radius = 0.2; % meter
PAR.wheel2.radius = 0.2;
PAR.length = 0.8;

% constraints on the system dynamics
CON.state.x.max = 5; 
CON.state.x.min = -5;
CON.state.y.max = 5;
CON.state.y.min = -5;
CON.state.theta.max = 2*pi;
CON.state.theta.min = -2*pi;
CON.wheel1.u.max = 200;
CON.wheel1.u.min = -200;
CON.wheel2.u.max = 200;
CON.wheel2.u.min = -200;
CON.wheel1.udot.max = 100;
CON.wheel1.udot.min = -100;
CON.wheel2.udot.max = 100;
CON.wheel2.udot.min = -100;

% Simulation Values 
% dimension of the x vector
SIM.dimx = 3;
% dimension of the control input
SIM.dimu = 2;
% time step h 
SIM.h = 0.1;
% noise and initial disturbance covariance
SIM.eps = 0.003;
SIM.eps_d = 0.005;
% integration method
SIM.int = 'Euler';

% cost structure
% only penalize positions
COST.Q = diag([1,1,0]);

% initialize model
TW = TwoWheeledCar(PAR,CON,COST,SIM);

%% Create desired trajectory 

h = SIM.h;
tfin = 1;
t = h:h:tfin;
% create a pre-nominal sine-curve 
s(1,:) = t;
s(2,:) = sin(2*pi*t);
s(3,1:end-1) = atan2(diff(s(2,:)),diff(s(1,:)));
s(3,end) = s(3,1); % since trajectory is periodical after t = 1

Traj = TW.trajectory(t,s);

%% Evolve system dynamics and animate the robot arm

x0 = s(:,1);
xact = TW.evolve(t,x0,Traj.unom);

% Plot the controls and animate the robot arm
TW.plot_controls(Traj);
TW.animate(xact,s(1:2,:));

%% Start learning with ILC

num_trials = 2;

% get the deviations
% add performance to trajectory
Traj.addPerformance(Traj.unom,xact,TW.COST,'Nominal');
%STR = construct_model(t,s,PAR,Traj.unom,CON);
ilc = aILC(TW,Traj);

for i = 1:num_trials
    % get next inputs
    u = ilc.feedforward(Traj,TW,Traj.PERF(end).err);
    % evolve system
    x_act = TW.evolve(t,x0,u);
    % add performance to trajectory
    Traj.addPerformance(u,x_act,TW.COST,ilc);
    % Plot the controls and animate the robot arm
    %RR.plot_controls(Traj);
    %RR.animateArm(q_act(1:2,:),s(1:2,:));
end

% Plot the controls and animate the robot arm
TW.plot_controls(Traj);
TW.animate(x_act,s(1:2,:));
