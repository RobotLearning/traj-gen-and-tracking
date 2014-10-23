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
CON.wheel1.u.max = 200;
CON.wheel1.u.min = -200;
CON.wheel2.u.max = 200;
CON.wheel2.u.min = -200;
CON.wheel1.udot.max = 100;
CON.wheel1.udot.min = -100;
CON.wheel2.udot.max = 100;
CON.wheel2.udot.min = -100;

% Simulation Values 
% system is continous
SIM.discrete = false;
% dimension of the x vector
SIM.dimx = 3;
% dimension of the control input
SIM.dimu = 2;
% time step h 
SIM.h = 0.01;
% noise and initial error
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

dim_u = SIM.dimu;
dim_x = SIM.dimx;
h = SIM.h;
tfin = 1;
t = h:h:tfin;
N = length(t);
Nu = N - 1;
% create a pre-nominal sine-curve 
s(1,:) = t;
s(2,:) = sin(2*pi*t);
s(3,1:end-1) = atan2(diff(s(2,:)),diff(s(1,:)));
s(3,end) = s(3,1); % since trajectory is periodical after t = 1

Traj = TW.trajectory(t,s);

%% Evolve system dynamics and animate the robot arm

x0 = s(:,1);
xact = TW.evolve(t,x0,Traj.unom);

% add performance to trajectory
Traj.addPerformance(Traj.unom,xact,TW.COST,'Nominal model inversion');

% Plot the controls and animate the robot arm
TW.plot_inputs(Traj);
TW.plot_outputs(Traj);
TW.animate(xact,s(1:2,:));

%% Iterative Learning Control

num_trials = 10;
ilc = aILC(TW,Traj);
% TODO: add C matrix and observe output with observe function
y = TW.evolve(t,x0,Traj.unom);
dev = y - s;

for i = 1:num_trials
    
    u = ilc.feedforward(Traj,TW,dev);    
    % get error (observed trajectory deviation)
    y = TW.evolve(t,x0,u);
    %y = TW.observe(t,x0,u);
    dev = y - s;
    Traj.addPerformance(u,y,TW.COST,ilc);
    
end

% Plot the controls and animate the robot arm
TW.plot_inputs(Traj);
TW.plot_outputs(Traj);
TW.animate(y,s(1:2,:));