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

% Uncoupled springs
PAR.m1 = 0.5;
PAR.k1 = 10;
PAR.b1 = 0.1;
PAR.m2 = 1.5;
PAR.k2 = 10;
PAR.b2 = 0.2;
        
PAR.state.init = zeros(4,1);

% constraints on the system dynamics
CON.state.x.max = 5; 
CON.state.x.min = -5;
CON.state.y.max = 5;
CON.state.y.min = -5;
CON.u1.max = 200;
CON.u1.min = -200;
CON.u2.max = 200;
CON.u2.min = -200;
CON.u1.dot.max = 100;
CON.u1.dot.min = -100;
CON.u2.dot.max = 100;
CON.u2.dot.min = -100;

% Simulation Values 
% dimension of the x vector
SIM.dimx = 4;
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
COST.Q = diag([1,0,1,0]);

% initialize model
lin = Linear2DDynamics(PAR,CON,COST,SIM);

%% Create desired trajectory 

dim_u = SIM.dimu;
dim_x = SIM.dimx;
h = SIM.h;
tfin = 1;
t = h:h:tfin;
N = length(t);
Nu = N - 1;
% make the linear system go around in a circle
% later increase the tfin and use rhythmic DMP
s(1,:) = cos(2*pi*t);
s(2,:) = sin(2*pi*t);

Traj = lin.trajectory(t,s);

%% Evolve system dynamics

x0 = PAR.state.init;
xact = lin.evolve(t,x0,Traj.unom);

% Plot the controls and animate the robot arm
lin.plot_controls(Traj);
lin.plot_states(xact([1,3],:),s);

% add performance to trajectory
Traj.addPerformance(Traj.unom,xact,lin.COST,'Nominal');

%% Iterative Learning Control

num_trials = 10;
ilc = bILC(lin,Traj);
% observe output
y = lin.observe(t,x0,Traj.unom);
dev = y - s;

for i = 1:num_trials
    
    u = ilc.feedforward(Traj,lin,dev);    
    % get error (observed trajectory deviation)
    xact = lin.evolve(t,x0,u);
    y = lin.observe(t,x0,u);
    dev = y - s;
    Traj.addPerformance(u,xact,lin.COST,ilc);
    
end

% Plot the controls and animate the robot arm
lin.plot_controls(Traj);
lin.plot_states(xact([1,3]),s);