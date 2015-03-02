%% Simulate trajectories for the pendulum
%

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

%% Define constants and parameters

% Simulation Values 
% system is continous
SIM.discrete = false;
% learn in cartesian space
SIM.cartesian = false;
% dimension of the x vector
SIM.dimx = 2;
% dimension of the output y
SIM.dimy = 2;
% dimension of the control input
SIM.dimu = 1;
% time step h 
SIM.h = 0.01;
% noise and initial error
SIM.eps = 3e-10;
SIM.eps_d = 3e-10;
% integration method
SIM.int = 'Euler';
% trajectory in joint space?
SIM.jref = false;

% constants
g = 9.81;
% joint parameters
m = 1; %mass of first link, kg
l = 0.50; %length of first link, m
l_c = 0.25; %distance of link's center of gravity to base, m
I = (1/12)*m*l^2; %assume thin rod moment of inertia around c.o.g.
% motor parameters
J_a = 0.100; % actuator inertia of link
J_g = 0.050; % gear inertia of link1
J_m = J_a + J_g;
% gear ratio typically ranges from 20 to 200 or more - comment from book
r = 20;
 
% pass it as a parameter structure
PAR.const.g = g;
PAR.link.mass = m;
PAR.link.length = l;
PAR.link.centre.dist = l_c;
PAR.link.inertia = I;
PAR.link.motor.inertia = J_m;
PAR.link.motor.gear_ratio = r;

% TODO: we should observe joint angles only
PAR.C = eye(SIM.dimy,SIM.dimx);

% form constraints
CON = [];

% cost structure
% only penalize positions
COST.Q = 100*diag([1,0]);
COST.R = 1 * eye(SIM.dimu);

% initialize model
r = R(PAR,CON,COST,SIM);

%% Generate inputs for a desired trajectory

% TODO: implement Jacobian and inverse computations of
% q, qd, qdd from x, xd, xdd
% Put Jacobian in Kinematics

h = SIM.h;
tin = 0; tfin = 1;
t = tin:h:tfin;
% trajectory consists of a compact gaussian
amp = pi/6;
sigma = 0.5;
mean = 0.5;
theta = amp * exp(-0.5*((t-mean).^2)/(sigma^2));
theta = theta - theta(1);
rad = l;
x_des = rad * cos(theta);
y_des = rad * sin(theta);
ref = [x_des; y_des]; % displacement profile 
traj = r.generateInputs(t,ref); % trajectory generated in cartesian space

%% Evolve system dynamics and animate the robot arm

% TODO: add a nonzero friction matrix B

q0 = traj.s(:,1);
% q0 = r.invKinematics(traj.s(:,1));
% q1 = r.invKinematics(traj.s(:,2));
% q0(2,1) = (q1 - q0)/h;
% add nonzero velocity
% q0(2) = q0(2) + 0.1*randn;
% observe output
qact = r.evolve(t,q0,traj.unom);
% get the cartesian coordinates
% y = r.kinematics(qact(1,:));
% % add cartesian velocities
% yd = diff(y')'/ h; yd(:,end+1) = yd(:,end);
% y = [y;yd];
% add performance to trajectory
traj.addPerformance(traj.unom,qact,r.COST,'Inverse Dynamics');

% Plot the controls and animate the robot arm
r.plot_inputs(traj);
r.plot_outputs(traj);
r.animateArm(qact(1,:),ref);

%% Start learning with ILC

num_trials = 10;

%ilc = aILC(r,traj);
ilc = mILC(r,traj); 

for i = 1:num_trials
    % get next inputs
    u = ilc.feedforward(traj,qact);
    % evolve complete system
    qact = r.evolve(t,q0,u);
    % get the cartesian coordinates
    %y = r.kinematics(qact(1,:));
    % add cartesian velocities
    %yd = diff(y')'/ h; yd(:,end+1) = yd(:,end);
    %y = [y;yd];
    % add performance to trajectory
    traj.addPerformance(u,qact,r.COST,ilc);
    % Plot the controls and animate the robot arm
    %r.animateArm(qact(1,:),ref);
end

% Plot the controls and animate the robot arm
r.plot_inputs(traj);
r.plot_outputs(traj);
r.animateArm(qact(1,:),ref);