%% Simulate trajectories for the planar RR arm
%
% TODO: Investigate the following
%
% Adding feedback 
% slow down with DMP 

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
SIM.dimx = 4;
% dimension of the output y
SIM.dimy = 4;
% dimension of the control input
SIM.dimu = 2;
% time step h 
SIM.h = 0.01;
% noise and initial error
SIM.eps = 3e-10;
SIM.eps_d = 3e-10;
% integration method
SIM.int = 'Euler';

% constants
g = 9.81;
% joint parameters
m1 = 1; %mass of first link, kg
m2 = 0.5; %mass of second link, kg
l1 = 0.50; %length of first link, m
l2 = 0.40; %length of second link, m
l_c1 = 0.25; %distance of first link's center of gravity to prev. joint, m
l_c2 = 0.20; %dist. of second link's c.oZ.g. to prev. joint, m
I1 = (1/12)*m1*l1^2; %assume thin rod moment of inertia around c.o.g.
I2 = (1/12)*m1*l2^2; %kg m^2
% motor parameters
J_a1 = 0.100; % actuator inertia of link 1
J_g1 = 0.050; % gear inertia of link1
J_m1 = J_a1 + J_g1;
J_a2 = 0.080; % actuator inertia of link 2
J_g2 = 0.040; % gear inertia of link2
J_m2 = J_a2 + J_g2;
% gear ratio typically ranges from 20 to 200 or more - comment from book
r_1 = 20;
r_2 = 20;
 
% pass it as a parameter structure
PAR.const.g = g;
PAR.link1.mass = m1;
PAR.link2.mass = m2;
PAR.link1.length = l1;
PAR.link2.length = l2;
PAR.link1.centre.dist = l_c1;
PAR.link2.centre.dist = l_c2;
PAR.link1.inertia = I1;
PAR.link2.inertia = I2;
PAR.link1.motor.inertia = J_m1;
PAR.link2.motor.inertia = J_m2;
PAR.link1.motor.gear_ratio = r_1;
PAR.link2.motor.gear_ratio = r_2;

% TODO: we should observe joint angles only
PAR.C = eye(SIM.dimy,SIM.dimx);

% form constraints
CON = [];

% cost structure
% only penalize positions
COST.Q = 100*diag([1,1,0,0]);
COST.R = 1 * eye(SIM.dimu);

% initialize model
rr = RR(PAR,CON,COST,SIM);

%% Generate inputs for a desired trajectory

% TODO: implement Jacobian and inverse computations of
% q, qd, qdd from x, xd, xdd
% Put Jacobian in Kinematics

h = SIM.h;
tin = 0; tfin = 1;
t = tin:h:tfin;
y_des = 0.4 + 0.2 * t;
x_des = 0.6 * ones(1,length(t));
ref = [x_des; y_des]; % displacement profile 
traj = rr.generateInputs(t,ref); % trajectory generated in joint space

%% Evolve system dynamics and animate the robot arm

% TODO: add a nonzero friction matrix B

q0 = traj.s(:,1);
% add nonzero velocity
%q0(3:4) = q0(3:4) + 0.1*rand(2,1);
% observe output
qact = rr.evolve(t,q0,traj.unom);
% add performance to trajectory
traj.addPerformance(traj.unom,qact,rr.COST,'Inverse Dynamics');

% Plot the controls and animate the robot arm
rr.plot_inputs(traj);
rr.plot_outputs(traj);
rr.animateArm(qact(1:2,:),ref);

%% Start learning with ILC

num_trials = 10;

%ilc = aILC(rr,traj);
ilc = mILC(rr,traj); 

for i = 1:num_trials
    % get next inputs
    u = ilc.feedforward(traj,qact);
    %u = ilc.feedforward(traj,rr,dev);
    % evolve system
    qact = rr.evolve(t,q0,u);
    % get the cartesian coordinates
    %[~,y] = rr.kinematics(qact(1:2,:));
    % add performance to trajectory
    traj.addPerformance(u,qact,rr.COST,ilc);
    % Plot the controls and animate the robot arm
    %rr.animateArm(qact(1:2,:),ref);
end

% Plot the controls and animate the robot arm
rr.plot_inputs(traj);
rr.plot_outputs(traj);
rr.animateArm(qact(1:2,:),ref);