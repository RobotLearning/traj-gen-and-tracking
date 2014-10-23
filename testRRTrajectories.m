%% Simulate trajectories for the planar RR arm
%
% TODO: Investigate the following
%
% Other models (RRR arm, ...)
% Other mismatches (parametric, friction, ...)
% Adding feedback 
% DMP 

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

% form constraints
CON.link1.q.max = Inf;
CON.link1.q.min = -Inf;
CON.link1.qd.max = Inf;
CON.link1.qd.min = -Inf;
CON.link1.qdd.max = Inf;
CON.link1.qdd.min = -Inf;
CON.link1.u.max = Inf;
CON.link1.u.min = -Inf;
CON.link1.udot.max = 100;
CON.link1.udot.min = -100;
CON.link2.q.max = Inf;
CON.link2.q.min = -Inf;
CON.link2.qd.max = Inf;
CON.link2.qd.min = -Inf;
CON.link2.qdd.max = Inf;
CON.link2.qdd.min = -Inf;
CON.link2.u.max = Inf;
CON.link2.u.min = -Inf;
CON.link2.udot.max = 100;
CON.link2.udot.min = -100;

% cost structure
% only penalize positions
COST.Q = diag([1,1,0,0]);

% simulation parameters
SIM.discrete = false;
SIM.h = 0.02;
SIM.eps = 3e-10;
SIM.eps_d = 3e-10;
SIM.int = 'Euler'; % or RK4

% initialize model
rr = RR(PAR,CON,COST,SIM);

%% Generate a desired trajectory

% TODO: implement Jacobian and inverse computations of
% q, qd, qdd from x, xd, xdd
% Put Jacobian in Kinematics

h = SIM.h;
y_des = 0.4:h:0.6;
yd_des = [0, diff(y_des)];
x_des = 0.6 * ones(1,length(y_des));
xd_des = [0, diff(x_des)];
t = h * 1:length(y_des);
s = [x_des; y_des; xd_des; yd_des]; % desired trajectory 
Traj = rr.trajectory(t,s);

%% Evolve system dynamics and animate the robot arm

% TODO: add a nonzero friction matrix B

q0 = [rr.q(:,1); rr.qd(:,1)];
% observe output
q_obs = rr.evolve(t,q0,Traj.unom);

% get the deviations
% TODO: xd should also be returned
[~,y] = rr.kinematics(q_obs(1:2,:));
yd = [zeros(2,1), diff(y')'];
% add performance to trajectory
Traj.addPerformance(Traj.unom,[y;yd],rr.COST,'Nominal model inversion');

% Plot the controls and animate the robot arm
rr.plot_inputs(Traj);
rr.plot_outputs(Traj);
rr.animateArm(q_obs(1:2,:),s(1:2,:));

%% Start learning with ILC

num_trials = 10;

% get deviation
dev = Traj.PERF(end).err;
ilc = aILC(rr,Traj);
%ilc = bILC(Traj);

for i = 1:num_trials
    % get next inputs
    u = ilc.feedforward(Traj,rr,dev);
    % evolve system
    q_obs = rr.evolve(t,q0,u);
    % get the cartesian coordinates
    [~,y] = rr.kinematics(q_obs(1:2,:));
    yd = [zeros(2,1), diff(y')'];
    % add performance to trajectory
    Traj.addPerformance(u,[y;yd],rr.COST,ilc);
    dev = Traj.PERF(end).err;
    % Plot the controls and animate the robot arm
    %RR.plot_controls(Traj);
    %RR.animateArm(q_act(1:2,:),s(1:2,:));
end

% Plot the controls and animate the robot arm
rr.plot_inputs(Traj);
rr.plot_outputs(Traj);
rr.animateArm(q_obs(1:2,:),s(1:2,:));