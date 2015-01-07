%% Test Putting with underactuated RR arm

clc; clear; close all;

%% Define constants and parameters

% Simulation Values 
% system is continous
SIM.discrete = false;
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

%% Calculate desired final velocity of the arm endeffector

% distance to the hole
dist = 0.95;
% hole radius
hole_rad = 0.05;
% total length of ball trajectory
l = dist + hole_rad;
% kinetic friction constant
mu_k = 0.6;
% mass of the ball
m_b = 1;
% gravity is already defined
% initial desired velocity of the ball
vin = sqrt(2*mu_k*g*l);
% calculate vf using conservation of momentum
vf = m_b * vin / (m1 + m2);

%% Generate inputs for a desired trajectory

% TODO: implement Jacobian and inverse computations of
% q, qd, qdd from x, xd, xdd
% Put Jacobian in Kinematics

h = SIM.h;
tin = 0; tfin = 1;
t = tin:h:tfin;
amp = pi/6;
sigma = 0.5;
mean = 0.5;
% trajectory consists of a compact gaussian
theta = amp * exp(-0.5*((t-mean).^2)/(sigma^2));
theta = theta - theta(1);
rad = l1 + l2;
y_des = -rad * sin(theta);
x_des = rad * cos(theta);
ref = [x_des; y_des]; % displacement profile 
% modify ref slightly to avoid singularity in inverse kinematics 
ref = ref - 1e-6;
% modify time profile to meet desired final velocity 
% in y_des due to 90 deg rotation
scale = vf/((y_des(end)-y_des(end-1))/h); 
t = t / scale;

% initialize model
SIM.h = SIM.h / scale;
rr = RR(PAR,CON,COST,SIM);
% create trajectory in cartesian space
traj = rr.generateInputs(t,ref,0);

%% Evolve system dynamics and animate the robot arm

%q0 = traj.s(:,1);
q0 = rr.invKinematics(traj.s(:,1));
% add nonzero velocity
%q0(3:4) = q0(3:4) + 0.1*rand(2,1);
% observe output
y = rr.evolve(t,q0,traj.unom);
% add performance to trajectory
traj.addPerformance(traj.unom,y,rr.COST,'Computed Torque - IDM');

% Plot the controls and animate the robot arm
rr.plot_inputs(traj);
rr.plot_outputs(traj);
% modify the animation
rr.animateArm(y(1:2,:),ref);