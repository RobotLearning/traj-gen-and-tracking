%% ILC on RR modifying DMPs now
%
clc; clear; close all;

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
% measurement noise covariance
SIM.eps_m = 3e-10;
% integration method
SIM.int = 'Euler';
% trajectory in joint space?
SIM.jref = false;

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
I2 = (1/12)*m2*l2^2; %kg m^2
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

% we should observe joint angles only
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

h = SIM.h;
tin = 0; tfin = 1;
t = tin:h:tfin;
y_des = 0.4 + 0.4 * t;
x_des = 0.6 - 0.2 * t;
ref = [x_des; y_des]; % displacement profile 

%traj = rr.generateInputs(t,ref);
bfs = 100;
% trajectory generated in joint space
[traj,dmp] = rr.generateInputsWithDMP(t,bfs,ref); 

%save dmp weights for later use
w_origin = zeros(length(dmp),bfs);
for i = 1:length(dmp)
    w_origin(i,:) = dmp(i).w;
end

% Generate feedback with LQR
rr.generateFeedback(traj);

%% Evolve system dynamics and animate the robot arm

q0 = traj.s(:,1);
% add zero velocity or perturb initial velocity
%q0(3:4) = 0;
% observe output
%qact = rr.observeWithFeedbackErrorForm(traj,q0);
qact = rr.observeWithFeedbackErrorForm(traj,q0,dmp);
% add performance to trajectory
traj.addPerformance([],qact,rr.COST,'ID + LQR');

% Plot the controls and animate the robot arm
rr.plot_outputs(traj);
%rr.animateArm(qact(1:2,:),ref);

%% Start learning with ILC

num_trials = 5;
%ilc = wILC(traj,rr,'t');
ilc = wILC(traj,rr,'dmp');

for i = 1:num_trials
    % get next inputs
    dmp = ilc.feedforward(traj,dmp,qact);
    %traj2 = ilc.feedforward(traj,[],qact);   
    % get the measurements
    qact = rr.observeWithFeedbackErrorForm(traj,q0,dmp);
    %qact = rr.observeWithFeedbackErrorForm(traj2,q0);
    traj.addPerformance([],qact,rr.COST,ilc);
    % Plot the controls and animate the robot arm
    %rr.animateArm(qact(1:2,:),ref);
end

% Plot the controls and animate the robot arm
rr.plot_outputs(traj);
rr.animateArm(qact(1:2,:),ref);

%% Change the hitting point slightly and see how the robot is doing

%{
% new reference in joint space
%rr.flag_ref_jsp = true;
% center of the region
g0 = traj.s(:,end);
% radius of region
r = 0.05;
% sample a point from that region
g = g0 + [r;r;0;0];
[dmpNew,jNew] = adaptDMP(q0,g,dmp,w_origin,length(t));
[~,refNew] = rr.kinematics(jNew);

trajNew = rr.generateInputs(t,refNew);
rr.generateFeedback(trajNew);

q0 = trajNew.s(:,1);
% observe output
%qact = rr.observeWithFeedbackErrorForm(trajNew,q0);
qact = rr.observeWithFeedbackErrorForm(trajNew,q0,dmpNew);
% add performance to trajectory
trajNew.addPerformance([],qact,rr.COST,'ID + LQR');

%ilc1 = wILC(trajNew,rr,'t');
ilc1 = wILC(trajNew,rr,'dmp');

for i = 1:num_trials
    % get next inputs
    dmpNew = ilc1.feedforward(trajNew,dmpNew,qact);
    %traj3 = ilc1.feedforward(trajNew,[],qact);   
    % get the measurements
    qact = rr.observeWithFeedbackErrorForm(trajNew,q0,dmpNew);
    %qact = rr.observeWithFeedbackErrorForm(traj3,q0);
    trajNew.addPerformance([],qact,rr.COST,ilc1);
    % Plot the controls and animate the robot arm
    %rr.animateArm(qact(1:2,:),ref);
end

% Plot the controls and animate the robot arm
rr.plot_outputs(trajNew);
rr.animateArm(qact(1:2,:),refNew);

% Lets see how well the learned dmp generalizes
%qact = rr.observeWithFeedbackErrorForm(traj2,q0);
qact = rr.observeWithFeedbackErrorForm(trajNew,q0,dmp);
trajNew.addPerformance([],qact,rr.COST,'OLD ILC');

% Plot the controls and animate the robot arm
rr.plot_outputs(trajNew);
rr.animateArm(qact(1:2,:),refNew);

%ilc2 = wILC(trajNew,rr,'t');
ilc2 = wILC(trajNew,rr,'dmp');

for i = 1:num_trials
    % get next inputs
    dmp = ilc2.feedforward(trajNew,dmp,qact);
    %traj3 = ilc2.feedforward(trajNew,[],qact);   
    % get the measurements
    qact = rr.observeWithFeedbackErrorForm(trajNew,q0,dmp);
    %qact = rr.observeWithFeedbackErrorForm(traj3,q0);
    trajNew.addPerformance([],qact,rr.COST,ilc2);
    % Plot the controls and animate the robot arm
    %rr.animateArm(qact(1:2,:),ref);
end

% Plot the controls and animate the robot arm
rr.plot_outputs(trajNew);
rr.animateArm(qact(1:2,:),refNew);

%}