%% Simulate trajectories for the planar RR arm
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
tin = h; tfin = 1;
t = tin:h:tfin;
y_des = 0.4 + 0.4 * t;
x_des = 0.6 - 0.2 * t;
ref = [x_des; y_des]; % displacement profile 

% downsample reference
freq = 50;
freq_original = 100;
rate = freq_original/freq;
idx = rate * (1:length(t)/rate);
ref = ref(:,idx);
t = t(idx);

traj = rr.generateInputs(t,ref); % trajectory generated in joint space

%% Evolve system dynamics and animate the robot arm

% add feedback K to traj
rr.generateFeedback(traj); % form feedback to stabilize

%q0 = [0.0325; 1.2901; 0.7065; -0.2238];
q0 = traj.s(:,1);
% add nonzero velocity
%q0(3:4) = 0.1*rand(2,1);
% observe output
%qact = rr.observe(t,q0,traj.unom);
% add performance to trajectory
%traj.addPerformance(traj.unom,qact,rr.COST,'ID');
% observe output
[qact,ufull] = rr.observeWithFeedbackErrorForm(traj,q0);
traj.addPerformance(ufull,qact,rr.COST,'ID + LQR');

% Plot the controls and animate the robot arm
rr.plot_inputs(traj);
rr.plot_outputs(traj);
%rr.animateArm(qact(1:2,:),ref);

%% Start learning feedforward with ILC

% num_trials = 10;
% 
% %ilc = aILC(rr,traj);
% ilc = mILC(rr,traj); 
% 
% for i = 1:num_trials
%     % get next inputs
%     u = ilc.feedforward(traj,qact);
%     %u = ilc.feedforward(traj,rr,dev);
%     % evolve system
%     qact = rr.observe(t,q0,u);
%     % get the cartesian coordinates
%     %[~,y] = rr.kinematics(qact(1:2,:));
%     % add performance to trajectory
%     traj.addPerformance(u,qact,rr.COST,ilc);
%     % Plot the controls and animate the robot arm
%     %rr.animateArm(qact(1:2,:),ref);
% end
% 
% % Plot the controls and animate the robot arm
% rr.plot_inputs(traj);
% rr.plot_outputs(traj);
% rr.animateArm(qact(1:2,:),ref);

%% Learn with feedback

num_trials = 5;

ilc = mILC(rr,traj); 
%ilc.inp_last = traj.unom;

for i = 1:num_trials
    
    us = ilc.feedforward(traj,qact);     
    traj.unom = us;
    % get the measurements
    [qact,ufull] = rr.observeWithFeedbackErrorForm(traj,q0);
    traj.addPerformance(ufull,qact,rr.COST,ilc);

end

% Plot the controls and animate the robot arm
rr.plot_inputs(traj);
rr.plot_outputs(traj);
%rr.animateArm(qact(1:2,:),ref);

%% Change the reference slightly and see how the robot is doing

% % a new trajectory
% p1 = 0.3;
% p2 = 0.1;
% y_des = 0.4 + p1 * t;
% x_des =  0.6 - p2 * t;
% refNew = [x_des; y_des]; % displacement profile 
% 
% trajNew = rr.generateInputs(t,refNew);
% rr.generateFeedback(trajNew);
% 
% q0 = trajNew.s(:,1);
% % observe output
% qact = rr.observeWithFeedbackErrorForm(trajNew,q0);
% % add performance to trajectory
% trajNew.addPerformance(trajNew.unom,qact,rr.COST,'ID + LQR');
% 
% % Plot the controls and animate the robot arm
% rr.plot_outputs(trajNew);
% rr.animateArm(qact(1:2,:),refNew);
% 
% % test learning starting from scratch
% num_trials = 10;
% ilc = mILC(rr,trajNew);
% 
% for i = 1:num_trials
%     
%     us2 = ilc.feedforward(trajNew,qact);     
%     trajNew.unom = us2;
%     % get the measurements
%     qact = rr.observeWithFeedbackErrorForm(trajNew,q0);
%     trajNew.addPerformance(us2,qact,rr.COST,ilc);
% 
% end
% 
% % Plot the controls and animate the robot arm
% rr.plot_inputs(trajNew);
% rr.plot_outputs(trajNew);
% rr.animateArm(qact(1:2,:),refNew);
% 
% % test with learned inputs of old reference trajectory
% trajNew.unom = us;
% qact = rr.observeWithFeedbackErrorForm(trajNew,q0);
% trajNew.addPerformance(us,qact,rr.COST,'OLD ILC');
% % Plot the controls and animate the robot arm
% rr.plot_outputs(trajNew);
% rr.animateArm(qact(1:2,:),refNew);
% 
% % test learning starting from old inputs
% num_trials = 10;
% ilc = mILC(rr,trajNew);
% for i = 1:num_trials
%     
%     us = ilc.feedforward(trajNew,qact);     
%     trajNew.unom = us;
%     % get the measurements
%     qact = rr.observeWithFeedbackErrorForm(trajNew,q0);
%     trajNew.addPerformance(us,qact,rr.COST,ilc);
% 
% end
% 
% % Plot the controls and animate the robot arm
% rr.plot_inputs(trajNew);
% rr.plot_outputs(trajNew);
% rr.animateArm(qact(1:2,:),refNew);