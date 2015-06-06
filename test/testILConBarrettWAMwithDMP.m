%% ILC on BarrettWAM dynamics modifying DMPs now

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

% folder that stores all data
prefs_folder = '../robolab/barrett/prefs/';
% folder that stores all config files, link parameters, etc.
config_folder = '../robolab/barrett/config/';

%% Define constants and parameters

N_DOFS = 7;
% Simulation Values 
% system is continous
SIM.discrete = false;
% learn in cartesian space
SIM.cartesian = false;
% dimension of the x vector
SIM.dimx = 2*N_DOFS;
% dimension of the output y
SIM.dimy = 2*N_DOFS;
% dimension of the control input
SIM.dimu = N_DOFS;
% time step h 
SIM.h = 0.002; % 500 Hz recorded data
% noise and initial error
SIM.eps = 3e-10;
SIM.eps_d = 3e-10;
% integration method
SIM.int = 'Symplectic Euler';
% reference trajectory in joint space?
SIM.jref = true;

% load (nominal) parameter values for robot dynamics
%loadActualBarrettValues;
loadNominalBarrettValues;

% We observe all joints and all joint velocities
PAR.links = links;
PAR.link0 = link0;
PAR.eff = eff;
PAR.basec = basec;
PAR.baseo = baseo;
PAR.uex = uex;
PAR.uex0 = uex0;
PAR.C = eye(SIM.dimy,SIM.dimx);

% form constraints
CON = [];

% Set FB gains yourself - e.g. PD control
% loading from the gains.cf file 
K(1,1) = 200.0;
K(1,8) = 7.0;
K(2,2) = 300.0;
K(2,9) = 15.0;
K(3,3) = 100.0;
K(3,10) = 5.0;
K(4,4) = 50.0;
K(4,11) = 2.5;
K(5,5) = 10.0;
K(5,12) = 0.3;
K(6,6) = 10.0;
K(6,13) = 0.3;
K(7,7) = 2.5;
K(7,14) = 0.075;
% TODO: load also u_max limits!

% cost structure
% only penalize positions
Q1 = 1*diag([ones(1,4),1*ones(1,3),1*ones(1,4),1*ones(1,3)]);
Q2 = 1*diag([ones(1,4),1*ones(1,3),0.1*ones(1,4),0.1*ones(1,3)]);
COST.Q = Q1;
COST.R = 0.01 * eye(SIM.dimu);

% initialize model
wam = BarrettWAM(PAR,CON,COST,SIM);

%% Generate inputs for a desired trajectory

% load percentage of trajectory from dmp file 
%file = [prefs_folder,'dmp_strike.txt'];
file = 'dmp.txt';
M = dlmread(file);
perc = 1.0; % learning on whole traj can be unstable unless LQR is used
len = size(M,1);
M = M(1:(len * perc),:);
t = SIM.h * (1:perc*len);
% order for refs in file: q1 qd1, ...
% switching to order: q1 ... q7, qd1, ..., qd7
q = M(:,2:2:2*N_DOFS);
qd = M(:,3:2:2*N_DOFS+1);
ref = [q';qd'];

% generate u_ff inputs with inverse dynamics
%traj = wam.generateInputs(t,ref); % trajectory generated in joint space

% basis functions of DMPs
bfs = 50;
[traj,dmp] = wam.generateInputsWithDMP(t,bfs,ref); % trajectory generated in joint space

% downsample reference
%traj = traj.downsample(10);

% Generate feedback with LQR
%wam.generateFeedback(traj);
% Load feedback in case trajectory is very large
%load('LQR.mat','FB');
% load initial LQR (LQR0)
load('LQR0.txt','LQR0');
for i = 1:traj.N-1, FB(:,:,i) = LQR0; end;
% PD control
%for i = 1:traj.N-1, FB(:,:,i) = -K; end;
%traj.K = FB;

%% Evolve system dynamics and animate the robot arm

q0 = traj.s(:,1);
% add zero velocity as disturbance
%q0(N_DOFS+1:end) = 0;
% observe output
%qact = wam.evolve(t,q0,traj.unom);
% observe with feedback
[qact,ufull] = wam.observeWithFeedbackErrorForm(traj,q0);
% add performance to trajectory
%traj.addPerformance(traj.unom,qact,wam.COST,'Inverse Dynamics');
traj.addPerformance(ufull,qact,wam.COST,'ID + FB');

% Plot the controls and animate the robot arm
wam.plot_inputs(traj);
wam.plot_outputs(traj);
%wam.animateArm(qact(1:2:2*N_DOFS-1,:),ref);
trajNew = Trajectory(traj.t,traj.s,traj.unom,traj.K);

%% Start learning with ILC

num_trials = 10;
ilc = mILC(wam,traj,10);
uadd = zeros(size(traj.unom));

for i = 1:num_trials
    % adapt the dmps accordingly
    % get the next inputs normally as in standard ILC
    uadd = ilc.feedforwardDMP(trajNew,qact,uadd);
    % change initial condition slightly
    q0new = q0 + 0.1 * randn(length(q0),1);
    trajModified = wam.generateInputsWithDMP(t,bfs,ref,q0new);
    trajNew = Trajectory(traj.t,trajModified.s,trajModified.unom,traj.K);
    % adjust for the IDM change
    trajNew.unom = uadd + trajNew.unom;
    [qact,ufull] = wam.observeWithFeedbackErrorForm(trajNew,q0new);
    % current iteration part of ILC
    uadd = ufull - trajNew.unom + uadd;
    traj.addPerformance(ufull,qact,wam.COST,ilc);
    % Plot the controls and animate the robot arm
    %wam.animateArm(qact(1:2:2*N_DOFS-1,:),ref);
end

% Plot the controls and animate the robot arm
wam.plot_inputs(trajNew);
wam.plot_outputs(trajNew);
%wam.animateArm(qact(1:2:2*N_DOFS-1,:),ref);