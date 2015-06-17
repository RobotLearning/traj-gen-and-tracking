%% Simulate trajectories for the Barrett WAM

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
save_folder = '../robolab/barrett/saveData/';
% folder that stores all config files, link parameters, etc.
config_folder = '../robolab/barrett/config/';

%% Define constants and parameters

initializeWAM;

%% Generate inputs for a desired trajectory

% load percentage of trajectory from dmp file 
%file = [save_folder,'dmp_strike.txt'];
%file = 'dmp_strike.txt';
file = 'dmp.txt';
M = dlmread(file);
perc = 1.0; % learning on whole traj can be unstable unless LQR is used
len = size(M,1);
M = M(1:(len * perc),:);
%t = M(:,1); 
t = SIM.h * (1:perc*len);
% order for refs in file: q1 qd1, ...
% switching to order: q1 ... q7, qd1, ..., qd7
q = M(:,2:2:2*N_DOFS);
qd = M(:,3:2:2*N_DOFS+1);
ref = [q';qd'];

% test the kinematics from SL
cart_ref = wam.kinematics(ref);

%[traj,dmp] = wam.generateInputsWithDMP(t,50,ref);
traj = wam.generateInputs(t,ref); % trajectory generated in joint space

% downsample reference
%traj = traj.downsample(10);

% Generate feedback with LQR
%wam.generateFeedback(traj);
% Set initial LQR matrix throughout
%for i = 2:traj.N-1, traj.K(:,:,i) = traj.K(:,:,1); end
% load initial LQR (LQR0)
%load('LQR0.txt','LQR0');
%for i = 1:traj.N-1, FB(:,:,i) = LQR0; end;
% PD control
for i = 1:traj.N-1, FB(:,:,i) = -K; end;
traj.K = FB;

%% Evolve system dynamics and animate the robot arm

%q0 = traj.s(:,1);
q0 = ref(:,1);
%q0(8:end) = 0.0;
% add disturbances around zero velocity
%q0(1:7) = q0(1:7) + 1e-3 * randn(7,1);
%q0(8:end) = 1e-3 * randn(7,1);
% observe output
%qact = wam.evolve(t,q0,traj.unom);
% observe with feedback
[qact,ufull] = wam.observeWithFeedbackErrorForm(traj,q0);
% add performance to trajectory
%traj.addPerformance(traj.unom,qact,wam.COST,'ID');
traj.addPerformance(ufull,qact,wam.COST,'ID + FB');

% Plot the controls and animate the robot arm
%wam.plot_inputs(traj);
wam.plot_outputs(traj);
%wam.animateArm(qact(1:2:2*N_DOFS-1,:),ref);

%% Start learning with ILC

num_trials = 10;
%ilc = aILC(wam,traj,10);
ilc = mILC(wam,traj,2); %downsample 10
%ilc = bILC(traj);

for i = 1:num_trials
    % get next inputs
    %u = ilc.feedforwardQN(traj,qact);
    u = ilc.feedforward(traj,qact);
    % evolve system
    %qact = wam.evolve(t,q0,u);
    % evolve system with feedback
    traj.unom = u;
    [qact,ufull] = wam.observeWithFeedbackErrorForm(traj,q0);
    % add performance to trajectory
    traj.addPerformance(ufull,qact,wam.COST,ilc);
    % Plot the controls and animate the robot arm
    %wam.animateArm(qact(1:2:2*N_DOFS-1,:),ref);
end

% Plot the controls and animate the robot arm
wam.plot_inputs(traj);
wam.plot_outputs(traj);
%wam.animateArm(qact(1:2:2*N_DOFS-1,:),ref);