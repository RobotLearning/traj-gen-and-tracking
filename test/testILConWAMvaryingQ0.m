%% ILC on BarrettWAM dynamics modifying discrete DMPs now

%clc; clear; close all;
% TEST SCRIPT FOR HUMANOIDS 2015

%% Initialize WAM: define constants and parameters

initializeWAM;

%% Generate ff inputs and feedback for a desired trajectory

if firsttime
    firsttime = 0;
    loadWAMTrajFromFile; 
end

% downsample reference
%traj = traj.downsample(10);

% Generate feedback with LQR
%wam.generateFeedback(traj);
% Load feedback in case trajectory is very large
% Set initial LQR matrix throughout
%for i = 2:traj.N-1, traj.K(:,:,i) = traj.K(:,:,1); end
%load('LQR.mat','FB');
% load initial LQR (LQR0)
load('LQR0.txt','LQR0');
for i = 1:traj.N-1, FB(:,:,i) = LQR0; end;
traj = wam.generateInputsForDMP(dmp,traj.N);
traj.K = FB;

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
%wam.plot_inputs(traj);
%wam.plot_outputs(traj);
%wam.animateArm(qact(1:2:2*N_DOFS-1,:),ref);

%% Start learning with ILC

num_trials = 5;
ilcClass = 'gILC';
switch ilcClass
    
    case 'mILC'
        ilc = mILC(wam,traj,10,1);
    case 'fbILC'
        ilc = mILC(wam,traj,10,0);
    case 'gILC'
        ilc = mILC(wam,traj,10,2);
    
end

uadd = zeros(size(traj.unom));

for i = 1:num_trials
    % change initial condition slightly
    [q0,e0] = perturbInitCondDMP(dmp,0.1);
    traj = wam.generateInputsForDMP(dmp,traj.N);
    traj.K = FB;
    % adapt the dmps accordingly
    % get the next inputs normally as in standard ILC
    uadd = ilc.feedforwardDMP(traj,qact,uadd,e0);
    % adjust for the IDM change
    traj.unom = uadd + traj.unom;
    [qact,ufull] = wam.observeWithFeedbackErrorForm(traj,q0);
    % current iteration part of ILC
    uadd = ufull - traj.unom + uadd;
    traj.addPerformance(ufull,qact,wam.COST,ilc);
    % Plot the controls and animate the robot arm
    %wam.animateArm(qact(1:2:2*N_DOFS-1,:),ref);
end

% Plot the controls and animate the robot arm
%wam.plot_inputs(traj);
%wam.plot_outputs(traj);
%wam.animateArm(qact(1:2:2*N_DOFS-1,:),ref);