%% Simulate trajectories for the Barrett WAM

clc; clear; close all;
robot_side = 'RIGHT';
% Define constants and parameters
initializeWAM;
% load initial LQR
LQR = load('LQR0.txt');
% load initial and final desired positions
q0 = load('q0.txt');
qf = load('qf.txt');
W = load('w_strike.txt');

%% Generate inputs for a desired trajectory

%{
tfin = 0.8;
h = SIM.h;
t = h:h:tfin;
N = length(t);
% canonical system
tau = 1/t(end);
alpha = 25;
beta = alpha/4;
ax = 1;
% number of basis functions
numbf = 50;
pat = 'd';
can = CAN(h,ax,tau,numbf,tfin,pat);

for i = 1:7
    dmp(i) = DDMP(can,alpha,beta,qf(i),[q0(i);0;0]);
    dmp(i).w = W(:,i);
end

% Generate inputs for DMP 
traj = wam.generateInputsForDMP(dmp,N);
%}
loadWAMTrajFromFile; 

% Generate feedback with LQR
%wam.generateFeedback(traj);
% Set initial LQR matrix throughout
%for i = 2:traj.N-1, traj.K(:,:,i) = traj.K(:,:,1); end
% For controlling only based on positions
% lqrPos = LQR; %traj.K(:,:,1); 
% lqrPos(:,8:end) = 0.0;
% for i = 1:traj.N-1, traj.K(:,:,i) = lqrPos; end
%for i = 1:traj.N-1, traj.K(:,:,i) = LQR; end;
% PD control
for i = 1:traj.N-1, traj.K(:,:,i) = PD; end;

%% Evolve system dynamics and animate the robot arm

% when loading traj from file
q0 = traj.s(:,1);
%q0 = ref(:,1);
q0(8:14) = 0.0;
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
wam.plot_inputs(traj);
wam.plot_outputs(traj);
%wam.animateArm(qact(1:2:2*N_DOFS-1,:),ref);

%% Start learning with ILC

num_trials = 5;
%ilc = aILC(wam,traj,10);
%ilc = mILC(wam,traj,10); %downsample 10
ilc = bILC(traj,wam);

for i = 1:num_trials
    % get next inputs
    u = ilc.feedforward(traj,qact);
    % evolve system
    %qact = wam.evolve(t,q0,u);
    % evolve system with feedback
    traj.unom = u;
    [qact,ufull] = wam.observeWithFeedbackErrorForm(traj,q0);
    % add performance to trajectory
    %traj.addPerformance(u,qact,wam.COST,ilc);
    traj.addPerformance(ufull,qact,wam.COST,ilc);
    % Plot the controls and animate the robot arm
    %wam.animateArm(qact(1:2:2*N_DOFS-1,:),ref);
end

% Plot the controls and animate the robot arm
wam.plot_inputs(traj);
wam.plot_outputs(traj);
%wam.animateArm(qact(1:2:2*N_DOFS-1,:),ref);