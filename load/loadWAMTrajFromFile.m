%% load percentage of trajectory from dmp file 
% assumes wam has been initialized

%file = [save_folder,'dmp_strike.txt']; % this is 1.5 sec traj
%file = 'dmp_strike.txt'; % about half a second dmp striking trajectory 
file = 'dmp.txt';% another 1.5 sec traj
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

% test the kinematics from SL
cart_ref = wam.kinematics(ref);

% generate u_ff inputs with inverse dynamics
traj = wam.generateInputs(t,ref); % trajectory generated in joint space

% basis functions of DMPs
%bfs = 50;
%[traj,dmp] = wam.generateInputsWithDMP(t,bfs,ref); % trajectory generated in joint space