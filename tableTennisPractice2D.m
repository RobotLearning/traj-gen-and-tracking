%% Table tennis practice using the table tennis class

% original seed was 1
clc; clear; close all; rng(2);
load('ballInitDist1.mat','mu','Sigma');
% load('ballTrain1');
initializeRRR;

% drawing related params
OPT.draw = true; % draw the simulation
OPT.record = false; % record the simulation
OPT.rotate = -pi/2;

% planning related flags and parameters
OPT.plan.vhp.flag = false;
OPT.plan.vhp.y = -0.6;
OPT.train = false; % train a lookup table using optimization results
OPT.lookup.flag = false; % use lookup table instead of optimizing online
OPT.lookup.mode = 'closest'; %'regress';
OPT.lookup.savefile = ['R3-LookupTable-', date, '.mat'];

% initial ball parameters
% OPT.distr.type = 'normal';
% OPT.distr.init.mean = mu;
% OPT.distr.init.cov = Sigma;
% OPT.distr.type = 'empirical';
% OPT.distr.data = ballTrain.bInit;
OPT.distr.type = 'landing';
OPT.distr.init.mean = mu(2:3);
OPT.distr.init.cov = 0.01*Sigma(2:3,2:3);
OPT.distr.land.mean = -1.5;
OPT.distr.land.cov = 0.01;
% OPT.distr.type = 'workspace';
% OPT.distr.init.mean = mu(1:3);
% OPT.distr.init.cov = 1e-6;
% OPT.distr.workspace.mean = [-0.5; -0.5; -1.0];
% OPT.distr.workspace.cov = 0.1*eye(3);

% measurement covariance
OPT.vision.cov = 0.0;
OPT.vision.filter = 'EKF';
OPT.vision.draw = true; %plot filtered state

tt = TableTennis2D(rrr,q0,OPT);
tt.practice(q0,1);