%% Table tennis practice using the table tennis class

clc; clear; close all;

% load('ballInitDist1.mat','mu','Sigma');
load('ballTrain1');
initializeWAM;

% drawing related params
OPT.draw = true; % draw the simulation
OPT.record = false; % record the simulation

% planning related flags and parameters
OPT.plan.vhp.flag = false;
OPT.plan.vhp.y = -0.6;
OPT.train = false; % train a lookup table using optimization results
OPT.lookup.flag = false; % use lookup table instead of optimizing online
OPT.lookup.savefile = 'LookupTable.mat';

% initial ball parameters
% OPT.distr.type = 'normal';
% OPT.distr.mean = mu;
% OPT.distr.cov = Sigma;
OPT.distr.type = 'empirical';
OPT.distr.data = ballTrain.bInit;

% measurement covariance
OPT.camera.cov = 0.0;

tt = TableTennis(wam,q0,OPT);
tt.practice(q0,5);