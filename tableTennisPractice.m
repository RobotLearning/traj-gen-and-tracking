%% Table tennis practice using the table tennis class

clc; clear; close all;

rng(1);
load('ballInitDist1.mat','mu','Sigma');
% load('ballTrain1');
initializeWAM;

% drawing related params
OPT.draw = true; % draw the simulation
OPT.record = false; % record the simulation

% planning related flags and parameters
OPT.plan.vhp.flag = false;
OPT.plan.vhp.y = -0.6;
OPT.train = false; % train a lookup table using optimization results
OPT.lookup.flag = false; % use lookup table instead of optimizing online
OPT.lookup.mode = 'regress';
OPT.lookup.savefile = 'LookupTableTest.mat';

% initial ball parameters
% OPT.distr.type = 'normal';
% OPT.distr.init.mean = mu;
% OPT.distr.init.cov = Sigma;
% OPT.distr.type = 'empirical';
% OPT.distr.data = ballTrain.bInit;
OPT.distr.type = 'landing';
OPT.distr.init.mean = mu(1:3);
OPT.distr.init.cov = Sigma(1:3,1:3);
OPT.distr.land.mean = [-0.2;-1.5];
OPT.distr.land.cov = 0.01*eye(2);

% measurement covariance
OPT.camera.cov = 0.0;

tt = TableTennis(wam,q0,OPT);
tt.practice(q0,5);