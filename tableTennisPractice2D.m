%% Table tennis practice using the table tennis class

% TODO: 
% Ball should reflect only when touching the racket. 
%    check contact model
% Check for torque limits as well in optimization
% Check for table in optimization
% Build local policy

% original seed was 2 for lookup table
clc; clear; close all; rng(4);
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
OPT.lookup.flag = true; % use lookup table instead of optimizing online
OPT.lookup.mode = 'local-policy'; %'regress'; %'closest';
OPT.lookup.savefile = 'R3-LookupTable-22-May-2016.mat'; %, date, '.mat'];

% initial ball parameters
% OPT.distr.type = 'normal';
% OPT.distr.init.mean = mu;
% OPT.distr.init.cov = Sigma;
% OPT.distr.type = 'empirical';
% OPT.distr.data = ballTrain.bInit;
OPT.distr.type = 'landing';
OPT.distr.init.mean = mu(2:3);
OPT.distr.init.cov = Sigma(2:3,2:3);
OPT.distr.land.mean = -1.85;
OPT.distr.land.cov = 0.05;
OPT.distr.time.mean = 0.5; %time to land
OPT.distr.time.cov = 0.01;

% measurement covariance
OPT.vision.cov = 0.0;
OPT.vision.filter = 'EKF';
OPT.vision.draw = true; %plot filtered state

tt = TableTennis2D(rrr,q0,OPT);
tt.practice(q0,5);