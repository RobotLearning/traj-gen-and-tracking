%% Table tennis practice using the table tennis class

clc; clear; close all;

load('ballInitDist1.mat','mu','Sigma');
initializeWAM;

OPT.draw = true; % draw the simulation
OPT.record = false; % record the simulation
OPT.train = false; % train a lookup table using optimization results
OPT.lookup.flag = false; % use lookup table instead of optimizing online
OPT.lookup.savefile = 'LookupTable.mat';

% initial ball parameters
OPT.distr.type = 'normal';
OPT.distr.mean = mu;
OPT.distr.cov = Sigma;
% measurement covariance
OPT.cov.camera = 0.0;

tt = TableTennis(wam,q0,OPT);
tt.practice(q0,1);