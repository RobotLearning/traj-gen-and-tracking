%% Table tennis practice using the table tennis class

clc; clear; close all;
load('train/ballInitDist1.mat','mu','Sigma');
initializeWAM;
wam2 = [];
OPT.draw = false; % draw the simulation
OPT.record = false; % record the simulation
OPT.train = true; % train a lookup table using optimization results
OPT.lookup = false; % use lookup table instead of optimizing online
OPT.vhp = false; % use vhp strategy
% initial ball mean
OPT.mean.ballinit = mu;
% initial ball variance
OPT.cov.ballinit = Sigma;
% measurement covariance
OPT.cov.camera = 0.0;
tt = TableTennis(wam,wam2,q0,OPT);
tt.practice(q0,5000);