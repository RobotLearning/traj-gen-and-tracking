%% Table tennis practice using the table tennis class

% original seed was 1
clc; clear; close all; rng(5);
load('ballInitDist1.mat','mu','Sigma');
load('ballTrain1.mat');
% load('ballTrain1');
initializeWAM;

% drawing related params
OPT.draw = true; % draw the simulation
OPT.record = false; % record the simulation

% planning related flags and parameters
OPT.plan.method = 'VHP'; % CENTRED , VHP, LAZY
OPT.plan.vhp.y = -0.6;
OPT.train = false; % train a lookup table using optimization results
OPT.lookup.flag = false; % use lookup table instead of optimizing online
OPT.lookup.mode = 'closest'; %'GP-regress'; %'lin-regress';
OPT.lookup.savefile = 'LookupTable-16-May-2016'; %['LookupTable-', date, '.mat'];

% initial ball parameters
% OPT.distr.type = 'normal';
% OPT.distr.init.mean = mu;
% OPT.distr.init.cov = 0; %Sigma;
% OPT.distr.type = 'empirical';
% OPT.distr.data = ballTrain.bInit;
OPT.distr.type = 'landing';
OPT.distr.init.mean = mu(1:3);
OPT.distr.init.cov = 0; %0.01*Sigma(1:3,1:3);
OPT.distr.land.mean = [-0.2;-1.5];
OPT.distr.land.cov = 0; %0.01*eye(2);
% OPT.distr.type = 'workspace';
% OPT.distr.init.mean = mu(1:3);
% OPT.distr.init.cov = 1e-6;
% OPT.distr.workspace.mean = [-0.5; -0.5; -1.0];
% OPT.distr.workspace.cov = 0.1*eye(3);

% measurement covariance
OPT.camera.cov = 0; %5e-4;
% sampling time
dt = 0.01; 
tt = TableTennis3D(wam,dt,q0,OPT);
tt.practice(q0,5);