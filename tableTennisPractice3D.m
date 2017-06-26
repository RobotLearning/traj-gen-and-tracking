%% Table tennis practice using the table tennis class

% original seed was 1
clc; clear; close all; %rng(5);
std = 0.05;
robot_side = 'CENTRE';
ballgun_side = 'CENTRE';
initializeWAM;
[b0, v0] = initBallGun(std,ballgun_side);

% drawing related params
OPT.draw = true; % draw the simulation
OPT.record = false; % record the simulation

% planning related flags and parameters
OPT.plan.method = 'LAZY'; % FOCUSED , VHP, LAZY
OPT.plan.vhp.y = -0.6;
OPT.train = false; % train a lookup table using optimization results
OPT.lookup.flag = false; % use lookup table instead of optimizing online
OPT.lookup.mode = 'closest'; %'GP-regress'; %'lin-regress';
OPT.lookup.savefile = 'LookupTable-16-May-2016'; %['LookupTable-', date, '.mat'];

% initial ball parameters
OPT.distr.type = 'normal';
OPT.distr.init.mean = [b0; v0];
OPT.distr.init.cov = std^2; 
%{
% OPT.distr.type = 'empirical';
% OPT.distr.data = ballTrain.bInit;
% OPT.distr.type = 'landing';
% OPT.distr.init.mean = 0.3 * randn(3,1) + mu(1:3);
% OPT.distr.init.cov = Sigma(1:3,1:3);
% OPT.distr.land.mean = [-0.2;-1.5];
% OPT.distr.land.scale = 0.4; % centered uniform distr 
%}

% measurement covariance
OPT.camera.cov = 0.0; %5e-4;
% sampling time
dt = 0.01; 
tt = TableTennis3D(wam,dt,q0,OPT);
tt.practice(q0,1);