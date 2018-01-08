%% Table tennis practice using the table tennis class

% original seed was 1
clc; clear; close all; rng(1);
std = 0.05;
robot_side = 'RIGHT';
initializeWAM;

% drawing related params
OPT.draw.setup = true; % draw the simulation
OPT.draw.ball.pred.in = false;
OPT.draw.ball.act.out = true;
OPT.draw.robot_traj = true;
OPT.record = false; % record the simulation

% planning related flags and parameters
OPT.plan.method = 'FOCUSED'; % FOCUSED , VHP, DEFENSIVE
OPT.plan.vhp.flag = false;
OPT.plan.vhp.y = -0.6;
OPT.train = true; % train a lookup table using optimization results
OPT.lookup.flag = false; % use lookup table instead of optimizing online
OPT.lookup.mode = 'closest'; %'GP-regress'; %'lin-regress';
OPT.lookup.savefile = 'left_lookup_260617'; %'LookupTable-16-May-2016'; 

% initial ball parameters
OPT.distr.type = 'normal';
OPT.distr.init.std = std;
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
