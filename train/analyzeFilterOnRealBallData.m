%% Investigate filtering and prediction performance for real ball data

% Load real ball data saved on 13/04/2016
% and look at filtering and prediction performance

clc; clear; close all;

% load table parameters
loadTennisTableValues;
% Initialize Barrett WAM
initializeWAM;

yCenter = dist_to_table - table_length/2;

file = '../Desktop/realBallData1';
M = dlmread([file,'.txt']);
% ball data
B = M(1:2:end,:);
% robot joint data
Q = M(2:2:end,1:2*N_DOFS);
scale = 0.001; % recorded in milliseconds
t = scale * B(:,11);

%% initialize EKF
dim = 6;
eps = 1e-6;
C = [eye(3),zeros(3)];

Cdrag = 0.1476;
gravity = -10.37;

params.C = Cdrag;
params.g = gravity;
params.zTable = table_z;
params.yNet = dist_to_table - table_y;
params.table_length = table_length;
params.table_width = table_width;
% coeff of restitution-friction vector
params.CFTX = CFTX;
params.CFTY = CFTY;
params.CRT = CRT;
params.ALG = 'RK4';

funState = @(x,u,dt) discreteBallFlightModel(x,dt,params);
% very small but nonzero value for numerical stability
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);
filter = EKF(dim,funState,mats);

%% Get the number of incoming balls 
% this is to throw away really bad observations as to not confuse the
% filter
zMax = 0.5;
% camera 1 time indices for reliable observations
idx1 = find(B(:,2) == 1 & B(:,5) >= table_z & B(:,5) < zMax);
t1 = t(idx1);
Mball1 = B(:,3:5);
bx1 = Mball1(idx1,1);
by1 = Mball1(idx1,2);
bz1 = Mball1(idx1,3);

% camera 3 time indices for reliable observations
idx3 = find(B(:,7) == 1 & B(:,10) >= table_z & B(:,10) < zMax);
t3 = t(idx3);

% order the ball data w.r.t. time and throw away same values
b3 = [t3,B(idx3,8:10)];
b3 = sortrows(b3);

j = 1; % time from ballgun to net 
tol = 1e-3;
idxDiffBallPos = []; % index for diff ball positions
for i = 1:size(b3,1)-1
    if norm(b3(i,2:4) - b3(i+1,2:4)) > tol
        idxDiffBallPos(j) = i;
        j = j+1;
    end
end

b3 = b3(idxDiffBallPos,:);

% if there is more than 1 second difference it means its a new trial
diffBall3 = diff(b3);
idxStart = find(diffBall3(:,1) > 1.0);
idxStart = [1;idxStart+1];
tStart = b3(diffBall3(:,1) > 1.0,1);
numTrials = length(idxStart);

%% Predict ball using estimate (SL uses it for lookup table)

trial = 1;

b3Pred = [t3,B(idx3,12:17)];
b3Pred = sortrows(b3Pred);
b3Pred = b3Pred(idxDiffBallPos,:);

ballPred = b3Pred(idxStart(trial):idxStart(trial+1)-1,:);
% remove first 12 entries
ballPred = ballPred(12+1:end,:);

toly = 0.4;

% particular rule implemented in lookup table in SL
idxPred = find(ballPred(:,3) > yCenter & ballPred(:,3) < yCenter + toly & ...
               ballPred(:,6) > 0.5,1);

ballLookUp = ballPred(idxPred,2:end);

% predict using ball estimate
dt = 0.01;
initVar = 1;
filter.initState(ballLookUp(:),initVar);
filter.linearize(dt,0);
predictHorizon = 0.8;
table.DIST = dist_to_table;
table.LENGTH = table_length;
table.Z = table_z;
table.WIDTH = table_width;
[ballPred,ballTime,numBounce,time2PassTable] = ...
            predictBallPath(dt,predictHorizon,filter,table);