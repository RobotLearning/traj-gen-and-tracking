%% Investigate filtering and prediction performance for real ball data

% Load real ball data saved on 13/04/2016
% and look at filtering and prediction performance

clc; clear; close all;

% load table parameters
loadTennisTableValues;

y_center = dist_to_table - table_length/2;

dataSet = 1;
file = ['../Desktop/realBallData',int2str(dataSet)];
M = dlmread([file,'.txt']);
% ball data
B = M(1:2:end,:);
% robot joint data
N_DOFS = 7;
Q = M(2:2:end,1:2*N_DOFS);
scale = 0.001; % recorded in milliseconds
t = scale * B(:,11);

%% initialize EKF

params.C = Cdrag;
params.g = gravity;
params.zTable = table_z;
params.yNet = dist_to_table - table_y;
params.table_length = table_length;
params.table_width = table_width;
params.radius = ball_radius;
% coeff of restitution-friction vector
params.CFTX = CFTX;
params.CFTY = CFTY;
params.CRT = CRT;
params.ALG = 'RK4';

funState = @(x,u,dt) discreteBallFlightModel(x,dt,params);
% very small but nonzero value for numerical stability
eps = 1e-6; dim = 6; 
C = [eye(3),zeros(3)];
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);
filter = EKF(dim,funState,mats);

%% Get reliable ball observations
% this is to throw away really bad observations as to not confuse the
% filter
zMax = 0.5;

% camera 1 time indices for reliable observations
idx1 = find(B(:,2) == 1 & B(:,5) >= table_z & B(:,5) < zMax);
t1 = t(idx1);
% order the ball data w.r.t. time and throw away same values
b1 = [t1,B(idx1,3:5)];
b1 = sortrows(b1);
j = 1; % time from ballgun to net 
tol = 1e-3;
idxDiffBallPos1 = []; % index for diff ball positions
for i = 1:size(b1,1)-1
    if norm(b1(i,2:4) - b1(i+1,2:4)) > tol
        idxDiffBallPos1(j) = i;
        j = j+1;
    end
end
b1 = b1(idxDiffBallPos1,:);

% camera 3 time indices for reliable observations
idx3 = find(B(:,7) == 1 & B(:,10) >= table_z & B(:,10) < zMax);
t3 = t(idx3);
% order the ball data w.r.t. time and throw away same values
b3 = [t3,B(idx3,8:10)];
b3 = sortrows(b3);
j = 1; % time from ballgun to net 
tol = 1e-3;
idxDiffBallPos3 = []; % index for diff ball positions
for i = 1:size(b3,1)-1
    if norm(b3(i,2:4) - b3(i+1,2:4)) > tol
        idxDiffBallPos3(j) = i;
        j = j+1;
    end
end
b3 = b3(idxDiffBallPos3,:);

%% Get the number of incoming balls 

% if there is more than 1 second difference it means its a new trial
diffBall3 = diff(b3);
idxStart3 = find(diffBall3(:,1) > 1.0);

% if second dataset then dont keep first trial
if dataSet == 1
    idxStart3 = [1;idxStart3+1];
else
    idxStart3 = idxStart3+1;
end

tStart3 = b3(idxStart3,1);
numTrials = length(idxStart3);

%% Predict ball using estimate (SL uses it for lookup table)

trial = 3;
% DATASET 1
% badExamples = [2,19,20];
% rebound coeffs dont fit well in 9,11,14,15,
% 25,26,33,38 due to bad estimation/lower spin?
% DATASET 2
% rebound doesnt happen in trial 3,40,53 due to bad ballPred
% and dist_to_table difference
% rebound didnt happen in trial 30,41,45,48,50?
% 10,12,14 has too many outliers in camera1
% 19 doesnt seem to fit well

b3est = [t3,B(idx3,12:17)];
b3est = sortrows(b3est);
b3est = b3est(idxDiffBallPos3,:);

ballEst = b3est(idxStart3(trial):idxStart3(trial+1)-1,:);
% remove first 12 entries
ballEst = ballEst(12+1:end,:);

toly = 0.4;

% particular rule implemented in lookup table in SL
idxPred = find(ballEst(:,3) > y_center & ballEst(:,3) < y_center + toly & ...
               ballEst(:,6) > 0.5,1);

ballLookUp = ballEst(idxPred,2:end);
tLookUp = ballEst(idxPred,1);

% get the indices for plotting
b3plot = b3(idxStart3(trial):idxStart3(trial+1)-1,2:4);
t3plot = b3(idxStart3(trial):idxStart3(trial+1)-1,1);
b1plot = b1(b1(:,1) >= tStart3(trial) & b1(:,1) < tStart3(trial+1),2:4);
t1plot = b1(b1(:,1) >= tStart3(trial) & b1(:,1) < tStart3(trial+1),1);

% remove outliers in b1plot after getting ball closest to robot
[~,idxClosest2Robot] = max(b1plot(:,2));
b1plot = b1plot(1:idxClosest2Robot,:);
t1plot = t1plot(1:idxClosest2Robot);

% use ransac to further prune outliers
% outlier detection again
outlierIdx = detectOutlierBalls(t1plot,b1plot);
inlierIdx = setdiff(1:length(t1plot),outlierIdx);
b1plot = b1plot(inlierIdx,:);
t1plot = t1plot(inlierIdx,:);
    

% find the index at bounce
[~,idxBounce] = min(b3plot(:,3));

dtPredTillBounce =  t3plot(idxBounce) - tLookUp;
dtPredTillLastBlob3 = t3plot(end) - tLookUp;
dtPredTillLastBlob1 = t1plot(end) - tLookUp;

% predict using ball estimate
dt = 1/60;
initVar = 1;
filter.initState(ballLookUp(:),initVar);
filter.linearize(dt,0);
predictHorizon = dtPredTillLastBlob3; % only predict till last blob3
table.DIST = dist_to_table;
table.LENGTH = table_length;
table.Z = table_z;
table.WIDTH = table_width;
[ballPred,ballTime,numBounce,time2PassTable] = ...
            predictBallPath(dt,predictHorizon,filter,table);
        
%% Plot predictions

figure;
s1 = scatter3(b1plot(:,1),b1plot(:,2),b1plot(:,3),'r');
hold on;
s3 = scatter3(b3plot(:,1),b3plot(:,2),b3plot(:,3),'b');
%predColor = [0.200 0.200 0.200];
sP = scatter3(ballPred(1,:),ballPred(2,:),ballPred(3,:),'k');

title('Predicted ball trajectory');
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
fill3(T(1:4,1),T(1:4,2),T(1:4,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
s1.MarkerEdgeColor = s1.CData; % due to a bug in MATLAB R2015b
s3.MarkerEdgeColor = s3.CData;
sP.MarkerEdgeColor = sP.CData;
legend('cam1','cam3','filter');
hold off;