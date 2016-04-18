%% Train bounce model using new data

clc; clear; close all;

% Detect when it bounces
% Use all data before bounce
% Use nonlinear least squares to optimize flight model parameters
% C and g and also estimate initial ball pos and velocities

%% Load table values

% load table parameters
loadTennisTableValues;

% Load real ball data saved on 13/04/2016
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

%% For each dataset estimate initial b0 values

set1 = 15:65; % dataset of demonstrations
dropSet1 = [2,19,20]; % bad recordings
set = setdiff(1:numTrials,dropSet1);
idxStart3(end+1) = size(b3,1)+1; % for the last dataset

for idx = 1:length(set)
    
    trial = set(idx);
    data = b3(idxStart3(trial):idxStart3(trial+1)-1,2:4);
    time = b3(idxStart3(trial):idxStart3(trial+1)-1,1);
    time = time - time(1);
    [~,idxBounce] = min(data(:,3));
    timePreBounce = time(1:idxBounce);
    ballPreBounce = data(1:idxBounce,:);
    
    data = b1(b1(:,1) >= tStart3(trial) & b1(:,1) < tStart3(trial+1),2:4);
    time = b1(b1(:,1) >= tStart3(trial) & b1(:,1) < tStart3(trial+1),1);
    % remove outliers in b1 after getting ball closest to robot
    [~,idxClosest2Robot] = max(data(:,2));
    data = data(1:idxClosest2Robot,:);
    time = time(1:idxClosest2Robot);
    % outlier detection again
    outlierIdx = detectOutlierBalls(time,data);
    inlierIdx = setdiff(1:length(time),outlierIdx);
    data = data(inlierIdx,:);
    time = time(inlierIdx,:);
    


end