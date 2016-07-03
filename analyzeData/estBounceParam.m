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

%% initialize EKF
dim = 6; eps = 1e-6;
C = [eye(3),zeros(3)];

params.C = Cdrag;
params.g = gravity;
params.zTable = table_z;
params.radius = ball_radius;
params.yNet = dist_to_table - table_y;
params.table_length = table_length;
params.table_width = table_width;
% coeff of restitution-friction vector [TO BE LEARNED!]
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

dropSet = [2,19,20,32]; % bad recordings
set = setdiff(1:numTrials,dropSet);
idxStart3(end+1) = size(b3,1)+1; % for the last dataset
tStart3(end+1) = b3(end,1) + 0.01;

for idx = 1:length(set)
    
    %% GET BALLS PRE AND POST BOUNCE
    trial = set(idx);
    data = b3(idxStart3(trial):idxStart3(trial+1)-1,2:4);
    time = b3(idxStart3(trial):idxStart3(trial+1)-1,1);
    time = time - time(1);
    [~,idxBounce] = min(data(:,3));
    timePreBounce = time(1:idxBounce);
    ballPreBounce = data(1:idxBounce,:);
    
    data = b1(b1(:,1) >= tStart3(trial) & b1(:,1) < tStart3(trial+1),2:4);
    time = b1(b1(:,1) >= tStart3(trial) & b1(:,1) < tStart3(trial+1),1);
    time = time - time(1);
    % remove outliers in b1 after getting ball closest to robot
    [~,idxClosest2Robot] = max(data(:,2));
    data = data(1:idxClosest2Robot,:);
    time = time(1:idxClosest2Robot);
    % outlier detection again using ransac
    outlierIdx = detectOutlierBalls(time,data,1);
    inlierIdx = setdiff(1:length(time),outlierIdx);
    data = data(inlierIdx,:);
    time = time(inlierIdx,:);
    
    % get bounce index for balls from camera 1
    idxBallBounce = findReboundIndex(data);
    timePostBounce = time(idxBallBounce:end);
    ballPostBounce = data(idxBallBounce:end,:);
    
    %% RUN KALMAN SMOOTHER TO GET VELOCITIES BEFORE AND AFTER REBOUND
    
    % Provide good initial ball estimates    
    % using polyfit on first 12 balls
    sampleSize = 12; 
    M = [ones(sampleSize,1),timePreBounce(1:sampleSize),timePreBounce(1:sampleSize).^2];
    Y = ballPreBounce(1:sampleSize,:);
    beta = M \ Y;
    ballInitPosEst = beta(1,:);
    ballInitVelEst = beta(2,:);
    ballInit = [ballInitPosEst,ballInitVelEst];
    
    % Run nonlinear least squares to estimate ballInit better
    x0 = ballInit(:);
    ballFun = @(b0,C,g) predictNextBall(b0,C,g,[timePreBounce,ballPreBounce],length(ballPreBounce));
    fnc = @(x) ballFun(x,Cdrag,gravity);
    options = optimoptions('lsqnonlin');
    options.Display = 'final';
    options.Algorithm = 'levenberg-marquardt';
    options.MaxFunEvals = 1500;
    [x,err] = lsqnonlin(fnc,x0,[],[],options);
    ballInit = x(1:6);
    
    % get the pre-rebound velocity using Kalman smoother
    velBeforeBounce(:,idx) = calcReboundVel(ballInit,ballPreBounce,timePreBounce,filter,1);
    
    % provide good initial ball estimate for the 2nd kalman smoother
    timePostBounce = timePostBounce - timePostBounce(1);
    sampleSize = min(12,length(timePostBounce));
    M = [ones(sampleSize,1),timePostBounce(1:sampleSize),timePostBounce(1:sampleSize).^2];
    Y = ballPostBounce(1:sampleSize,:);
    beta = M \ Y;
    ballInitPosEst = beta(1,:);
    ballInitVelEst = beta(2,:);
    ballInit = [ballInitPosEst,ballInitVelEst];
    
    % Run nonlinear least squares to estimate ballInit better
    x0 = ballInit(:);
    ballFun = @(b0,C,g) predictNextBall(b0,C,g,[timePostBounce,ballPostBounce],length(ballPostBounce));
    fnc = @(x) ballFun(x,Cdrag_post,gravity_post);
    [x,err] = lsqnonlin(fnc,x0,[],[],options);
    ballInit = x(1:6);
    
    % there might be a spin change after rebound
    params.C = Cdrag_post;
    params.g = gravity_post;
    funState = @(x,u,dt) discreteBallFlightModel(x,dt,params);
    filter.f = funState;
    velAfterBounce(:,idx) = calcReboundVel(ballInit,ballPostBounce,timePostBounce,filter,0);

end

%% ESTIMATE REBOUND MATRIX

% WITH RANSAC + LS
% beta1 = ransac_ls(velBeforeBounce(1,:)',velAfterBounce(1,:)',5,1000,0.5,0.5);
% beta2 = ransac_ls(velBeforeBounce(2,:)',velAfterBounce(2,:)',5,1000,0.5,0.5);
% beta3 = ransac_ls(velBeforeBounce(3,:)',velAfterBounce(3,:)',5,1000,0.2,0.5);
% beta = ransac_ls(velBeforeBounce',velAfterBounce',10,1000,0.2,0.5);
% beta = beta';

% WITH LLS
% beta1 = velBeforeBounce(1,:)' \ velAfterBounce(1,:)';
% beta2 = velBeforeBounce(2,:)' \ velAfterBounce(2,:)';
% beta3 = velBeforeBounce(3,:)' \ velAfterBounce(3,:)';
% B = diag([beta1,beta2,beta3]);
B = velBeforeBounce' \ velAfterBounce';
B = B';

% LLS WITH A CONSTANT TERM 
% M = [velBeforeBounce',ones(size(velBeforeBounce,2),1)];
% B = M \ velAfterBounce';
% B = B';

%velAfterFit = beta * velBeforeBounce;
velAfterFit = B * velBeforeBounce;
% velAfterFit = B * M';
figure;
title('Velocity after vs. before');
subplot(3,1,1);
scatter(velBeforeBounce(1,:),velAfterBounce(1,:));
hold on;
scatter(velBeforeBounce(1,:),velAfterFit(1,:));
hold off;
subplot(3,1,2);
scatter(velBeforeBounce(2,:),velAfterBounce(2,:));
hold on;
scatter(velBeforeBounce(2,:),velAfterFit(2,:));
hold off;
subplot(3,1,3);
scatter(velBeforeBounce(3,:),velAfterBounce(3,:));
hold on;
scatter(velBeforeBounce(3,:),velAfterFit(3,:));
hold off;

figure;
scatter3(velAfterBounce(1,:),velAfterBounce(2,:),velAfterBounce(3,:));
hold on;
scatter3(velAfterFit(1,:),velAfterFit(2,:),velAfterFit(3,:));