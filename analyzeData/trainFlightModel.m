%% Train flight model using new data

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

numTrials = length(idxStart3);

%% For each dataset estimate initial b0 values

set1 = 15:65; % dataset of demonstrations
dropSet1 = [2,19,20]; % bad recordings
set = setdiff(1:numTrials,dropSet1);
tStart3 = b3(idxStart3,1);
idxStart3(end+1) = size(b3,1)+1; % for the last dataset
tStart3(end+1) = b3(end,1) + 0.01;

for idx = 1:length(set)
    
    %% PRE BOUNCE ESTIMATION
    trial = set(idx);
    ballPreBounce = b3(idxStart3(trial):idxStart3(trial+1)-1,2:4);
    timePreBounce = b3(idxStart3(trial):idxStart3(trial+1)-1,1);
    timePreBounce = timePreBounce - timePreBounce(1);
    [~,idxBounce] = min(ballPreBounce(:,3));
    timePreBounce = timePreBounce(1:idxBounce);
    ballPreBounce = ballPreBounce(1:idxBounce,:);
    ballData = [timePreBounce,ballPreBounce];
    
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
    ballFun = @(b0,C,g) predictNextBall(b0,C,g,ballData,length(ballData));
    fnc = @(x) ballFun(x,Cdrag,gravity);
    options = optimoptions('lsqnonlin');
    options.Display = 'final';
    options.Algorithm = 'levenberg-marquardt';
    options.MaxFunEvals = 1500;
    [x,err] = lsqnonlin(fnc,x0,[],[],options);
    ballInit = x(1:6);
    
    % data used for optimization is saved here
    ball.pre.bInit(:,idx) = ballInit(:);
    ball.pre.data{idx} = ballData;
    ball.pre.error{idx} = err;
    
    %% POST BOUNCE PARAMETER ESTIMATION
    
    ballPostBounce = b1(b1(:,1) >= tStart3(trial) & b1(:,1) < tStart3(trial+1),2:4);
    timePostBounce = b1(b1(:,1) >= tStart3(trial) & b1(:,1) < tStart3(trial+1),1);
    timePostBounce = timePostBounce - timePostBounce(1);

    % remove outliers in b1plot after getting ball closest to robot
    [~,idxClosest2Robot] = max(ballPostBounce(:,2));
    ballPostBounce = ballPostBounce(1:idxClosest2Robot,:);
    timePostBounce = timePostBounce(1:idxClosest2Robot);
 
    % use ransac to further prune outliers
    % outlier detection again
    outlierIdx = detectOutlierBalls(timePostBounce,ballPostBounce,1);
    inlierIdx = setdiff(1:length(timePostBounce),outlierIdx);
    ballPostBounce = ballPostBounce(inlierIdx,:);
    timePostBounce = timePostBounce(inlierIdx);
    
    % get the post bounce balls
    idxBounce = findReboundIndex(ballPostBounce);
    ballPostBounce = ballPostBounce(idxBounce+1:end,:);
    timePostBounce = timePostBounce(idxBounce+1:end);
    ballData = [timePostBounce,ballPostBounce];
    
    % Provide good initial ball estimates    
    % using polyfit on first 12 balls
    sampleSize = min(length(timePostBounce),12); 
    M = [ones(sampleSize,1),timePostBounce(1:sampleSize),timePostBounce(1:sampleSize).^2];
    Y = ballPostBounce(1:sampleSize,:);
    beta = M \ Y;
    ballInitPosEst = beta(1,:);
    ballInitVelEst = beta(2,:);
    ballInit = [ballInitPosEst,ballInitVelEst];
    
    % Run nonlinear least squares to estimate ballInit better
    x0 = ballInit(:);
    ballFun = @(b0,C,g) predictNextBall(b0,C,g,ballData,length(ballData));
    fnc = @(x) ballFun(x,Cdrag,gravity);
    [x,err] = lsqnonlin(fnc,x0,[],[],options);
    ballInit = x(1:6);
    
    % data used for optimization is saved here
    ball.post.bInit(:,idx) = ballInit(:);
    ball.post.data{idx} = ballData;
    ball.post.error{idx} = err;
end

%% Nonlinear least squares to estimate flight model parameters

b0 = ball.pre.bInit(:);
x0 = [b0; Cdrag; gravity];
ballPreBounce = [];
for i = 1:length(set)
    ballPreBounce = [ballPreBounce;ball.pre.data{idx}];
    lenSet(i) = size(ball.pre.data{idx},1);
end
ballFun = @(b0,C,g) predictNextBall(b0,C,g,ballPreBounce,lenSet);
fun = @(x) ballFun(x(1:end-2),x(end-1),x(end));
options = optimoptions('lsqnonlin','Display','iter','MaxFunEvals',3000);
x = lsqnonlin(fun,x0,[],[],options);
C_pre = x(end-1)
g_pre = x(end)

b0 = ball.post.bInit(:);
x0 = x;
% x0 = [b0; C_pre; g_pre];
%x0 = [x(1:end-2); Cdrag; gravity];
ballPostBounce = [];
for i = 1:length(set)
    ballPostBounce = [ballPostBounce;ball.post.data{idx}];
    lenSet(i) = size(ball.post.data{idx},1);
end
ballFun = @(b0,C,g) predictNextBall(b0,C,g,ballPostBounce,lenSet);
fun = @(x) ballFun(x(1:end-2),x(end-1),x(end));
x = lsqnonlin(fun,x0,[],[],options);
C_post = x(end-1)
g_post = x(end)

%}

%% Estimate initial distribution as a gaussian

ballInit = ballPreBounce.bInit;
N = size(ballInit,2);
mu_hat = sum(ballInit,2)/N;
res = ballInit - repmat(mu_hat,1,N);
Sigma_hat = res * res' / (N-1);

save('train/ballTrain1.mat','ball');