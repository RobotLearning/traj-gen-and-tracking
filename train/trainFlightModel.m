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
    time = time(1:idxBounce);
    data = data(1:idxBounce,:);
    ballData = [time,data];
    
    % Provide good initial ball estimates    
    % using polyfit on first 12 balls
    sampleSize = 12; 
    M = [ones(sampleSize,1),time(1:sampleSize),time(1:sampleSize).^2];
    Y = data(1:sampleSize,:);
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
    ballTrain.bInit(:,idx) = ballInit(:);
    ballTrain.data{idx} = ballData;
    ballTrain.error{idx} = err;
end

save('train/ballTrain.mat','ballTrain');

%% Nonlinear least squares to estimate flight model parameters

b0 = ballTrain.bInit(:);
x0 = [b0; Cdrag; gravity];
data = [];
for i = 1:length(set)
    data = [data;ballTrain.data{idx}];
    lenSet(i) = size(ballTrain.data{idx},1);
end
ballFun = @(b0,C,g) predictNextBall(b0,C,g,data,lenSet);
fun = @(x) ballFun(x(1:end-2),x(end-1),x(end));
options = optimoptions('lsqnonlin','Display','iter','MaxFunEvals',3000);
x = lsqnonlin(fun,x0,[],[],options);
C = x(end-1)
g = x(end)
%}