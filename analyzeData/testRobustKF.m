%% Test robust filter on actual ball data

clc; clear; close all;
loadTennisTableValues; 
dataSet = 1; 
badExamples = [2,19,20];
numTrials = 53;
goodExamples = setdiff(1:numTrials,badExamples);
trial = 6; iter = 1;% Get the data corresponding to trial of interest
[t,B] = loadBallData(dataSet);

% initialize Kalman Filter with no spin
spin.flag = false;
[ekfFilter,ball_func] = initFilterEKF(spin);

% Just to get the time index to which we will apply the filter
% we want to get the prune the ball data (apply pre-processing/initial
% outlier detection)
% Get reliable ball observations from trial and remove outliers
[B1,B3,Ball_est] = pruneBallData(t,B);
[b1,b3,~,~] = getTrialData(B1,B3,trial,dataSet);

% get all of ball data till trial's end
t = t(t < b1(end,1));
cam1isValid = B(1:length(t),2);
cam3isValid = B(1:length(t),7);
B = B(1:length(t),[3:5,8:10]);

%% start filtering

% initialize the EKF filter without initial state estimation
v0 = [rand; 4.0; 2.0];
p0 = [0; 0;0]; %dist_to_table-table_length-0.2; table_height+0.2];
x0 = [p0; v0];
P0 = 1e3*eye(length(x0));
ekfFilter.initState(x0,P0);
ekfFilter.linearize(0.01,0);
Q = eye(length(x0));
R = 1e3*eye(length(3));
ekfFilter.setCov(Q,R);

t_last = 0.0;
idxValidBall = 1;
for j = 1:length(t)
    diff_t = t(j) - t_last;
    %filter.linearize(dt,0);
    ekfFilter.predict(diff_t,0);
    
    if cam3isValid(j)
        ekfFilter.robust_update(B(j,4:6)',0);
        ballEsts(idxValidBall,:) = [t(j),ekfFilter.x(:)'];   
        idxValidBall = idxValidBall + 1;
    %elseif cam1isValid(j)
    %    ekfFilter.robust_update(B(j,1:3)',0);
    %    ballEsts(idxValidBall,:) = [t(j),ekfFilter.x(:)'];   
    %    idxValidBall = idxValidBall + 1;
    end
    
    t_last = t(j);
end

%% plot filtering results

figure;
plot3(ballEsts(:,2),ballEsts(:,3),ballEsts(:,4));