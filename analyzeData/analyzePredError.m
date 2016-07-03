%% Analyze prediction error on ball data

% Load real ball data saved on 13/04/2016
% and look at filtering and prediction performance
% as we get more and more ball data
% DATASET 1
% badExamples = [2,19,20];
% in example2 ballgun throws 2 balls
% rebound coeffs dont fit well in 9,11,14,15,
% 25,26,33,38 due to bad estimation/lower spin?
% DATASET 2
% rebound doesnt happen in trial 3,40,53 due to bad ballPred
% and dist_to_table difference
% rebound didnt happen in trial 30,41,45,48,50?
% 10,12,14 has too many outliers in camera1
% 19 doesnt seem to fit well

clc; clear; close all;
loadTennisTableValues; 
dataSet = 1; trial = 1; % Get the data corresponding to trial of interest
[t,B] = loadBallData(dataSet);

% initialize EKF
spin = true;
[filter,ball_func] = initFilterEKF(spin);

% Get reliable ball observations from trial and remove outliers
[b1,b3,ball_est] = pruneBallData(t,B);
[b1,b3,ball_est,numTrials] = getTrialData(b1,b3,trial,dataSet,ball_est);
[t1,b1] = removeOutliers(b1);

% Merge balls
[tMerge,ballMerge] = mergeBallData(t1,b1,b3(:,1),b3(:,2:4));

% Check SL filter results with MATLAB polynomial filter class
ball_est = filterBallsPoly(ball_func,tMerge,ballMerge);

% Get the SL filter ball estimate corresponding to ball trial
[ball_lookup,idx_lookup,t_lookup] = getFilterLookup(ball_est);

% Predict using ball estimate till last blob from camera 1
ballPred = predictTillLastBlob(filter,tMerge,t_lookup,ball_lookup);

% Plot predictions
plotPredictionResult(b1,b3(:,2:4),ballPred);

% Calculate RMS prediction error as we get more ball data till bounce
% update till ball hits table
[~,idx_bounce] = min(b3(:,4));% find the index at bounce
rms_pred = calculatePredErrors(filter,ball_est,idx_lookup,idx_bounce,tMerge,ballMerge);      
figure; plot(rms_pred);