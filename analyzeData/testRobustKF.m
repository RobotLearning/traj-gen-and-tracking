%% Test reseting on robust filter for multiple incoming ball trials

clc; clear; close all;
loadTennisTableValues; 
dataSet = 1; 
badExamples = [2,19,20];
numTrials = 53;
goodExamples = setdiff(1:numTrials,badExamples);
range = 3:5; % specify range of trials
trial = range(end); % Get the data corresponding to trial of interest
[t,B] = loadBallData(dataSet);

% initialize Kalman Filter with no spin
spin.flag = false;
[ekfFilter,ball_func] = initFilterEKF(spin);
x0 = ones(6,1);
P0 = 1e3*eye(length(x0));
ekfFilter.initState(x0,P0);
ekfFilter.linearize(0.01,0);

% Just to get the time index to which we will apply the filter
% we want to get the prune the ball data (apply pre-processing/initial
% outlier detection)
% Get reliable ball observations from trial and remove outliers
[B1,B3,~] = pruneBallData(t,B);
[b1,b3,~,~] = getTrialData(B1,B3,range(1),dataSet);
tbegin = b3(1,1);
[b1,b3,~,~] = getTrialData(B1,B3,trial,dataSet);
tend = b1(end,1);

% get all of ball data till trial's end
tIdx = find(t >= tbegin & t <= tend);
t = t(tIdx);
cam1isValid = B(tIdx,2);
cam3isValid = B(tIdx,7);
B = B(tIdx,:);

%% start filtering

t_last = t(1);
diff_t = 0.0;
idxValidBall = 1;
update_flag = false;

for j = 1:length(t)
    
    diff_t = t(j) - t_last;
    ekfFilter.predict(diff_t,0);
    [update_flag,reset_flag] = ekfFilter.robust_update(B(j,:));
    
    if update_flag
        ball_est_RKF(idxValidBall,:) = [t(j),ekfFilter.x(:)'];   
        idxValidBall = idxValidBall + 1;
    end
    
    if reset_flag
        disp('Resetting collected ball data');
        ball_est_RKF = [];
        idxValidBall = 1;
    end
    
    t_last = t(j);

end

%% to compare we look at a noncausal filter
% that can peak into the future to remove outliers better

[t1,b1] = removeOutliers(b1);
% Merge balls - camera1 offset is also removed
[tMerge,ballMerge,b1,offset] = mergeBallData(t1,b1,b3(:,1),b3(:,2:4));

ball_est_off = filterBallsEKF(tMerge,ballMerge,ekfFilter,spin);

%% plot filtering results

figure;
title('Robust filtering of real ball observations. Single trial');
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
hold on;
sRKF = scatter3(ball_est_RKF(:,2),ball_est_RKF(:,3),ball_est_RKF(:,4),'r');
%sOFF = scatter3(b1(:,1),b1(:,2),b1(:,3),'b');
sOFF = scatter3(ball_est_off(:,2),ball_est_off(:,3),ball_est_off(:,4),'b');
sRKF.MarkerEdgeColor = sRKF.CData;
sOFF.MarkerEdgeColor = sOFF.CData;
legend('robust filter','offline filter');
fill3(T(1:4,1),T(1:4,2),T(1:4,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
hold off;