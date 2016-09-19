%% Test reseting on robust filter for multiple incoming ball trials

clc; clear; close all;
loadTennisTableValues; 
dataSet = 1; 
badExamples = [2,19,20];
numTrials = 53;
goodExamples = setdiff(1:numTrials,badExamples);
range = 1; % specify range of trials
trial = range(end); % Get the data corresponding to trial of interest
[t,B] = loadBallData(dataSet);

% initialize Kalman Filter with no spin
w0 = [-50*2*pi;0;0]; %3000rpm topspin
spin.flag = false;
spin.known = true;
spin.Clift = Clift;
spin.est = w0;
[ekfFilter,ball_func] = initFilterEKF(spin);

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
    ekfFilter.robust_update(diff_t,B(j,:));
    
    if ekfFilter.update_flag
        ball_est_RKF(idxValidBall,:) = [t(j),ekfFilter.x(:)'];
        idxValidBall = idxValidBall + 1;
        t_last = t(j);
    end
    
    if ekfFilter.reset_flag && ekfFilter.update_flag && size(ekfFilter.reset_balls,2) == 1
        disp('Resetting collected ball data');
        ball_est_RKF = [];
        idxValidBall = 1;
        t_last = t_last - 0.014;
    end
    

end

%% to compare we look at a noncausal filter
% that can peak into the future to remove outliers better
% 
% [t1,b1] = removeOutliers(b1);
% % Merge balls - camera1 offset is also removed
% [tMerge,ballMerge,b1,offset] = mergeBallData(t1,b1,b3(:,1),b3(:,2:4));
% 
% ball_est_off = filterBallsEKF(tMerge,ballMerge,ekfFilter,spin);

%% look at prediction performance

startPredIdx = 20;
tstart = ball_est_RKF(startPredIdx,1);
ballStart = ball_est_RKF(startPredIdx,2:end);
tPred = t(t > tstart & t < ball_est_RKF(end,1));
tPred = tPred(1):0.016:tPred(end);
ballPred = predictTillLastBlob(ekfFilter,tPred,tstart,ballStart);

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
%sOFF = scatter3(ball_est_off(:,2),ball_est_off(:,3),ball_est_off(:,4),'b');
sPred = scatter3(ballPred(1,:),ballPred(2,:),ballPred(3,:),'b');
sRKF.MarkerEdgeColor = sRKF.CData;
%sOFF.MarkerEdgeColor = sOFF.CData;
legend('robust filter','pred'); %'offline filter');
fill3(T(1:4,1),T(1:4,2),T(1:4,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
hold off;