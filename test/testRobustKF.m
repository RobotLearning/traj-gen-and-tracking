%% Test robust filter on actual ball data

clc; clear; close all;
loadTennisTableValues; 
dataSet = 1; 
badExamples = [2,19,20];
numTrials = 53;
goodExamples = setdiff(1:numTrials,badExamples);
trial = 1; iter = 1;% Get the data corresponding to trial of interest
[t,B] = loadBallData(dataSet);

% initialize Kalman Filter with no spin
spin.flag = false;
[ekfFilter,ball_func] = initFilterEKF(spin);

% Just to get the time index to which we will apply the filter
% we want to get the prune the ball data (apply pre-processing/initial
% outlier detection)
% Get reliable ball observations from trial and remove outliers
[B1,B3,~] = pruneBallData(t,B);
[b1,b3,~,~] = getTrialData(B1,B3,trial,dataSet);

% get all of ball data till trial's end
t = t(t >= b1(1,1) & t <= b1(end,1));
cam1isValid = B(1:length(t),2);
cam3isValid = B(1:length(t),7);
B = B(1:length(t),[3:5,8:10]);

%% start filtering

% initialize the EKF filter without initial state estimation
v0 = [rand; 4.0; 2.0];
p0 = [0; dist_to_table-table_length-0.2; table_height+0.2];
x0 = [p0; v0];
P0 = 1e3*eye(length(x0));
ekfFilter.initState(x0,P0);
ekfFilter.linearize(0.01,0);
Q = 1e-3*eye(length(x0));
R = 1e-3*eye(length(3));
ekfFilter.setCov(Q,R);

t_last = 0.0;
diff_t = 0.0;
idxValidBall = 1;
update_flag = true;

for j = 1:length(t)
    
    diff_t = t(j) - t_last;
    
    if cam3isValid(j)
        ekfFilter.predict(diff_t,0);
        update_flag = ekfFilter.robust_update(B(j,4:6)',0);
    elseif cam1isValid(j)
        ekfFilter.predict(diff_t,0);
        update_flag = ekfFilter.robust_update(B(j,1:3)',0);
    end
    
    if update_flag
        ball_est_RKF(idxValidBall,:) = [t(j),ekfFilter.x(:)'];   
        idxValidBall = idxValidBall + 1;
    end
    t_last = t(j);

end

%% to compare we look at a noncausal filter
% that can peak into the future to remove outliers better

[t1,b1] = removeOutliers(b1);
% Merge balls - camera1 offset is also removed
[tMerge,ballMerge,b1] = mergeBallData(t1,b1,b3(:,1),b3(:,2:4));

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
%sOFF = scatter3(ball_est_off(:,2),ball_est_off(:,3),ball_est_off(:,4),'b');
sRKF.MarkerEdgeColor = sRKF.CData;
%sOFF.MarkerEdgeColor = sOFF.CData;
legend('robust filter','offline filter');
fill3(T(1:4,1),T(1:4,2),T(1:4,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
hold off;