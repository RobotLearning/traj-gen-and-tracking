%% Analyze prediction error on ball data

% Load real ball data saved on 13/04/2016
% and look at filtering and prediction performance
% as we get more and more ball data

clc; clear; close all;
loadTennisTableValues;
dataSet = 1;
[t,B] = loadBallData(dataSet);

%% initialize EKF

spin = true;
[filter,ball_func] = initFilterEKF(spin);

%% Get reliable ball observations

[b1,b3,filterDataSL] = pruneBallData(t,B);

%% Get the data corresponding to trial of interest

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
trial = 1;
[b1,b3,filtSL,numTrials] = getTrialData(b1,b3,trial,dataSet,filterDataSL);

%% Get the SL filter ball estimate corresponding to ball trial

[ball_lookup,idx_lookup,t_lookup] = getSLFilterLookup(filtSL);
%[ball_lookup_matlab] = get_poly_filter_at_lookup();

%% Remove outliers in b1

[t1,b1] = removeOutliers(b1);

%% Predict using ball estimate till last blob from camera 1

t3 = b3(:,1);
b3 = b3(:,2:4);
ballPred = predictTillLastBlob(filter,t1,t3,b1,b3,t_lookup,ball_lookup);

%% Plot predictions

figure;
s1 = scatter3(b1(:,1),b1(:,2),b1(:,3),'r');
hold on;
s3 = scatter3(b3(:,1),b3(:,2),b3(:,3),'b');
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
        
%% Calculate RMS prediction error as we get more ball data till bounce

% merge balls
ballMerge = mergeBallData(t1,b1,t3,b3);

% update till ball hits table
N_updates = idx_bounce - idx_filter_matlab;
RMS_pred = zeros(N_updates,1);

for i = 1:N_updates
    t_new_lookup = filtSL(idx_lookup+i,1);
    ball_new_lookup = filtSL(idx_lookup+i,2:end);
    [ballPredNew] = predictTillLastBlob(filter,t1,t3,b1,b3,t_new_lookup,ball_new_lookup);
    x_pred = ballPredNew(1:3,:)';
    x_meas = ballMerge(idx_lookup+i:end,:);
    differ = x_pred - x_meas;
    RMS_pred(i) = sqrt(sum(diag((differ)*(differ)'))/N);
    
    filter_poly.update(dt,b3plot(idx_filter_matlab+i,:));
end
        
