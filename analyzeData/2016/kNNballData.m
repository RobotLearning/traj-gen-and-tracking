%% Investigate the performance of kNN in prediction on ball data

% Load real ball data saved on 13/04/2016
% and look at prediction performance of kNN
% as we get more and more ball data
% DATASET 1
% badExamples = [2,19,20];
% in example2 ballgun throws 2 balls

clc; clear; close all;
loadTennisTableValues; 
dataSet = 1;
[t_act,B] = loadBallData(dataSet);
% Get reliable ball observations from trial and remove outliers
[ball1,ball3,~] = pruneBallData(t_act,B);
% get number of examples
[~,~,~,numTrials] = getTrialData(ball1,ball3,1,dataSet,[]);

%% Build kNN on DATASET 1
badExamples = [2,19,20];
goodExamples = setdiff(1:numTrials,badExamples);
testSize = floor(length(goodExamples)/5);
testSet = goodExamples(end-testSize+1:end);
trainSet = goodExamples(1:end-testSize);
% time = cell(1,length(trainSet));
% ball = cell(1,length(trainSet));
% cnt = 1;
% for trial = trainSet
%     [b1,b3,~,~] = getTrialData(ball1,ball3,trial,dataSet,[]);
%     [t1,b1] = removeOutliers(b1);
% 
%     % Merge balls - camera1 offset is also removed
%     [tMerge,ballMerge,b1] = mergeBallData(t1,b1,b3(:,1),b3(:,2:4));
%     time{cnt} = tMerge - tMerge(1);
%     ball{cnt} = ballMerge;
%     cnt = cnt + 1;
% end
load('knn_ball_data.mat');
%save('knn_ball_data.mat','time','ball');

%% Test performance of kNN with k = 1

balls_revealed = 20;
idx = 1;
num_neighbor = 5; % k value for knn
for trial = testSet
    [b1,b3,~,~] = getTrialData(ball1,ball3,trial,dataSet,[]);
    [t1,b1] = removeOutliers(b1);

    % Merge balls - camera1 offset is also removed
    [t_act,b_act,b1] = mergeBallData(t1,b1,b3(:,1),b3(:,2:4));
    
    % get only the balls revealed and predict
    [t_pred, ball_pred] = predictKNN(time,ball,...
        t_act(1:balls_revealed),b_act(1:balls_revealed,:),num_neighbor);
                      
    err_rms(idx) = calcPredError(t_pred,ball_pred,t_act,b_act);
    idx = idx + 1;
end

figure; plot(err_rms);