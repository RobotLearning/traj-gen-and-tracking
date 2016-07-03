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
[t,B] = loadBallData(dataSet);