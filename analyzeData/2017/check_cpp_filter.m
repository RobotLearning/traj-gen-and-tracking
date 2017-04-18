%% Checking the performance of EKF filter impemented in C++ 

clc; clear; close all;
loadTennisTableValues;
file = '~/Dropbox/data/realBallData_030516.txt';
M = dlmread(file);
% ball data
B = M(:,1:10); % cam id and status cols

[B1,B3,~] = prune_ball_data(B);
[balls1,balls3] = get_trials(B1,B3); % cell structure

trial = 1;
% remove outliers from camera 1
balls1{trial} = remove_outliers(balls1{trial});
t1 = balls1{trial}(:,1);
b1 = balls1{trial}(:,2:4);
t3 = balls3{trial}(:,1);
b3 = balls3{trial}(:,2:4);
plot_trial(b1,b3);

% ADD CPP FILTERED BALL VALUES
file_filt = '~/Dropbox/data/realBallData_filtered.txt';
B_filt = dlmread(file_filt);
B_filt = process_filter_values(trial);

