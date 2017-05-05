%% Check prediction accuracy of spin model on new data

clc; clear; close all;
loadTennisTableValues;
file = '~/Dropbox/data/balls_0317.txt';
M = dlmread(file);
% ball data
B = M(:,1:10); % cam id and status cols

[B1,B3,~] = prune_ball_data(B);
[balls1,balls3] = get_trials(B1,B3); % cell structure

% plot one to see if it looks sensible
num_trials = length(balls1);
%trial = 3;
%plot_trial(balls1{trial}(:,2:4),balls3{trial}(:,2:4));

% initialize EKF
% initialize ball with high topspin (3000rpm or more)
w0 = [-30*2*pi;0;0]; %3000rpm topspin

for trial = 1:10 %1:num_trials
    fprintf('Trial %d\n',trial);
    [ekfFilter,ball_func] = init_const_spin_filter(w0);

    % remove outliers from camera 1
    balls1{trial} = remove_outliers(balls1{trial});
    t1 = balls1{trial}(:,1);
    b1 = balls1{trial}(:,2:4);
    t3 = balls3{trial}(:,1);
    b3 = balls3{trial}(:,2:4);
    %plot_trial(b1,b3);

    % merge the balls
    [t_m,ball_m,b1] = mergeBallData(t1,b1,t3,b3);

    spin.flag = false; % do not keep count of spin as a separate state
    ball_est = filterBallsEKF(t_m,ball_m,ekfFilter,spin);
    %plot_trial(b1,b3,ball_est(:,2:4));

    % Predict using ball estimate till last blob from camera 1
    t_start = ball_est(1,1);
    ball_start = ball_est(1,2:7);
    ball_pred = predictTillLastBlob(ekfFilter,t_m,t_start,ball_start);
    plot_trial(b1,b3,ball_pred(1:3,:)');

    % Calculate RMS prediction error as we get more ball data till bounce
    % update till ball hits table
    %[~,idx_bounce] = min(b3(:,3));% find the index at bounce
    idx_end = size(b3,1); % end of camera3 data not bounce % idx_bounce; 
    idx_start = 12; % start from first estimate not lookup
    rms_pred{trial} = calcModelPredErrors(ekfFilter,ball_est,...
                              idx_start,idx_end,t_m,ball_m,false);
end

%plotPredErrors(rms_pred);