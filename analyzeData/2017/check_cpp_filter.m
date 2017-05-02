%% Checking the performance of EKF filter implemented in C++ 

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

% ADD CPP FILTERED BALL VALUES
file_filt = '~/Dropbox/data/realBallData_filtered.txt';
B_filt = dlmread(file_filt);
% get valid ball indices for trial
b_filt = get_valid_balls(B_filt,balls1,balls3,trial);
%plot_trial(b1,b3,b_filt);

% initialize EKF
% initialize ball with high topspin (3000rpm or more)
w0 = [-50*2*pi;0;0]; %3000rpm topspin
spin.flag = false;
%[ekf,ball_func] = initFilterEKF(spin);
[ekf,ball_func] = init_const_spin_filter(w0);

% Filter with MATLAB filter for comparison
net_y = dist_to_table - table_length/2;
balls3_pre_net = balls3{trial}(balls3{trial}(:,3) < net_y,2:4);
t_filt = 0.002 * (1:length(b_filt));
t_net = 0.002 * size(balls3_pre_net,1);
ball_est_matlab = filterBallsEKF(t_filt,balls3_pre_net,ekf,spin);
ball_pred_net = predictTillLastBlob(ekf,t_filt,t_net,...
                            ball_est_matlab(end,2:end));
% % plot predicted balls
ball_preds_net = [ball_est_matlab;ball_pred_net'];
plot_trial(b1,b3,ball_preds_net);

% Plot predicted values at net
% net_y = dist_to_table - table_length/2;
% b_filt_pre_net = b_filt(b_filt(:,2) < net_y,:);
% t_filt = 0.002 * (1:length(b_filt));
% t_net = 0.002 * size(b_filt_pre_net,1);
% ball_pred_net = predictTillLastBlob(ekf,t_filt,t_net,b_filt_pre_net(end,:));
% % plot predicted balls
% ball_preds_net = [b_filt_pre_net;ball_pred_net'];
% plot_trial(b1,b3,ball_preds_net);
% 
% % Plot predicted values around bounce
% [~,idx_bounce] = min(b_filt(:,3));
% b_filt_pre_bounce = b_filt(1:idx_bounce,:);
% t_bounce = 0.002 * size(b_filt_pre_bounce,1);
% ball_pred_bounce = predictTillLastBlob(ekf,t_filt,t_bounce,b_filt_pre_bounce(end,:));
% ball_preds_bounce = [b_filt_pre_bounce;ball_pred_bounce'];
% plot_trial(b1,b3,ball_preds_bounce);