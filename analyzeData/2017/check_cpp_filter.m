%% Checking the performance of EKF filter implemented in C++ 

clc; clear; close all;
loadTennisTableValues;
file = '~/Dropbox/data/real_ball_data_0317/balls_30.txt';  
%'~/Dropbox/data/realBallData_030516.txt';
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
file_filt = '~/Dropbox/data/real_ball_data_0317/balls_30_filtered.txt';
%'~/Dropbox/data/realBallData_filtered.txt';
B_filt = dlmread(file_filt);
% get valid ball indices for trial
b_filt = get_valid_balls(B_filt,balls1,balls3,trial); %,3300);
%plot_trial(b1,b3,b_filt);

% initialize EKF
% initialize ball with high topspin (3000rpm or more)
w0 = [-50*2*pi;0;0]; %3000rpm topspin
spin.flag = false;
%[ekf,ball_func] = initFilterEKF(spin);
[ekf,ball_func] = init_const_spin_filter(w0);

figure;
% PREDICT WITH MATLAB FILTER VALUES FOR COMPARISON
%%{
net_y = dist_to_table - table_length/2;
balls3_pre_net = balls3{trial}(balls3{trial}(:,3) < net_y,:);
t_pre_net = balls3_pre_net(:,1) - balls3_pre_net(1,1);
t_end = balls1{trial}(end,1) - balls3_pre_net(1,1);
t_pred = t_pre_net(end):0.002:t_end; 
balls3_pre_net = balls3_pre_net(:,2:4);
ball_est_matlab = filterBallsEKF(t_pre_net,balls3_pre_net,ekf,spin);
ball_pred_net = predictTillLastBlob(ekf,t_pred,t_pred(1),...
                             ball_est_matlab(end,2:end));
% plot predicted balls
ball_preds_net = [ball_est_matlab(:,2:end);ball_pred_net'];
subplot(2,2,1);
plot_trial(b1,b3,ball_preds_net);

% % Plot predicted values around bounce
[~,idx_bounce] = min(balls3{trial}(:,4));
balls3_pre_bounce = balls3{trial}(1:idx_bounce,:);
t_pre_bounce = balls3_pre_bounce(:,1) - balls3_pre_bounce(1,1);
t_pred = t_pre_bounce(end):0.002:t_end;
balls3_pre_bounce = balls3_pre_bounce(:,2:4);
ball_est_matlab = filterBallsEKF(t_pre_bounce,balls3_pre_bounce,ekf,spin);
ball_pred_bounce = predictTillLastBlob(ekf,t_pred,t_pred(1),...
                                       ball_est_matlab(end,2:end));
ball_preds_bounce = [ball_est_matlab(:,2:end);ball_pred_bounce'];
subplot(2,2,2)
plot_trial(b1,b3,ball_preds_bounce);
%}

% PREDICT WITH CPP FILTER VALUES
%%{
% Plot predicted values starting from net
net_y = dist_to_table - table_length/2;
b_filt_pre_net = b_filt(b_filt(:,2) < net_y,:);
t_filt = 0.002 * (1:length(b_filt));
t_net = 0.002 * size(b_filt_pre_net,1);
ball_pred_net = predictTillLastBlob(ekf,t_filt,t_net,b_filt_pre_net(end,:));
% plot predicted balls
ball_preds_net = [b_filt_pre_net;ball_pred_net'];
subplot(2,2,3);
plot_trial(b1,b3,ball_preds_net);
% 
% % Plot predicted values around bounce
b_filt_pre_table_end = b_filt(b_filt(:,2) < dist_to_table,:);
[~,idx_bounce] = min(b_filt_pre_table_end(:,3));
b_filt_pre_bounce = b_filt(1:idx_bounce,:);
%idx_bounce = idx_bounce - 1;
t_bounce = 0.002 * size(b_filt_pre_bounce,1);
ball_pred_bounce = predictTillLastBlob(ekf,t_filt,t_bounce,b_filt_pre_bounce(end,:));
ball_preds_bounce = [b_filt_pre_bounce;ball_pred_bounce'];
subplot(2,2,4);
plot_trial(b1,b3,ball_preds_bounce);
%}