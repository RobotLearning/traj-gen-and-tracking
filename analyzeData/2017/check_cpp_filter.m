%% Checking the performance of EKF filter implemented in C++ 

clc; clear; close all;
loadTennisTableValues;
trial = 4;
num_init_balls = 12;
%file = ['~/Dropbox/data/real_ball_data_0317/balls_', int2str(trial), '.txt'];
file = ['~/Desktop/data/real_ball_data_050817/balls_',int2str(trial),'.txt'];
M = dlmread(file);
B = M(:,1:10); % ball data including cam id and status cols

[balls1full,balls3,~] = prune_ball_data(B);
%[balls1,balls3] = get_trials(B1,B3); % cell structure

% remove outliers from camera 1
balls1 = remove_outliers(balls1full);
t1 = balls1(:,1);
b1 = balls1(:,2:4);
t3 = balls3(:,1);
b3 = balls3(:,2:4);

% initialize EKF
% initialize ball with high topspin (3000rpm or more)
w0 = [-50*2*pi;0;0]; %3000rpm topspin
spin_state.flag = false; % estimation with spin doesnt work yet!
% spin_state.est = w0;
% spin_state.Clift = Clift;
% spin_state.known = false;
%[ekf,ball_func] = initFilterEKF(spin);
[ekf,ball_func] = init_const_spin_filter(w0);

%% GET CPP FILTERED BALL VALUES AND INSPECT OUTLIER REJECTION
%file_filt = ['~/Dropbox/data/real_ball_data_0317/balls_',int2str(trial),'_filtered.txt'];  
b_filt = M(:,11:end);
%B_filt = dlmread(file_filt);
% get valid ball indices for trial
%b_filt = get_valid_balls(B_filt,balls1,balls3);
% figure;
% subplot(2,1,1);
% plot_trial(b1,b3,b_filt);
% title('Balls WITHOUT outliers'); % outliers removed surgically in matlab
% subplot(2,1,2);
% b_filt_full = get_valid_balls(B_filt,balls1full,balls3);
% plot_trial(balls1full(:,2:4),b3,b_filt_full);
% title('Balls WITH outliers');

%% PREDICT WITH MATLAB FILTER VALUES FOR COMPARISON
%%{
net_y = dist_to_table - table_length/2;
balls3_pre_net = balls3(balls3(:,3) < net_y,:);
t_pre_net = balls3_pre_net(:,1) - balls3_pre_net(1,1);
t_end = balls1(end,1) - balls3_pre_net(1,1);
t_pred = t_pre_net(end):0.002:t_end; 
balls3_pre_net = balls3_pre_net(:,2:4);
ball_est_mat_net = filterBallsEKF(t_pre_net,balls3_pre_net,ekf,num_init_balls,spin_state);
ball_pred_net = predictTillLastBlob(ekf,t_pred,t_pred(1),...
                             ball_est_mat_net(end,2:end));
% plot predicted balls
figure;
subplot(2,2,1);
plot_trial(b1,b3,ball_pred_net');

% % Plot predicted values around bounce
[~,idx_bounce] = min(balls3(:,4));
balls3_pre_bounce = balls3(1:idx_bounce,:);
t_pre_bounce = balls3_pre_bounce(:,1) - balls3_pre_bounce(1,1);
t_pred = t_pre_bounce(end):0.002:t_end;
balls3_pre_bounce = balls3_pre_bounce(:,2:4);
ball_est_mat_bounce = filterBallsEKF(t_pre_bounce,balls3_pre_bounce,ekf,num_init_balls,spin_state);
ball_pred_bounce = predictTillLastBlob(ekf,t_pred,t_pred(1),...
                                       ball_est_mat_bounce(end,2:end));
subplot(2,2,2)
plot_trial(b1,b3,ball_pred_bounce');
%}

%% PREDICT WITH CPP FILTER VALUES
%%{
% Plot predicted values starting from net
b_filt_pre_net = b_filt(b_filt(:,2) < net_y,:);
t_filt = 0.002 * (1:length(b_filt));
t_net = 0.002 * size(b_filt_pre_net,1);
ball_pred_net = predictTillLastBlob(ekf,t_filt,t_net,b_filt_pre_net(end,:));
% plot predicted balls
subplot(2,2,3);
plot_trial(b1,b3,ball_pred_net');
% 
% % Plot predicted values around bounce
b_filt_pre_table_end = b_filt(b_filt(:,2) < dist_to_table,:);
[~,idx_bounce] = min(b_filt_pre_table_end(:,3));
b_filt_pre_bounce = b_filt_pre_table_end(1:idx_bounce-1,:);
%idx_bounce = idx_bounce - 1;
t_bounce = 0.002 * size(b_filt_pre_bounce,1);
ball_pred_bounce = predictTillLastBlob(ekf,t_filt,t_bounce,b_filt_pre_bounce(end,:));
subplot(2,2,4);
plot_trial(b1,b3,ball_pred_bounce');
%}