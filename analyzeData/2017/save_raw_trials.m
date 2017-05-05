%% Save raw trials from one big file into single trials in dropbox folder

clc; clear; close all;
loadTennisTableValues;
file = '~/Dropbox/data/balls_0317.txt';
M = dlmread(file);
% ball data
B = M(:,1:10); % cam id and status cols

[B1,B3,~] = prune_ball_data(B);

% adapted from get_trials function
% if there is more than threshold second difference it means its a new trial
t_threshold = 0.8;
diffBall3 = diff(B3);
idxStart3 = find(diffBall3(:,1) > t_threshold);
idxStart3 = [1;idxStart3+1];
tStart3 = B3(idxStart3,1);

idx_start = floor(500 * tStart3);
num_trials = length(idx_start);

balls = cell(1,num_trials);

for trial = 1:num_trials-1
    balls{trial} = B(idx_start(trial):idx_start(trial+1)-1,:);
end
balls{num_trials} = B(idx_start(num_trials):end,:);

for trial = 1:num_trials
    % save trial
    filename = sprintf('~/Dropbox/data/real_ball_data_0317/balls_%d.txt',trial);
    dlmwrite(filename,balls{trial},'delimiter','\t','precision',4)
end

% % plot one to see the result
trial = 1;
% remove outliers from camera 1
[B1,B3,~] = prune_ball_data(balls{trial});
[balls1,balls3] = get_trials(B1,B3); % cell structure
balls1{trial} = remove_outliers(balls1{trial});
t1 = balls1{trial}(:,1);
b1 = balls1{trial}(:,2:4);
t3 = balls3{trial}(:,1);
b3 = balls3{trial}(:,2:4);
plot_trial(b1,b3);