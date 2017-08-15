%% Save raw trials from one big file into single trials in dropbox folder

clc; clear; close all;
loadTennisTableValues;
file = '~/Desktop/data/balls_050817.txt';
B = dlmread(file);

%[B1,B3,~] = prune_ball_data(B);
% if there is more than threshold y-difference it means its a new trial
thresh = -1.0;
diffBall3 = diff(B(:,9));
idx_start = [1; find(diffBall3 < thresh) + 1];
num_trials = length(idx_start);

balls = cell(1,num_trials);

for trial = 1:num_trials-1
    % add some vals from previous trial for testing cpp filter
    start_idx = max(idx_start(trial)-5,1);
    balls{trial} = B(start_idx:idx_start(trial+1)-1,:);
end
balls{num_trials} = B(idx_start(num_trials):end,:);

for trial = 1:num_trials
    % save trial
    filename = sprintf('~/Desktop/data/real_ball_data_050817/balls_%d.txt',trial);
    dlmwrite(filename,balls{trial},'delimiter','\t','precision',4)
end

% % plot one to see the result
trial = 1;
% remove outliers from camera 1
[b1,b3,~] = prune_ball_data(balls{trial});
[balls1,balls3] = get_trials(b1,b3); % cell structure
balls1{trial} = remove_outliers(balls1{trial});
t1 = balls1{trial}(:,1);
b1 = balls1{trial}(:,2:4);
t3 = balls3{trial}(:,1);
b3 = balls3{trial}(:,2:4);
plot_trial(b1,b3);