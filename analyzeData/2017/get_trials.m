% Estimate the number of trials within the data

function [balls1,balls3] = get_trials(b1,b3)

% if there is more than 1 second difference it means its a new trial
diffBall3 = diff(b3);
idxStart3 = find(diffBall3(:,1) > 1.0);
idxStart3 = [1;idxStart3+1];

tStart3 = b3(idxStart3,1);
num_trials = length(idxStart3);

balls3 = cell(1,num_trials);
balls1 = cell(1,num_trials);

for trial = 1:num_trials-1
    balls3{trial} = b3(idxStart3(trial):idxStart3(trial+1)-1,:);
    balls1{trial} = b1(b1(:,1) >= tStart3(trial) & b1(:,1) < tStart3(trial+1),:);
end
balls3{num_trials} = b3(idxStart3(num_trials):end,:);
balls1{num_trials} = b1(b1(:,1) >= tStart3(num_trials),:);