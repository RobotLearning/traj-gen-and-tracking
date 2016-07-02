% Estimate the number of trials within the data

function [b1,b3,ballEst,numTrials] = getTrialData(b1,b3,trial,dataset,filtData)

% if there is more than 1 second difference it means its a new trial
diffBall3 = diff(b3);
idxStart3 = find(diffBall3(:,1) > 1.0);

% if second dataset then dont keep first trial because it is bad
if dataset == 1
    idxStart3 = [1;idxStart3+1];
else
    idxStart3 = idxStart3+1;
end

tStart3 = b3(idxStart3,1);
numTrials = length(idxStart3);

% get the indices for plotting
b3 = b3(idxStart3(trial):idxStart3(trial+1)-1,:);
b1 = b1(b1(:,1) >= tStart3(trial) & b1(:,1) < tStart3(trial+1),:);

try
    ballEst = filtData(idxStart3(trial):idxStart3(trial+1)-1,:);
    % remove first 12 entries
    ballEst = ballEst(12+1:end,:);  
catch
    warning('Filter data from SL seems to be unavailable!');
    ballEst = [];
end