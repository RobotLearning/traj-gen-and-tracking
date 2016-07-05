% Predict with kNN

function [t_pred,ball_pred] = predictKNN(t_set,ball_set,t_star,ball_star,num_neigh)

% assume t_star is identical to the values in t_set 
t_star = t_star - t_star(1);
numRevealed = length(t_star);

% get the balls from the ball_set
numSet = length(ball_set);
balls = zeros(numRevealed*3,numSet);
for i = 1:numSet
    balls_at_set = ball_set{i}(1:numRevealed,:)';
    balls(:,i) = balls_at_set(:);
end

% find the closest matches and take average
b_star = ball_star';
b_star = b_star(:);
diff = repmat(b_star,1,numSet) - balls;
%[~,idx] = min(diag(diff'*diff));

[distances,idxs] = sortrows(diag(diff'*diff));
k = num_neigh;
[t_pred,ball_pred] = matchPredictions(idxs(1:k),t_set,ball_set);

% t_pred = t_set{idx};
% ball_pred = ball_set{idx};
end

% Function that tries to take average 
function [t_pred,ball_pred] = matchPredictions(idxs,t_set,ball_set)

    t_sync = t_set{idxs(1)};
    idx_sync = 1:length(t_sync);
    for i = 2:length(idxs)
        [~,t_sync,idx_sync] = synchronizeBalls(t_sync,t_set{idxs(i)},ball_set{idxs(i)},1e-2);
    end
    t_pred = t_sync;
    ball_pred = zeros(length(t_pred),3);
    for i = 1:length(idxs)
        ball_pred = ball_pred + ball_set{idxs(i)}(idx_sync,:);
    end
    ball_pred = ball_pred/length(idxs);

end