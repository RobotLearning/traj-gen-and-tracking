% Predict with kNN

function [tpred,ballpred] = predictKNN(t_set,ball_set,t_star,ball_star)

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

% find the closest match
% find the closest point among Xs
b_star = ball_star';
b_star = b_star(:);
diff = repmat(b_star,1,numSet) - balls;
[~,idx] = min(diag(diff'*diff));

tpred = t_set{idx};
ballpred = ball_set{idx};