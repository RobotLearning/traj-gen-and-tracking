% Remove outliers in camera 1 with RANSAC method

function balls = remove_outliers(b1)

[~,idxClosest2Robot] = max(b1(:,3));
t1 = b1(1:idxClosest2Robot,1);
b1 = b1(1:idxClosest2Robot,2:4);

% use ransac to further prune outliers
% outlier detection again
outlierIdx = detectOutlierBalls(t1,b1,1);
inlierIdx = setdiff(1:length(t1),outlierIdx);
b1 = b1(inlierIdx,:);
t1 = t1(inlierIdx,:);
balls = [t1,b1];

end