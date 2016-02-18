% Check if incoming ball will bounce once
% Note: ball must be incoming
function numBounce = estimateNumBounce(ballObs,ballPred,Z)

if ballPred(5,1) < 0
    warning('Ball may not be incoming!');
end

% look at ball observations and predictions and find num of bounces 
ball = [ballObs,ballPred(1:3,:)];
tol = 1e-2;
bounceIdx = ball(3,:) < (Z + tol) && ball(2,:) 
numBounces 


% we dont take serve into account
serve = false;
% if ball has already bounced on opponents court 
if ~serve && 

