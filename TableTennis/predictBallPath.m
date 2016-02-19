function [ballPred,ballTime,numBounce] = predictBallPath(dt,filter)

loadTennisTableValues();
robotTableCenterY = dist_to_table - table_length/4;
% comingDown variable is used because z-derivative estimate might not
% be valid to predict number of bounces correctly
numBounce = 0;
comingDown = true;
tol = 2e-2;
maxPredictHorizon = 0.8;
predictHorizon = maxPredictHorizon;
predictLen = floor(predictHorizon / dt);
ballPred = zeros(6,predictLen);       
for j = 1:predictLen
    %filter.linearize(dt,0);
    filter.predict(dt,0);
    ballPred(:,j) = filter.x;
    
    % This part is for estimating num of bounces
    if filter.x(3) < table_z + tol && ...
       abs(filter.x(1)) < table_width/2 && ...
       abs(filter.x(2) - robotTableCenterY) < table_length/4 && comingDown
   
        numBounce = numBounce + 1;
        comingDown = false;
    else if filter.x(6) < 0 && ~comingDown
        comingDown = true;
        end
    end
end
% for now only considering the ball positions after table
ballTime = (1:predictLen) * dt;