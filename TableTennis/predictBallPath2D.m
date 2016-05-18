% Predict the ball path for a fixed seconds into the future
% maxPredictHorizon indicates the fixed time
% includes also bounce prediction (num of estimated bounces)

function [ballPred,ballTime,numBounce,time2PassTable] = ...
            predictBallPath2D(dt,predictHorizon,filter,table)

dist_to_table = table.DIST;
table_length = table.LENGTH;
table_z = table.Z;
robotTableCenterY = dist_to_table - table_length/4;

% init necessary variables
tol = 2e-2;
predictLen = floor(predictHorizon / dt);
ballPred = zeros(4,predictLen);       
ballTime = (1:predictLen) * dt;
time2PassTable = Inf;

% save filter state before prediction
xSave = filter.x;
PSave = filter.P;

% comingDown variable is used because z-derivative estimate might not
% be valid to predict number of bounces correctly
numBounce = 0;
comingDown = true;

for j = 1:predictLen
    %filter.linearize(dt,0);
    filter.predict(dt,0);
    ballPred(:,j) = filter.x;
    
    % This part is for estimating num of bounces
    if filter.x(2) < table_z + tol && ...
       abs(filter.x(1) - robotTableCenterY) < table_length/4 && comingDown
   
        numBounce = numBounce + 1;
        comingDown = false;
    else if filter.x(4) < 0 && ~comingDown
        comingDown = true;
        end
    end
    
    % update the time it takes to pass table
    if filter.x(1) > dist_to_table && time2PassTable == Inf
        time2PassTable = j * dt;
    end
    
end

% revert back to saved state
filter.initState(xSave,PSave);