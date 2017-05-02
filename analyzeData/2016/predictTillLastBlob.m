
function ballPred = predictTillLastBlob(filter,t,tstart,ball_lookup)

t = t(t > tstart) - tstart;
initVar = 1;
filter.initState(ball_lookup(:),initVar);
if isempty(filter.A)
    filter.linearize(t(1),0);
end
ballPred = predictBallPath(t,filter);
        
% Predict the ball path for a fixed seconds into the future
% maxPredictHorizon indicates the fixed time
% includes also bounce prediction (num of estimated bounces)

function ballPred = predictBallPath(t,filter)

ballPred = zeros(length(filter.x),length(t));       

% save filter state before prediction
xSave = filter.x;
PSave = filter.P;

t_last = 0.0;
for j = 1:length(t)
    diff_t = t(j) - t_last;
    %filter.linearize(dt,0);
    filter.predict(diff_t,0);
    ballPred(:,j) = filter.x;   
    t_last = t(j);
end

% revert back to saved state
filter.initState(xSave,PSave);        