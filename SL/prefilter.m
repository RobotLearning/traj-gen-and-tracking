% Prefilter balls 
% Remove balls that are the same location and/or time

function [filter,balls,numObs] = prefilter(filter,balls,ballObs,ballTime)
    
    numObs = 0;
    assert(size(ballObs,2) == length(ballTime), 'time and obs num doesnt match!');
    tolTime = 1e-3;
    tolPos = 1e-2;
    lastBallPos = balls(1:3,end);
    lastBallTime = balls(4,end);
    for i = 1:length(ballTime)
        if ballTime(i) - lastBallTime > tolTime && ...
           ballObs(2,i) - lastBallPos(2) > tolPos
            % update the stack
            lastBallPos = ballObs(1:3,i);
            lastBallTime = ballTime(i);
            lastBall = [lastBallPos;lastBall];
            balls = [balls,lastBall];
            numObs = numObs + 1;
        end
    end

    % check for different ball pos
    if size(balls,2) >= 2 && balls(4,end) ~= balls(4,end-1)
        ballObsDer = (balls(1:3,end) - balls(1:3,end-1)) ./ ...
                       (balls(4,end) - balls(4,end-1));
        tol = 1e-2;
        eps = 0.001;
        % if ball appears to be coming towards the robot
        if ballObsDer(2) > tol
            filter.initState([balls(1:3,1); ballObsDer(:)],eps);
        end
    end