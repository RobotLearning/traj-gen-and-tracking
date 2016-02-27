% Prefilter balls 
% Remove balls that are the same location and/or time

function [ballCleanObs,ballCleanTime,lastBallPos,lastBallTime] = ...
                  prefilter(ballObs,ballTime,lastBallPos,lastBallTime,table)
    
    % this is to throw away really bad observations as to not confuse the
    % filter
    zMax = 0.5;
    zMin = table.Z;
    assert(size(ballObs,2) == length(ballTime), 'time and obs num doesnt match!');
    tolTime = 1e-3;
    tolPos = 1e-2;
    ballCleanObs = [];
    ballCleanTime = [];
    for i = 1:length(ballTime)
        if (ballTime(i) - lastBallTime) > tolTime && ...
           abs(ballObs(2,i) - lastBallPos(2)) > tolPos && ...
            ballObs(3,i) < zMax && ballObs(3,i) > zMin
            % update the stack
            lastBallPos = ballObs(1:3,i);
            lastBallTime = ballTime(i);
            ballCleanObs = [ballCleanObs,ballObs(:,i)];
            ballCleanTime = [ballCleanTime,ballTime(i)];
        end
    end

%     % check for different ball pos
%     if size(balls,2) >= 2 && balls(4,end) ~= balls(4,end-1)
%         ballObsDer = (balls(1:3,end) - balls(1:3,end-1)) ./ ...
%                        (balls(4,end) - balls(4,end-1));
%         tol = 1e-2;
%         eps = 0.001;
%         % if ball appears to be coming towards the robot
%         if ballObsDer(2) > tol
%             filter.initState([balls(1:3,1); ballObsDer(:)],eps);
%         end
%     end