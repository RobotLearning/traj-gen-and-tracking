% Predict with a flight model the velocity just before/after rebound
% before = 1, after = 0 for flagbefore

function vel = calcReboundVel(ballInit,ball,time,filter,flagbefore)
    
    eps = 1e-6;
    C = [eye(3),zeros(3)];
    filter.initState(ballInit(:),eps);
    u = zeros(1,size(ball,1));
    [xEKFSmooth, ~] = filter.smooth(time,ball',u);
    
    % Predict if necessary to get velocity before bounce
    dt = 0.001;
    
    if flagbefore % predict ahead 
        % make sure that state is before rebound
        ind = find(xEKFSmooth(6,:) < 0, 1, 'last');
        filter.initState(xEKFSmooth(:,ind),eps);
        while filter.x(6) < 0
            vel = filter.x(4:6);
            filter.predict(dt,0);
        end    
    else
        % make sure that state is ahead of rebound
        ind = find(xEKFSmooth(6,:) > 0, 1);
        stateFirst = xEKFSmooth(:,ind);
        filter.initState(stateFirst,eps);
        while filter.x(6) > 0
            vel = filter.x(4:6);
            filter.predict(-dt,0);
        end
    end
        
    % Return velocity before/after bounce
    vel = vel(:);
    
end