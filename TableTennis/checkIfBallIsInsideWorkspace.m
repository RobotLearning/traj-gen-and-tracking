% Checks if the predicted ball trajectory happens to 
% intersect with the robot workspace
% Using linear programming with a bunch of workspace points

function intersect = checkIfBallIsInsideWorkspace(robot,ballPred)

    pts = robot.workspace;
    num = size(pts,2);
    % for speed choose a subset
    N = 1000;
    pts = pts(:,randi(num,1,N));

    f = zeros(N,1);
    lb = zeros(N,1);
    ub = ones(N,1);
    intersect = false;
    
    % first get the bounding box
    maxs = max(pts,[],2);
    mins = min(pts,[],2);
    xmax = maxs(1);
    xmin = mins(1);
    ymax = maxs(2);
    ymin = mins(2);
    zmax = maxs(3);
    zmin = mins(3);
    
    % no intersection if none of the trj pts are inside the bounding box
    if sum(ballPred(1,:) < xmax & ballPred(1,:) > xmin & ...
            ballPred(2,:) < ymax & ballPred(2,:) > ymin & ...
            ballPred(3,:) < zmax & ballPred(3,:) > zmin) ~= 0
        % check for each point
        M = size(ballPred,2);
        for i = 1:M
            if ballPred(1,i) < xmax && ballPred(1,i) > xmin && ...
                ballPred(2,i) < ymax && ballPred(2,i) > ymin && ...
                ballPred(3,i) < zmax && ballPred(3,i) > zmin
                % lets check workspace with linprog
                p = ballPred(1:3,i);   
                beq = [p(:);1];
                Aeq = [pts;ones(1,N)];
                options = optimoptions('linprog','Display','off');
                x0 = (1/N) * ones(N,1); 
                [~,~,exitflag,~] = linprog(f,[],[],Aeq,beq,lb,ub,x0,options);
    
                if exitflag ~= -2
                    % solution was found
                    intersect = true;
                    return;
                end
            end
        end
    end

end

