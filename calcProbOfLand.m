%% Calculate probability of landing the ball on opponents court
% x - centre of racket trajectory 
% o - orientation of racket trajectory

function probLand = calcProbOfLand(obj,ballPopulation,q,qd)

verbose = false;

loadTennisTableValues();
% rotate everything by 90 degrees
R = [0 1; -1 0];
dt = 0.01;
tol = ball_radius;
M = size(ballPopulation,2);

% calculate probability of landing the ball  
numLands = 0;

% ball flight model
params.C = Cdrag;
params.g = gravity;
params.zTable = table_z;
params.yNet = dist_to_table - table_y;
params.table_length = table_length;
% coeff of restitution-friction vector
params.CFTY = CFTY;
params.CRT = CRT;
net_width_2d = 0.01;

funState = @(x,u,dt) symplecticFlightModel2D(x,dt,params);

% simulation-related variables
maxTime = 2.0;

for j = 1:M
    if verbose,fprintf('Simulating ball %d\n',j); end        
    % initialize ball state
    ball(:,1) = ballPopulation(:,j);
    k = 1;
    totalTime = 0.0;
    land = false;
    hit = false;                    
    while ~land && totalTime < maxTime  
        % compute hitting time
        if ~hit && k <= size(q,2) % before racket trj is completed            
            % get racket pos, vel and orientation (on plane)
            [xRacket,velRacket,mats] = obj.getEndEffectorState(q(:,k),qd(:,k));
            xRacket = R * xRacket;
            velRacket = R * velRacket;
            racketPlane = R * squeeze(mats(1:2,2,3,:)); % orientation along racket
            vecFromRacketToBall = ball(1:2,k) - xRacket;
            projPlane = racketPlane*racketPlane'/(racketPlane'*racketPlane);
            projOrth = eye(2) - projPlane;
            distToRacketPlane = norm(projOrth * vecFromRacketToBall);
            distOnRacketPlane = norm(projPlane * vecFromRacketToBall);
            if distToRacketPlane < tol && distOnRacketPlane < racket_radius
                hit = true;
                if verbose,fprintf('Hit at y = %f, z = %f\n',...
                        ball(1,k),ball(2,k));end
                % calculate outgoing velocity
                %timeHit = ballTime(k);
                % Change ball velocity based on contact model
                % get ball velocity
                velIn = ball(3:4,k);
                velInAlongNormal = projOrth * velIn;
                velRacketAlongNormal = projOrth * velRacket;
                % this is kept the same in mirror law
                velInAlongRacket = projPlane * velIn; 
                velOutAlongNormal = velRacketAlongNormal + ...
                    CRR * (velRacketAlongNormal - velInAlongNormal);
                velOut = velOutAlongNormal + velInAlongRacket;
                ball(3:4,k) = velOut;                                
            end
        end
        if hit % hit happened
            tol = 2*net_width_2d;
            % check contact with net
            if abs(ball(1,k) - (dist_to_table - table_y)) <= tol 
                if ball(2,k) < (table_z + net_height)
                    if verbose,disp('Hit the net! Resetting...');end
                    break;
                else
                    if verbose,disp('Passing over the net');end
                end
            end
            % check for landing
            tol = 0.01;
            if ball(2,k) <= table_z + ball_radius + tol && ball(4,k) < 0             
                land = true;
                if verbose,fprintf('Land at y = %f\n', ball(1,k));end
                if abs(ball(1,k) - (dist_to_table - 3*table_y/2)) < table_y/2
                    if verbose,disp('Ball landed! Amazing!');end
                    numLands = numLands + 1;                
                    %else
                    %disp('Ball is not inside the court. Lost a point!')
                end
            end
        end
        % evolve ball according to flight model
        ball(:,k+1) = funState(ball(:,k),0,dt);
        k = k + 1;    
        %fprintf('Idx = %d\n',k);
        totalTime = totalTime + dt;
    end    
end

probLand = numLands/M;