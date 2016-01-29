%% Test 3D table tennis with Barrett WAM

clc; clear; close all; %dbstop if error;

% 1. Estimate pos and vel of ball with Kalman filter up to net/bounce
% 2. Predict the ball trajectory after net/bounce
% 3. Fix the center of the table as the desired goal position
%    Solve BVP problems (bvp4c) by sampling from the predicted traj.
% 4. Using the contact model estimate racket vels. and feasible
%    orientations and interpolate through whole traj.
%
% 5. Using minimum principle
%    State the initial position of the robot
%    State the input penalties
%    State the transversality conditions
%
% Transversality conditions:
% - positions of the ball should touch the racket (tubular region around
%   the path)
% - velocities of the racket are fixed
% - orientation of the racket could be fixed (although there is some
% freedom)
% - only using kinematics and jacobian
%
% 6. Compare with existing virtual hitting plane method
% - construct VHP
% - estimate pos and vel. as well as desired vel.
% - use inverse kinematics and polynomials to construct traj.
%
% 7. Evolve both with the nominal model (inverse dynamics) and compare
%   probability of returning the ball given actual ball velocities 
%   and/or control errors

%% Load table values

% load table parameters
loadTennisTableValues;

%% Initialize Barrett WAM

initializeWAM;

% initialize ball on the ball cannon with a sensible vel
ball(1:3,1) = ball_cannon;
ball(4:6,1) = [-0.9 4.000 3.2];% + 0.05 * randn(1,3);
ballPred(1:6,1) = ball(1:6,1); % to initialize drawing of predicted balls

% land the ball on the centre of opponents court
desBall(1) = 0.0;
desBall(2) = dist_to_table - 3*table_y/2;
desBall(3) = table_z + ball_radius;
time2reach = 0.5; % time to reach desired point on opponents court

%% initialize EKF
dim = 6;
eps = 1e-6; %1e-3;
C = [eye(3),zeros(3)];

params.C = Cdrag;
params.g = gravity;
params.zTable = table_z;
params.yNet = dist_to_table - table_y;
params.table_length = table_length;
params.table_width = table_width;
% coeff of restitution-friction vector
params.CFTX = CFTX;
params.CFTY = CFTY;
params.CRT = CRT;

ballFlightFnc = @(x,u,dt) symplecticFlightModel(x,dt,params);
% very small but nonzero value for numerical stability
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);
filter = EKF(dim,ballFlightFnc,mats);

%% Prepare the animation

% get joints, endeffector pos and orientation
[joint,ee,racket] = wam.drawPosture(q0);
q = q0; qd = qd0;

figure;
%uisetcolor is useful here
orange = [0.9100 0.4100 0.1700];
gray = [0.5020    0.5020    0.5020];
numPoints = 100;
[ballMeshX,ballMeshY,ballMeshZ] = sphere(numPoints);
ballMeshX = ball_radius * ballMeshX;
ballMeshY = ball_radius * ballMeshY;
ballMeshZ = ball_radius * ballMeshZ;
ballSurfX = ball(1,1) + ballMeshX;
ballSurfY = ball(2,1) + ballMeshY;
ballSurfZ = ball(3,1) + ballMeshZ;
% transform into base coord.
h1 = surf(ballSurfX,ballSurfY,ballSurfZ);
set(h1,'FaceColor',orange,'FaceAlpha',1,'EdgeAlpha',0);
hold on;
h10 = plot3(joint(:,1),joint(:,2),joint(:,3),'k','LineWidth',10);
endeff = [joint(end,:); ee];
h11 = plot3(endeff(:,1),endeff(:,2),endeff(:,3),'Color',gray,'LineWidth',5);
h12 = fill3(racket(1,:), racket(2,:), racket(3,:), 'r');

h2 = scatter3(ballPred(1,1),ballPred(2,1),ballPred(3,1),20,'b','filled');
h3 = scatter3(X0(1,1),X0(2,1),X0(3,1),20,'k','filled');
title('Ball-robot interaction');
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
tol_x = 0.1; tol_y = 0.1; tol_z = 0.3;
xlim([-table_x - tol_x, table_x + tol_x]);
ylim([dist_to_table - table_length - tol_y, tol_y]);
zlim([table_z - tol_z, table_z + 2*tol_z]);
legend('ball','robot');
fill3(T(:,1),T(:,2),T(:,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
net_width = 0.01;
%scatter3(ball(1,:),ball(2,:),ball(3,:),'r');
%scatter3(X(1,:),X(2,:),X(3,:),'k');
%scatter3(ballPred(1,:),ballPred(2,:),ballPred(3,:),'b');
%scatter3(ballEKFSmooth(1,:),ballEKFSmooth(2,:),ballEKFSmooth(3,:));

%% Main control and estimation loop

% flags for the main loop
contact = false;
numTrials = 0;
numHits = 0;
numLands = 0;
predict = false;
hit = false;
landOnTable = false;
touchNet = false;
outside = false;
cnt = 0.0; % counter
% maximum horizon to predict
time2PassTable = 1.0;
maxTime2Hit = 0.6;
maxPredictHorizon = 0.8;
maxWaitAfterHit = 1.0;

% initialize the filters state with sensible values
filter.initState([ball(1:3,1);ball(4:6,1) + sqrt(eps)*randn(3,1)],eps);

% initialize indices and time
ballIdx = 1; robotIdx = 1; t = 0.0;
dt = 0.01;
ballNoisyPos = ball(1:3,1);

while numTrials < 50
    
    %% Evolve ball flight, table and racket contact interactions
    
    % Reset criterion
    if ball(3,ballIdx) <= floor_level || touchNet || outside || cnt > maxWaitAfterHit %|| landOnTable
        numTrials = numTrials + 1;
        fprintf('Lands: %d/%d. Resetting...\n',numLands,numTrials);
        t = 0.0; ballIdx = 1; robotIdx = 1; cnt = 0.0;
        predict = false;
        hit = false;
        landOnTable = false;
        touchNet = false;  
        outside = false;
        time2PassTable = 1.0;
        % resetting the drawing
        set(h2,'Visible','off');
        set(h3,'Visible','off');
        [joint,ee,racket] = wam.drawPosture(q0);
        endeff = [joint(end,:); ee];
        % resetting the ball generation
        ball(1:3,1) = ball_cannon;
        ball(4:6,1) = [-0.9 4.000 3.2];% + 0.05 * randn(1,3);
        % reset the filter
        filter.initState([ball(1:3,1);ball(4:6,1) + sqrt(eps)*randn(3,1)],eps);
    end
    
    
    % check contact with racket    
    tol = ball_radius + 1e-2;
    l = min(robotIdx,size(q,2));
    % get racket pos, vel and orientation (on plane)
    [xRacket,velRacket,quatRacket] = wam.calcRacketState([q(:,l);qd(:,l)]);
    racketRot = quat2Rot(quatRacket);
    racketNormal = racketRot(1:3,3); % orientation along racket
    vecFromRacketToBall = ball(1:3,ballIdx) - xRacket;
    %racketPlane = racket_dir(:,l);
    projNormal = racketNormal*racketNormal'/(racketNormal'*racketNormal);
    projRacket = eye(3) - projNormal;
    distToRacketPlane = racketNormal'*vecFromRacketToBall;
    distOnRacketPlane = sqrt(vecFromRacketToBall'*vecFromRacketToBall - distToRacketPlane^2);
    if distToRacketPlane < tol && distOnRacketPlane < racket_radius && ~hit            
        %disp('A hit! Well done!');
        hit = true;        
        % TODO: find a precise hitting time with bisection                
        fprintf('Hit at x = %f, y = %f, z = %f\n',...
            ball(1,ballIdx),ball(2,ballIdx),ball(3,ballIdx));
        numHits = numHits + 1;
        % Change ball velocity based on contact model
        % get ball velocity
        velIn = ball(4:6,ballIdx);
        speedInAlongNormal = racketNormal'*velIn;
        speedRacketAlongNormal = racketNormal'*velRacket;
        velInAlongNormal = speedInAlongNormal.*racketNormal;
        velRacketAlongNormal = speedRacketAlongNormal.*racketNormal;
        % this is kept the same in mirror law
        velInAlongRacket = velIn - velInAlongNormal;
        velOutAlongNormal = velRacketAlongNormal + ...
            CRR * (velRacketAlongNormal - velInAlongNormal);
        velOut = velOutAlongNormal + velInAlongRacket;
        ball(4:6,ballIdx) = velOut;
    end
    
    % check contact with net and landing
    if hit
        cnt = cnt + dt; % run counter to terminate game
        tol = 2*net_width;
        % check contact with net
        if abs(ball(2,ballIdx) - (dist_to_table - table_y)) <= tol 
            if abs(ball(1,ballIdx)) > table_width/2
                disp('Outside of the net! Resetting...');
                outside = true;
            else if ball(3,ballIdx) < (table_z + net_height) 
                disp('Hit the net! Resetting...');
                touchNet = true;
                end
            end
        end
        % check for landing
        tol = 1e-2;
        if ~landOnTable && ball(3,ballIdx) <= table_z + ball_radius + tol && ball(6,ballIdx) < 0             
            landOnTable = true;
            %fprintf('Land on table at x = %f, y = %f\n', ball(1,ballIdx), ball(2,ballIdx));
            if abs(ball(1,ballIdx)) < table_width/2 && ...
               abs(ball(2,ballIdx) - (dist_to_table - 3*table_y/2)) < table_y/2
                %disp('Ball landed! Amazing!'); 
                fprintf('Land on table at x = %f, y = %f\n', ball(1,ballIdx), ball(2,ballIdx));
                numLands = numLands + 1;                
                %else
                %disp('Ball is not inside the court. Lost a point!')
            end
        end
    end

    % evolve ball according to flight model
    ball(:,ballIdx+1) = ballFlightFnc(ball(:,ballIdx),0,dt);
    
    %% ESTIMATE BALL STATE
    if time2PassTable >= maxTime2Hit
        % add noise
        ballNoisyPos = ball(1:3,ballIdx+1) + sqrt(eps) * randn(3,1);
        filter.linearize(dt,0);
        filter.update(ballNoisyPos,0);
        filter.predict(dt,0);
        xSave = filter.x;
        PSave = filter.P;
        % update the time it takes to pass table
        yBallEst = filter.x(2);
        tPredIncrement = dt;
        time2PassTable = 0.0;
        while yBallEst <= dist_to_table
             %filter.linearize(tPredIncrement,0);
             filter.predict(tPredIncrement,0);
             yBallEst = filter.x(2);
             time2PassTable = time2PassTable + tPredIncrement;
        end
        % revert back to saved state
        filter.initState(xSave,PSave);
    else
        %% PREDICT BALL TRAJECTORY
        if ~predict            
            predictHorizon = maxPredictHorizon;
            predictLen = floor(predictHorizon / dt);
            ballPred = zeros(6,predictLen);
            % FOR DEBUGGING
            filter.initState(ball(:,ballIdx),PSave);
            for j = 1:predictLen
                %filter.linearize(dt,0);
                filter.predict(dt,0);
                ballPred(:,j) = filter.x;
            end
            predict = true;        

            % for now only considering the ball positions after table
            tol = 5e-2;
            idxAfterTable = find(ballPred(2,:) > dist_to_table + tol);
            %ballPred(:,idxAfterTable);
            ballTime = (1:predictLen) * dt; %idxAfterTable * dt;
            minTimeToHit = ballTime(1);

            % Calculate ball outgoing velocities attached to each ball pos
            %%{
            tic
            fast = true; % compute outgoing vel with linear model for speed
            for j = 1:size(ballPred,2)
                
                velOut(:,j) = calcBallVelOut3D(desBall,ballPred(1:3,j),time2reach,fast);              
                % Use the inverse contact model to compute racket vels and normal
                % at every point                
                [rp,rv,ro] = calcDesRacketState(ballPred(1:3,j),velOut(:,j),ballPred(4:6,j));
                racketDes.time(j) = ballTime(j);
                racketDes.pos(:,j) = rp;
                racketDes.normal(:,j) = ro;
                racketDes.vel(:,j) = rv;
                
            end
            elapsedTimeForCalcDesRacket = toc;
            fprintf('Elapsed time for racket computation: %f sec.\n',...
                elapsedTimeForCalcDesRacket);
            %}
            
            %% COMPUTE TRAJECTORY HERE
                      
            [q,qd,qdd] = wam.generate3DTTTwithVHP(ballPred,ballTime,q0); 
            %[q,qd,qdd] = wam.generateOptimalTTT(racketDes,ballPred,ballTime,q0);
            [x,xd,o] = wam.calcRacketState([q;qd]);
            
            %{
            % first computing in cartesian space assuming a cartesian robot
            solve_method = 'BVP';            
            [tx,X,u,J] = mp(X0,timeIncoming,ballIncoming,racketVel,solve_method);
            % interpolate to get equal dt increments            
            teq = dt:dt:tx(end);
            Xeq = interp1(tx,X',teq)';
            %}
            
            % Debugging the trajectory generation 
            h3 = scatter3(x(1,:),x(2,:),x(3,:)); 
            h2 = scatter3(ballPred(1,:),ballPred(2,:),ballPred(3,:));
            
        end % end predict        

        %% Initiate the robot
        [joint,ee,racket] = wam.drawPosture(q(:,robotIdx));
        endeff = [joint(end,:);ee];
        robotIdx = min(robotIdx+1,size(q,2));
        
    end 
    
    %% ANIMATE BOTH THE ROBOT AND THE BALL
    
    % ball
    set(h1,'XData',ball(1,ballIdx) + ballMeshX);
    set(h1,'YData',ball(2,ballIdx) + ballMeshY);
    set(h1,'ZData',ball(3,ballIdx) + ballMeshZ);
    
    % robot joints
    set(h10,'XData',joint(:,1));
    set(h10,'YData',joint(:,2));
    set(h10,'ZData',joint(:,3));
    % robot endeffector
    set(h11,'XData',endeff(:,1));
    set(h11,'YData',endeff(:,2));
    set(h11,'ZData',endeff(:,3));
    % robot racket
    set(h12,'XData',racket(1,:));
    set(h12,'YData',racket(2,:));
    set(h12,'ZData',racket(3,:));
    
    drawnow;
    pause(0.002);
        
    t = t + dt;
    ballIdx = ballIdx + 1;
    
    
end

hold off;