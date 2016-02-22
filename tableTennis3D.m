%% 3D table tennis : Barrett WAM practicing

clc; clear; close all; %dbstop if error;
draw = true;

%% Load table values
loadTennisTableValues;

%% Initialize Barrett WAM and the ball

initializeWAM;
% get joints, endeffector pos and orientation
[x,xd,o] = wam.calcRacketState([q0;qd0]);
[joint,ee,racket] = wam.drawPosture(q0);
q = q0; qd = qd0;
% initialize ball and and all interaction models
stdPos = 0.1;
stdVel = 0.1;
ball = Ball(stdPos,stdVel);

%% Initialize EKF
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
params.ALG = 'RK4'; %'Euler'

ballFlightFnc = @(x,u,dt) discreteBallFlightModel(x,dt,params);
% very small but nonzero value for numerical stability
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);
filter = EKF(dim,ballFlightFnc,mats);

ballPred(1:6,1) = [ball.pos;ball.vel];

%% Prepare the animation

if draw
figure;
%uisetcolor is useful here
orange = [0.9100 0.4100 0.1700];
gray = [0.5020    0.5020    0.5020];
numPoints = 100;
[ballMeshX,ballMeshY,ballMeshZ] = sphere(numPoints);
ballMeshX = ball_radius * ballMeshX;
ballMeshY = ball_radius * ballMeshY;
ballMeshZ = ball_radius * ballMeshZ;
ballSurfX = ball.pos(1) + ballMeshX;
ballSurfY = ball.pos(2) + ballMeshY;
ballSurfZ = ball.pos(3) + ballMeshZ;
% transform into base coord.
h1 = surf(ballSurfX,ballSurfY,ballSurfZ);
set(h1,'FaceColor',orange,'FaceAlpha',1,'EdgeAlpha',0);
hold on;
h10 = plot3(joint(:,1),joint(:,2),joint(:,3),'k','LineWidth',10);
endeff = [joint(end,:); ee];
h11 = plot3(endeff(:,1),endeff(:,2),endeff(:,3),'Color',gray,'LineWidth',5);
h12 = fill3(racket(1,:), racket(2,:), racket(3,:), 'r');

%h2 = scatter3(ballPred(1,1),ballPred(2,1),ballPred(3,1),20,'b','filled');
%h3 = scatter3(X0(1,1),X0(2,1),X0(3,1),20,'k','filled');
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
end

%% Main control and estimation loop

% initialize indices and time
WAIT = 0;
PREDICT = 1;
HIT = 2;
FINISH = 3; % only when practicing solo
robotIdx = 1; 
stage = 0; % WAIT
tSim = 0.0;
dt = 0.01;

% flags for the main loop
numBounce = 0;
numTrials = 0;
maxTime2Hit = 0.6;
maxWait = 3.0;
numLands = 0;

% initialize the filters state with sensible values
filter.initState([ball.pos;ball.vel + sqrt(eps)*randn(3,1)],eps);
ballNoisyPos = ball.pos;

while numTrials < 50
    
    %% Reset criterion
    if tSim > maxWait 
        numTrials = numTrials + 1;
        tSim = 0.0; 
        robotIdx = 1; 
        stage = 0;
        clc;
        % resetting the drawing
        if draw && numBounce == 1
            set(h2,'Visible','off');
            set(h3,'Visible','off');
        end
        % reset the ball
        if ball.isLANDED
            numLands = numLands+1;
        end
        ball = Ball(stdPos,stdVel);
        % reset the filter
        filter.initState([ball.pos;ball.vel],eps);
    end
    
    %% Evolve ball flight, table and racket contact interactions
    
    racketStr.pos = x(:,robotIdx);
    racketStr.vel = xd(:,robotIdx);
    racketRot = quat2Rot(o(:,robotIdx));
    racketStr.normal = racketRot(:,3);
    ball.evolve(dt,racketStr);
    tSim = tSim + dt;
    
    %% WAIT IF BALL IS INCOMING
    
    % Estimate the ball state
    filter.linearize(dt,0);
    filter.predict(dt,0);
    filter.update(ball.pos,0);     
    
    % If it is coming towards the robot consider moving
    velEst = filter.x(4:6);
    if velEst(2) > 0 && stage == WAIT
        stage = PREDICT;
    end    
    
    %% PREDICT BALL PATH
    if stage == PREDICT
        xSave = filter.x;
        PSave = filter.P;
        % if filter state goes below threshold bounce is highly likely
        tol = 1e-2;
        if filter.x(3) < table_z + tol && ...
           abs(filter.x(2) - (dist_to_table - table_length/2)) < table_length/2 && ...
           abs(filter.x(1)) < table_width/2           
            % if it bounces do not hit
            disp('Ball bounces on opp table! Not hitting!');
            stage = FINISH;
        end
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
    end
    
    %% GENERATE TRAJECTORY AND MOVE
    if stage == PREDICT && time2PassTable <= maxTime2Hit         
        
        % FOR DEBUGGING INIT TO ACTUAL POS AND VEL
        filter.initState([ball.pos;ball.vel],PSave);
        TABLE.DIST = dist_to_table;
        TABLE.LENGTH = table_length;
        TABLE.Z = table_z;
        TABLE.WIDTH = table_width;
        [ballPred,ballTime,numBounce,~] = predictBallPath(dt,filter,TABLE);
        filter.initState([ball.pos;ball.vel],PSave);
        if numBounce ~= 1
            disp('Ball is not valid! Not hitting!');
            stage = FINISH;
        else
            robotIdx = 0;
            stage = HIT;       
            % Calculate ball outgoing velocities attached to each ball pos
            tic
            fast = true;
            % land the ball on the centre of opponents court
            desBall(1) = 0.0;
            desBall(2) = dist_to_table - 3*table_y/2;
            desBall(3) = table_z + ball_radius;
            time2reach = 0.5; % time to reach desired point on opponents court
            racketDes = calcRacketStrategy(desBall,ballPred,ballTime,time2reach,fast);
            elapsedTimeForCalcDesRacket = toc;
            fprintf('Elapsed time for racket computation: %f sec.\n',...
                elapsedTimeForCalcDesRacket);

            % Compute traj here
            % define virtual hitting plane (VHP)
            VHP = -0.4;
            [q,qd,qdd] = generate3DTTTwithVHP(wam,VHP,ballPred,ballTime,q0);
            %[q,qd,qdd] = generateOptimalTTT(wam,racketDes,dt,q0);            
            [q,qd,qdd] = wam.checkJointLimits(q,qd,qdd);
            [x,xd,o] = wam.calcRacketState([q;qd]);

            if draw
            % Debugging the trajectory generation 
            h3 = scatter3(x(1,:),x(2,:),x(3,:));
            h2 = scatter3(ballPred(1,:),ballPred(2,:),ballPred(3,:));
            end
        end

    end % end predict        

    %% ANIMATE BOTH THE ROBOT AND THE BALL
    
    % Move the robot
    if stage == HIT 
        robotIdx = robotIdx+1;
    end
    
    % if movement finished revert to waiting
    if robotIdx > size(q,2)
        robotIdx = size(q,2);
        stage = FINISH;
    end
    
    [joint,ee,racket] = wam.drawPosture(q(:,robotIdx));
    endeff = [joint(end,:);ee];
    
    if draw
    
    % ball
    set(h1,'XData',ball.pos(1) + ballMeshX);
    set(h1,'YData',ball.pos(2) + ballMeshY);
    set(h1,'ZData',ball.pos(3) + ballMeshZ);
    
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
    pause(0.001);
    end
    
    
end

hold off;