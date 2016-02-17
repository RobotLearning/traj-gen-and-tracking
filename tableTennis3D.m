%% 3D table tennis with 2 Barrett WAMs

clc; clear; close all; dbstop if error;

% For a match run
% 
% initializeWAM;
% % change base and initialize wam2
% tt = TableTennis(wam,wam2);
% tt.match();


% VHP method vs. optimal control approach

%% Load table values

% load table parameters
loadTennisTableValues;

% land the ball on the centre of opponents court
desBall(1) = 0.0;
desBall(2) = dist_to_table - 3*table_y/2;
desBall(3) = table_z + ball_radius;
time2reach = 0.5; % time to reach desired point on opponents court

%% Initialize Barrett WAM and the ball

initializeWAM;
ball = Ball();

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

% get joints, endeffector pos and orientation
[x,xd,o] = wam.calcRacketState([q0;qd0]);
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

%% Main control and estimation loop

% flags for the main loop
numTrials = 0;
ballPredicted = false;
% maximum horizon to predict
time2PassTable = 1.0;
maxTime2Hit = 0.6;
maxPredictHorizon = 0.8;
maxWait = 4.0;

% initialize the filters state with sensible values
filter.initState([ball.pos;ball.vel + sqrt(eps)*randn(3,1)],eps);

% initialize indices and time
robotIdx = 1; 
tSim = 0.0;
dt = 0.01;
ballNoisyPos = ball.pos;

while numTrials < 50
    
    %% Reset criterion
    if tSim > maxWait 
        numTrials = numTrials + 1;
        tSim = 0.0; 
        robotIdx = 1; 
        ballPredicted = false;
        time2PassTable = 1.0;
        % resetting the drawing
        set(h2,'Visible','off');
        set(h3,'Visible','off');
        [joint,ee,racket] = wam.drawPosture(q0);
        endeff = [joint(end,:); ee];
        % reset the ball
        ball = Ball();
        % reset the filter
        filter.initState([ball.pos;ball.vel + sqrt(eps)*randn(3,1)],eps);
    end
    
    %% Evolve ball flight, table and racket contact interactions
    
    racketStr.pos = x(:,robotIdx);
    racketStr.vel = xd(:,robotIdx);
    racketRot = quat2Rot(o(:,robotIdx));
    racketStr.normal = racketRot(:,3);
    ball.evolve(dt,racketStr);
    tSim = tSim + dt;
    
    
    %% ESTIMATE BALL STATE
    if time2PassTable >= maxTime2Hit
        % add noise
        ballNoisyPos = ball.pos + sqrt(eps) * randn(3,1);
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
        if ~ballPredicted            
            predictHorizon = maxPredictHorizon;
            predictLen = floor(predictHorizon / dt);
            ballPred = zeros(6,predictLen);
            % FOR DEBUGGING
            filter.initState([ball.pos;ball.vel],PSave);
            for j = 1:predictLen
                %filter.linearize(dt,0);
                filter.predict(dt,0);
                ballPred(:,j) = filter.x;
            end
            ballPredicted = true;        

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
                      
            % define virtual hitting plane (VHP)
            %VHP = -0.6;
            %[q,qd,qdd] = wam.generate3DTTTwithVHP(VHP,ballPred,ballTime,q0);
            [q,qd,qdd] = wam.generateOptimalTTT(racketDes,ballPred,ballTime,q0);            
            [q,qd,qdd] = wam.checkJointLimits(q,qd,qdd);
            [x,xd,o] = wam.calcRacketState([q;qd]);
            
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

hold off;