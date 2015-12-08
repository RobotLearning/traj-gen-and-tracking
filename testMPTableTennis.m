%% Testing minimum principle in table tennis

clc; clear; close all

% 1. Estimate initial pos and vel with kalman smoother up to net/bounce
% 2. Predict the ball trajectory after bounce
% 3. Fix the center of the table as the goal position
%    Solve BVP problems (bvp4c) by sampling from the predicted traj.
% 4. Using the contact model estimate racket vels. and feasible
% orientations and interpolate through whole traj.
%
% 5. Using minimum principle
% State the initial position of the robot
% State the input penalties
% State the transversality conditions
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
% initialize the arm with zero velocity on the right hand side
q0 = [1.8; -0.2; -0.1; 1.8; -1.57; 0.1; 0.3];
qd0 = zeros(N_DOFS,1);
Q0 = [q0; qd0];
[x0,xd0] = wam.kinematics(Q0);
X0 = [x0;xd0];
robot(:,1) = X0;
% initialize ball on the ball cannon with a sensible vel
ball(1:3,1) = ball_cannon;
ball(4:6,1) = [-0.9 4.000 3.2] + 0.05 * randn(1,3);
ballPred(1:6,1) = ball(1:6,1); % to initialize drawing of predicted balls
% land the ball on the centre of opponents court
desBall(1) = 0.0;
desBall(2) = dist_to_table - 3*table_y/2;
desBall(3) = table_z + ball_radius;
time2reach = 0.5; % time to reach desired point on opponents court

%% initialize EKF
dim = 6;
eps = 1e-3;
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

funState = @(x,u,dt) symplecticFlightModel(x,dt,params);
% very small but nonzero value for numerical stability
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);
filter = EKF(dim,funState,mats);

%% Prepare the animation

% get joints, endeffector pos and orientation
[j,e,r] = wam.drawPosture(q0);

figure;
%uisetcolor is useful here
orange = [0.9100 0.4100 0.1700];
gray = [0.5020    0.5020    0.5020];
h1 = scatter3(ball(1,1),ball(2,1),ball(3,1),20,orange,'filled');
hold on;

h10 = plot3(j(:,1),j(:,2),j(:,3),'k','LineWidth',10);
endeff = [j(end,:);e];
h11 = plot3(endeff(:,1),endeff(:,2),endeff(:,3),'Color',gray,'LineWidth',5);
h12 = fill3(r(1,:), r(2,:), r(3,:), 'r');

h2 = scatter3(ballPred(1,1),ballPred(2,1),ballPred(3,1),20,'b','filled');
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
%scatter3(ball(1,:),ball(2,:),ball(3,:),'r');
%scatter3(X(1,:),X(2,:),X(3,:),'k');
%scatter3(ballPred(1,:),ballPred(2,:),ballPred(3,:),'b');
%scatter3(ballEKFSmooth(1,:),ballEKFSmooth(2,:),ballEKFSmooth(3,:));

%% Main control and estimation loop

% flags for the main loop
finished = false;
contact = false;
collectData = 0;
predict = false;
% maximum horizon to predict
time2PassTable = 1.0;
maxTime2Hit = 0.6;
maxPredictHorizon = 0.8;

% initialize the filters state with sensible values
filter.initState([ball(1:3,1);ball(4:6,1) + sqrt(eps)*randn(3,1)],eps);

% initialize indices and time
i = 1; k = 1; t = 0.0;
dt = 0.01;
ballNoisyPos = ball(1:3,1);

while ~finished
    
    %% Evolve ball flight, table and racket contact interactions
    if ball(3,i) <= floor_level
        disp('Ball hit the floor. Resetting...');
        t = 0.0; i = 1; k = 1;
        predict = false;
        time2PassTable = 1.0;
        ball(1:3,1) = ball_cannon;
        ball(4:6,1) = [-0.9 4.000 3.2] + 0.05 * randn(1,3);
        filter.initState([ball(1:3,1);ball(4:6,1) + sqrt(eps)*randn(3,1)],eps);
    end
    
    ball(:,i+1) = funState(ball(:,i),0,dt);
    
    %TODO: check contact with racket
    if contact
        disp('A hit! Well done!');
        % change ball velocity based on contact model
    end
    
    %TODO: check simulation criterion
    if collectData == 50
        disp('Finishing up...');
        finished = true;
    end
    
    %% Get noisy ball positions till ball passes the net
    if time2PassTable >= maxTime2Hit
        % add noise
        ballNoisyPos = ball(1:3,i+1) + sqrt(eps) * randn(3,1);
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
        %% predict ball trajectory once
        if ~predict            
            predictHorizon = maxPredictHorizon;
            predictLen = floor(predictHorizon / dt);
            ballPred = zeros(6,predictLen);
            for j = 1:predictLen
                %filter.linearize(dt,0);
                filter.predict(dt,0);
                ballPred(:,j) = filter.x;
            end
            predict = true;        

            %% Calculate the intersection of ball path with robot workspace

            % for now only considering the ball positions after table
            tol = 5e-2;
            idxAfterTable = find(ballPred(2,:) > dist_to_table + tol);
            ballIncoming = ballPred; %ballPred(:,idxAfterTable);
            timeIncoming = (1:predictLen) * dt; %idxAfterTable * dt;
            minTimeToHit = timeIncoming(1);

            %% Calculate ball outgoing velocities attached to each ball pos

            for j = 1:size(ballIncoming,2)
                % desired pos is in the centre of opponents court
                velBallOut(1) = (desBall(1) - ballIncoming(1,j))/time2reach;
                velBallOut(2) = (desBall(2) - ballIncoming(2,j))/time2reach;
                velBallOut(3) = (desBall(3) - ballIncoming(3,j) - ...
                                0.5*gravity*time2reach^2)/time2reach;
                % initialize using a linear model (no drag)
                linFlightTraj = @(t) [ballIncoming(1:3,j) + velBallOut(:)*t;
                                      velBallOut(:)] + ...
                                      [0;0;0.5*gravity*t^2;0;0;gravity*t];
                flightModel = @(t,x) [x(4);
                                      x(5);
                                      x(6);
                              -Cdrag * x(4) * sqrt(x(4)^2 + x(5)^2 + x(6)^2);
                              -Cdrag * x(5) * sqrt(x(4)^2 + x(5)^2 + x(6)^2);
                      gravity - Cdrag * x(6) * sqrt(x(4)^2 + x(5)^2 + x(6)^2)];
                % boundary value condition
                bc = @(x0,xf) [x0(1) - ballIncoming(1,j);
                               x0(2) - ballIncoming(2,j);
                               x0(3) - ballIncoming(3,j);
                               xf(1) - desBall(1);
                               xf(2) - desBall(2);
                               xf(3) - desBall(3)];
                meshpoints = 50;
                solinit = bvpinit(linspace(0,time2reach,meshpoints),...
                                  linFlightTraj);
                sol = bvp4c(flightModel,bc,solinit);
                ballOut = deval(sol,0);
                velOut(:,j) = ballOut(4:6);

                %% Use the inverse contact model to compute racket vels and normal
                % at every point

                % FOR NOW USING ONLY THE MIRROR LAW
                normal(:,j) = (velOut(:,j) - ballIncoming(4:6,j)) ...
                              ./ norm(velOut(:,j) - ballIncoming(4:6,j),2);
                velOutAlongNormal = velOut(:,j)' * normal(:,j);
                velInAlongNormal = ballIncoming(4:6,j)' * normal(:,j);
                racketVelAlongNormal(j) = velOutAlongNormal + ... 
                    CRR * velInAlongNormal / (1 + CRR);
                % TO BE CHANGED 
                racketVel(:,j) = racketVelAlongNormal(j) * normal(:,j);
            end
            
            %% COMPUTE OPTIMAL TRAJECTORY HERE
            % first computing in cartesian space assuming a cartesian robot
            T0 = 0.5;
            B0 = ballIncoming(:,1);
            solve_method = 'BVP';
            
            [tx,X,u,J] = mp(X0,timeIncoming,ballIncoming,racketVel,solve_method);
            % interpolate to get equal dt increments
            teq = dt:dt:tx(end);
            Xeq = interp1(tx,X',teq)';
        end % end predict        

        %% Initiate the robot
        robot(:,k+1) = Xeq(:,min(k+1,size(Xeq,2)));
        k = k + 1;
        
    end 
    
    %% ANIMATE BOTH THE ROBOT AND THE BALL
    set(h1,'XData',ball(1,i));
    set(h1,'YData',ball(2,i));
    set(h1,'ZData',ball(3,i));
    set(h3,'XData',robot(1,k));
    set(h3,'YData',robot(2,k));
    set(h3,'ZData',robot(3,k));
    drawnow;
    pause(0.002);
        
    t = t + dt;
    i = i + 1;
    
    
end

hold off;