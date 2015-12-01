%% Testing minimum principle for table tennis

clc; clear; close all;

% 1. Estimate initial pos and vel with kalman smoother up to net/bounce
% 2. Predict the ball trajectory after bounce
% 3. Fix the center of the table as the goal position
%    Solve BVP problems (bvp4c) by sampling from the predicted traj.
% 4. Using the contact model estimate racket vels. and feasible
% orientations and interpolate through whole traj.
%
% 5. Using an optimal control toolbox
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

% for the second demo only
%dist_to_table = dist_to_table - 0.25;

table_z = floor_level - table_height;
table_x = table_center + table_width/2;
table_y = table_length/2;

T1 = [table_center - table_x; 
    dist_to_table - table_length; 
    table_z];
T2 = [table_center + table_x;
    dist_to_table - table_length;
    table_z];
T3 = [table_center + table_x;
    dist_to_table;
    table_z];
T4 = [table_center - table_x;
    dist_to_table;
    table_z];
T = [T1,T2,T3,T4]';

net1 = [table_center - table_x;
        dist_to_table - table_y;
        table_z];
net2 = [table_center + table_x;
        dist_to_table - table_y;
        table_z];
net3 = [table_center + table_x;
        dist_to_table - table_y;
        table_z + net_height];
net4 = [table_center - table_x;
        dist_to_table - table_y;
        table_z + net_height];

net = [net1,net2,net3,net4]';

%% Initialize Barrett WAM

initializeWAM;
% initialize the ball gun
ball_cannon(1) = table_center + 0.4;
ball_cannon(2) = dist_to_table - table_length - 0.2; % since dist_to_table is negative
ball_cannon(3) = floor_level - table_height + 0.15;

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
% coeff of restitution-friction vector [TO BE LEARNED!]
params.CFTX = CFTX;
params.CFTY = CFTY;
params.CRT = CRT;

funState = @(x,u,dt) symplecticFlightModel(x,dt,params);
% very small but nonzero value for numerical stability
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);
filter = EKF(dim,funState,mats);

%% Main control and estimation loop

finished = false;
predict = false;
predict_horizon = 0.8;
ball(1:3,1) = ball_cannon;
ball(4:6,1) = [-0.9 4.000 3.2] + 0.05 * randn(1,3);
filter.initState([ball(1:3,1);ball(4:6,1) + sqrt(eps)*randn(3,1)],eps);

% initialize indices and time
i = 1;
dt = 0.01;
ballNoisyPos = ball(1:3,1);
timeSimTotal = 1.2;
t = 0.0;

while t < timeSimTotal && ~finished
    
    ball(:,i+1) = funState(ball(:,i),0,dt);
    
    %% Get noisy ball positions till ball passes the net
    if ballNoisyPos(2) <= dist_to_table - table_y - ball_radius;
        % add noise
        ballNoisyPos = ball(1:3,i+1) + sqrt(eps) * randn(3,1);
        filter.linearize(dt,0);
        filter.update(ballNoisyPos,0);
        filter.predict(dt,0);
        % robot just waits here
    else
        if ~predict
            %% predict ball trajectory
            predictHorizon = min(predict_horizon,timeSimTotal - t);
            predictLen = floor(predictHorizon / dt);
            ballPred = zeros(6,predictLen);
            for j = 1:predictLen
                filter.linearize(dt,0);
                filter.predict(dt,0);
                ballPred(:,j) = filter.x;
            end
            predict = true;
        end
        
        %% Calculate the intersection of ball path with robot workspace
        
        % for now only considering the ball positions after table
        tol = 5e-2;
        idxAfterTable = find(ballPred(2,:) > dist_to_table + tol);
        ballIncoming = ballPred(:,idxAfterTable);
        timeIncoming = idxAfterTable * dt;
        minTimeToHit = timeIncoming(1);
        ballTrjPred = [timeIncoming; ballIncoming];
        
        %% Calculate ball outgoing velocities attached to each ball pos
        
        % desired pos is in the centre of opponents court
        time2reach = 0.5;
        solinit = bvpinit(linspace(0,4,5),[-1 0]);
        sol = bvp4c(@func,@bc,solinit);
        
        %% OPTIMAL CONTROL HERE
        
        
    end
            
        
        
    t = t + dt;
    i = i + 1;
end

%% Plot the balls path

figure(4);
scatter3(ball(1,:),ball(2,:),ball(3,:),'r');
hold on;
%scatter3(r1,r2,r3,'b');
scatter3(ballPred(1,:),ballPred(2,:),ballPred(3,:),'b');
%scatter3(ballEKFSmooth(1,:),ballEKFSmooth(2,:),ballEKFSmooth(3,:));
title('Ball-robot interaction');
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
legend('ball', 'ballPred');
%legend('ball','robot');
fill3(T(:,1),T(:,2),T(:,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
hold off;
