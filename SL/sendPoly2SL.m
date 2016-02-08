%% Sending polynomials to SL with the VHP method

% First get the ball and time and status 
% Start running the EKF and predict time2reach VHP
% Continue getting ball estimates until time2reach < maxTime

clc; clear all; close all;

%% Load table values

% load table parameters
loadTennisTableValues;

% land the ball on the centre of opponents court
desBall(1) = 0.0;
desBall(2) = dist_to_table - 3*table_y/2;
desBall(3) = table_z + ball_radius;
time2reach = 0.5; % time to reach desired point on opponents court

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

% initialize the filters state with sensible values
guessBallInitVel = [-0.80; 7.0; 2.0];
filter.initState([ball_cannon(:); guessBallInitVel],eps);

%% Initialize Barrett WAM

initializeWAM;

%% Create the socket
host = 'localhost'; 
port = '7646';
address = sprintf('tcp://%s:%s',host,port);
context = zmq.core.ctx_new();
socket  = zmq.core.socket(context, 'ZMQ_REQ');
zmq.core.connect(socket, address);

%% Get initial positioning of robot
msg = [uint8(3), typecast(uint32(1),'uint8'), uint8(0)];
data = typecast(msg, 'uint8'); 
zmq.core.send(socket, data);
response = zmq.core.recv(socket);

% get q,q0
STR = decodeResponseFromSL(response);
q0 = STR.robot.traj.q;
qd0 = STR.robot.traj.qd;
t0 = STR.robot.traj.time;

%% Trajectory generation

%cleanup = onCleanup(@() disconnectFromSL(socket,address,context));
ballPredicted = false;
ballPreTime = 0.0;
maxTime2Hit = 0.6;
maxPredictHorizon = 0.8;
time2PassTable = 1.0;
numTrials = 0;

while numTrials < 5
    
    % Get the ball positions
    msg = [uint8(4), typecast(uint32(1),'uint8'),uint8(0)];
    data = typecast(msg,'uint8');
    zmq.core.send(socket,data);
    response = zmq.core.recv(socket);
    STR = decodeResponseFromSL(response);
    ballPos = STR.ball.pos;
    ballTime = STR.ball.time;
    camId = STR.ball.cam.id;
    
    %% ESTIMATE BALL STATE
    if time2PassTable >= maxTime2Hit
        % predict balls current state
        dt = ballTime - ballPreTime;
        ballPreTime = ballTime;        
        filter.linearize(dt,0);
        filter.predict(dt,0);
        filter.update(ballPos,0);
        xSave = filter.x;
        PSave = filter.P;
        % update the time it takes to pass table
        yBallEst = filter.x(2);
        tPredIncrement = 0.02;
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
            dt = 0.02;
            predictLen = floor(predictHorizon / dt);
            ballPred = zeros(6,predictLen);
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
            velOut = zeros(3,size(ballPred,2));
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
            
            %% COMPUTE TRAJECTORY HERE AND SEND TO SL
                      
            [q,qd,qdd] = wam.generate3DTTTwithVHP(ballPred,ballTime,q0); 
            %[q,qd,qdd] = wam.generateOptimalTTT(racketDes,ballPred,ballTime,q0);
            timeSteps = size(q,2);
            ts = repmat(-1,1,timeSteps); % start immediately
            poly = [q;qd;qdd;ts];
            poly = poly(:);
            poly = typecast(poly,'uint8');
            % 1 is for clear
            % 2 is for push back
            N = typecast(uint32(timeSteps),'uint8');
            poly_zmq = [uint8(1), uint8(2), N, poly', uint8(0)];
            data = typecast(poly_zmq, 'uint8');
            zmq.core.send(socket, data);
            % response = zmq.core.recv(socket);

            numTrials = numTrials + 1;
        end % end predict        
        
    end
end
 
% disconnect from zmq and SL
disconnectFromSL(socket,address,context)