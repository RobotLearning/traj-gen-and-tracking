%% Sending polynomials to SL with the VHP method

% First get the ball and time and status 
% Start running the EKF and predict time2reach VHP
% Continue getting ball estimates until time2reach < maxTime

obj = onCleanup(@() disconnectFromSL(socket,address,context));
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
qInit = STR.robot.traj.q;
qdInit = STR.robot.traj.qd;
tInit = STR.robot.traj.time;

% Send to desired starting posture
Qinit = [qInit;qdInit];
dt = 0.002;
tf = 2.0;
p = generatePoly3rd(Qinit,Q0,dt,tf);
timeSteps = size(p,2);
ts = repmat(-1,1,timeSteps); % start immediately
poly = [p;ts];
poly = poly(:);
poly = typecast(poly,'uint8');
% 1 is for clear
% 2 is for push back
N = typecast(uint32(timeSteps),'uint8');
poly_zmq = [uint8(1), uint8(2), N, poly', uint8(0)];
data = typecast(poly_zmq, 'uint8');
zmq.core.send(socket, data);
response = zmq.core.recv(socket);
pause(tf);

%% Trajectory generation

waiting = true;
bufferLength = 1e6; %bytes
ballTime = [];
ballPreTime = [];
ballObs = [];
ballPreObs = [];
maxTime2Hit = 0.6;
maxPredictHorizon = 0.8;
time2PassTable = 1.0;
numTrials = 0;

% clear the ball positions
msg = [uint8(5), uint8(0)];
data = typecast(msg,'uint8');
zmq.core.send(socket,data);
response = zmq.core.recv(socket,bufferLength);

while numTrials <= 5
    
    %% GET BALL POSITIONS
    msg = [uint8(4), typecast(uint32(1),'uint8'),uint8(0)];
    data = typecast(msg,'uint8');
    zmq.core.send(socket,data);
    response = zmq.core.recv(socket,bufferLength);
    STR = decodeResponseFromSL(response);
    ballPreObs = ballObs;
    ballObs = STR.ball.pos;
    ballPreTime = ballTime;
    ballTime = STR.ball.time;
    
    %% WAITING STAGE
    
    if ~isempty(ballPreObs) && ~isempty(ballObs) && waiting
        ballObsDer = (ballObs(:,end) - ballPreObs(:,end)) ./ ...
                       (ballTime(end) - ballPreTime(end));
        %fprintf('Derivative of ball = %f\n', ballObsDer);
        tol = 1e-2;    
        % if ball appears on opponent side coming towards the robot
        if ballObs(2,end) > dist_to_table - table_y && ...
                  ballObsDer(2) > tol
            waiting = false;
            numTrials = numTrials + 1;
            time2PassTable = 1.0;
            % reset the filter
            filter.initState([ballObs(:,end);guessBallInitVel],eps);
        end    
    end

    %% ESTIMATE BALL STATE
    if time2PassTable >= maxTime2Hit && ~waiting
        % predict balls current state
        dt = ballTime(end) - ballPreTime(end);
        filter.linearize(dt,0);
        filter.predict(dt,0);
        filter.update(ballObs(:,end),0);
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
    end
    
    %% PREDICT BALL TRAJECTORY
    if ~waiting           
        predictHorizon = maxPredictHorizon;
        dt = 0.002;
        predictLen = floor(predictHorizon / dt);
        ballPred = zeros(6,predictLen);
        for j = 1:predictLen
            %filter.linearize(dt,0);
            filter.predict(dt,0);
            ballPred(:,j) = filter.x;
        end    
        ballTime = (1:predictLen) * dt;

        % Calculate ball outgoing velocities attached to each ball pos
        %{
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

        % define virtual hitting plane (VHP)
        VHP = -0.2;
        [q,qd,qdd] = wam.generate3DTTTwithVHP(VHP,ballPred,ballTime,q0); 
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
        response = zmq.core.recv(socket);

        tf = timeSteps * 0.002;
        %pause(tf);

        waiting = true;

    end % end predict        
        
    %}
end
 
% disconnect from zmq and SL
disconnectFromSL(socket,address,context);