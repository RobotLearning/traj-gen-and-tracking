%% Get ball positions from SL via matlab-zmq

obj = onCleanup(@() disconnectFromSL(socket,address,context));
clc; clear all; close all;

%% Create the socket

% wam or localhost
host = 'localhost';
port = '7646';
address = sprintf('tcp://%s:%s',host,port);
context = zmq.core.ctx_new();
socket  = zmq.core.socket(context, 'ZMQ_REQ');
zmq.core.connect(socket, address);

%% Load table values

% load table parameters
loadTennisTableValues;
title('Ball observations');
hold on;
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
tol_x = 0.1; tol_y = 0.1; tol_z = 0.3;
xlim([-table_x - tol_x, table_x + tol_x]);
ylim([dist_to_table - table_length - tol_y, tol_y]);
zlim([table_z - tol_z, table_z + 3*tol_z]);
fill3(T(1:4,1),T(1:4,2),T(1:4,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
net_width = 0.01;

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
guessBallInitVel = [-1.08; 4.80; 3.84]; %1.0 * [-0.80; 7.0; 2.0];
filter.initState([ball_cannon(:); guessBallInitVel],eps);

%% Clear the ball positions
bufferLength = 1e6; %bytes
msg = [uint8(5), uint8(0)];
data = typecast(msg,'uint8');
zmq.core.send(socket,data);
response = zmq.core.recv(socket,bufferLength);

%% GET BALL POSITIONS

t = [];
ballObs = [];
ballTime = [];
ballRaw = [];
ballFilt = [];
lastBallPos = zeros(3,1);
lastBallTime = 0.0;
j = 1;
reset = false;
firsttime = true;
predicted = false;
maxBallSize = 100;
minBall2Predict = 5;
predictTime = 1.0;

table.DIST = dist_to_table;
table.LENGTH = table_length;
table.Z = table_z;
table.WIDTH = table_width;

while ~reset && (size(ballRaw,2) < maxBallSize)

    msg = [uint8(4), typecast(uint32(1),'uint8'),uint8(0)];
    data = typecast(msg,'uint8');
    zmq.core.send(socket,data);
    response = zmq.core.recv(socket,bufferLength);
    STR = decodeResponseFromSL(response);
    ballObs = STR.ball.pos;
    ballTime = STR.ball.time;
    ballCam = STR.ball.cam.id;
    
    if ~isempty(ballTime)
        % get number of observations
        disp('Processing data...');
        [ballObs,ballTime,lastBallPos,lastBallTime] = ...
            prefilter(ballObs,ballTime,lastBallPos,lastBallTime,table);
        numObs = length(ballTime);
        
        if firsttime
            curTime = ballTime(1);
            %filter.initState([ballObs(:,1); guessBallInitVel],eps);           
            ballFilt(:,1) = ballObs(:,1);
            ballRaw(:,1) = ballObs(:,1);
            cam(:,1) = ballCam(:,1);
            firsttime = false;
            numObs = numObs - 1;
        end
         
        % keep the observed balls
        for i = 1:numObs
            ballRaw(:,j+i) = ballObs(:,i);
            cam(:,j+i) = ballCam(:,i);
            t(j+i) = ballTime(i);
        end
        
        % filter up to a point
        for i = 1:numObs
            dt = ballTime(i) - curTime;
            filter.linearize(dt,0);
            filter.predict(dt,0);
            filter.update(ballObs(:,i),0);
            curTime = ballTime(i);
            ballFilt(:,j+i) = filter.x(1:3);
        end 
        
        if size(ballRaw,2) > minBall2Predict && ...
                filter.x(2) > dist_to_table - table_length/2    
        % otherwise predict
            if ~predicted
                dtPred = 0.01;
                [ballPred,~,numBounce,time2PassTable] = ...
                    predictBall(dtPred,predictTime,filter,table);
                predicted = true;
            end
        end

        j = j + numObs;
    end
    
    % if suddenly there's a jump backwards stop
    tol = 1.0;
    if (size(ballRaw,2) > 2) && (abs(ballRaw(2,end) - ballRaw(2,end-1)) > tol)
        reset = true;
        % kick the last ball away
        ballRaw = ballRaw(:,1:end-1);
    end        


end

scatter3(ballRaw(1,:),ballRaw(2,:),ballRaw(3,:),'r');
scatter3(ballFilt(1,:),ballFilt(2,:),ballFilt(3,:),'y');
scatter3(ballPred(1,:),ballPred(2,:),ballPred(3,:),'b');
hold off;

% disconnect from zmq and SL
disconnectFromSL(socket,address,context);