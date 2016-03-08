%% Sending polynomials to SL with the VHP method

% First get the ball and time and status 
% Start running the EKF and predict time2reach VHP
% Continue getting ball estimates until time2reach < maxTime

obj = onCleanup(@() disconnectFromSL(socket,address,context));
clc; clear all; close all;

%% Load table values

% load table parameters
loadTennisTableValues;
tennisTable = T;

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

% wam or localhost
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
tf = 1.0;
p = generatePoly3rd(Qinit,Q0,dt,tf);
% change last velocity and acc command to zero
% p(8:end,end) = 0.0;
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
pause(2.0); % wait 5sec to decr to qd0 to 10^-7 range

%% Clear the ball positions
bufferLength = 1e6; %bytes
msg = [uint8(5), uint8(0)];
data = typecast(msg,'uint8');
zmq.core.send(socket,data);
response = zmq.core.recv(socket,bufferLength);

%% Load lookup table

% load the savefile
savefile = 'OfflineTrajGenLookupTable.mat';
load(savefile,'X','Y');

%% Trajectory generation

% initialize indices and time
WAIT = 0;
PREDICT = 1;
FINISH = 2;
stage = WAIT;
dt = 0.002;

table.DIST = dist_to_table;
table.LENGTH = table_length;
table.Z = table_z;
table.WIDTH = table_width;

% flags for the main loop
ballTime = [];
ballRaw = [];
numBounce = 0;
numTrials = 0;
minTime2Hit = 0.6;
lastBallPos = zeros(3,1);
lastBallTime = 0.0;
j = 1;
firsttime = true;
maxBallSize = 100;
minBall2Predict = 5;
predictTime = 0.6;
%numLands = 0; % not used for now

%* q0 should be set somewhere
%* how to do prefilter
%* correct finite state machine
%* check joint limits was commented - run table tennis again with lookup 100 times 

while numTrials < 1
     
    msg = [uint8(4), typecast(uint32(1),'uint8'),uint8(0)];
    data = typecast(msg,'uint8');
    zmq.core.send(socket,data);
    response = zmq.core.recv(socket,bufferLength);
    STR = decodeResponseFromSL(response);
    ballObs = STR.ball.pos;
    ballTime = STR.ball.time;
    ballCam = STR.ball.cam.id;
    
    [ballObs,ballTime,lastBallPos,lastBallTime] = ...
        prefilter(ballObs,ballTime,lastBallPos,lastBallTime,table);
    
    if ~isempty(ballTime)
        % get number of observations
        disp('Processing data...');
        numObs = length(ballTime);
        
        % if suddenly there's a jump backwards stop
        tol = 1.0;
        if size(ballRaw,2) > 2 && abs(ballRaw(2,end) - ballRaw(2,end-1)) > tol
            firsttime = true;
            stage = WAIT;
            ballRaw = []; ballFilt = []; cam = [];
            j = 1;
        end        

        if firsttime
            curTime = ballTime(1);
            filter.initState([ballObs(:,1); guessBallInitVel],eps);           
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

        j = j + numObs;
    end

    tol = 0.4;
    if size(ballRaw,2) > minBall2Predict && stage == WAIT && ...
            filter.x(2) > dist_to_table - table_length/2  && ...
            filter.x(2) < dist_to_table - table_length/2 + tol;
    % otherwise predict
        %dtPred = 0.01;
        %[ballPred,~,numBounce,time2PassTable] = ...
        %    predictBall(dtPred,predictTime,filter,table);
        %checkBounceOnOppTable(filter,table);
        stage = PREDICT;
    end


    % HIT THE BALL IF VALID
    if stage == PREDICT
        tic
        disp('Sending trj data');
        stage = FINISH;
        numTrials = numTrials + 1;
        % If we're training an offline model save optimization result
        %b0 = filter.x';
        v0 = (ballRaw(:,end) - ballRaw(:,end-1)) ./ dt;
        b0 = [ballRaw(:,end);v0(:)]';
        
        %%{
        msg = [uint8(3), typecast(uint32(1),'uint8'), uint8(0)];
        data = typecast(msg, 'uint8'); 
        zmq.core.send(socket, data);
        response = zmq.core.recv(socket,bufferLength);

        % get q,q0
        STR = decodeResponseFromSL(response);
        qInit = STR.robot.traj.q;
        qdInit = STR.robot.traj.qd;
        tInit = STR.robot.traj.time;
        %}
        
        N = size(X,1);
        % find the closest point among Xs
        dif = repmat(b0,N,1) - X;
        [~,idx] = min(diag(dif*dif'));
        val = Y(idx,:);
        qf = val(1:7)';
        qfdot = val(7+1:2*7)';
        T = val(end);
        q0dot = zeros(7,1);
        Tret = 1.0;
        dt = 0.002;
        %[q,qd,qdd] = generateSpline(0.002,qInit,qdInit,qf,qfdot,T,Tret);
        [q,qd,qdd] = generateSpline(dt,q0,q0dot,qf,qfdot,T,Tret);
        %[q,qd,qdd] = wam.checkJointLimits(q,qd,qdd);
        %q = qInit * ones(1,length(q));  
        %qd = qdInit * ones(1,length(q));
        %qdd = zeros(7,1) * ones(1,length(q));
        
        
        %{
        figure;
        subplot(3,2,1);
        plot(q');
        subplot(3,2,3);
        plot(qd');
        subplot(3,2,5);
        plot(qdd');
        subplot(3,2,2);
        plot(q');
        subplot(3,2,4);
        plot(diff(q')./dt);
        subplot(3,2,6);
        plot(diff(diff(q'))./(dt)^2);
        %}

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
        response = zmq.core.recv(socket,bufferLength);
        toc
        pause(4.0);
        %pause(T+Tret);
        
        disp('Finished sending trj');
        
        msg = [uint8(5), uint8(0)];
        data = typecast(msg,'uint8');
        zmq.core.send(socket,data);
        response = zmq.core.recv(socket,bufferLength);

    end       
    
end

figure;
scatter3(ballRaw(1,:),ballRaw(2,:),ballRaw(3,:),'r');
hold on;
scatter3(ballFilt(1,:),ballFilt(2,:),ballFilt(3,:),'y');
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
zlim([table_z - tol_z, table_z + 2*tol_z]);
fill3(tennisTable(1:4,1),tennisTable(1:4,2),tennisTable(1:4,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
net_width = 0.01;

 
% disconnect from zmq and SL
disconnectFromSL(socket,address,context);
%}