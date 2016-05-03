%% Strike ball from SL

% First get the ball and time and status 
% Start running the EKF and predict time2reach VHP
% Continue getting ball estimates until time2reach < maxTime

obj = onCleanup(@() disconnectFromSL(socket,address,context));
clc; clear all; close all;

%% Initialize Barrett WAM

initializeWAM;
[joint,ee,racket] = wam.drawPosture(q0);
endeff = [joint(end,:); ee];

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

% some useful colors
gray = [0.5020 0.5020 0.5020];
red = [1.0000    0.2500    0.2500];
plot3(joint(:,1),joint(:,2),joint(:,3),'k','LineWidth',10);
plot3(endeff(:,1),endeff(:,2),endeff(:,3),'Color',gray,'LineWidth',5);
fill3(racket(1,:), racket(2,:), racket(3,:),red);

%% Initialize EKF

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
dim = 6;
eps = 1e-6; %1e-3;
mats.O = eps * eye(dim);
mats.C = [eye(3),zeros(3)];
mats.M = eps * eye(3);
filter = EKF(dim,ballFlightFnc,mats);

% initialize the filters state with sensible values
guessBallInitVel = [-1.08; 4.80; 3.84]; %1.0 * [-0.80; 7.0; 2.0];
filter.initState([ball_cannon(:); guessBallInitVel],eps);

%% Get initial positioning of robot
% msg = [uint8(3), typecast(uint32(1),'uint8'), uint8(0)];
% data = typecast(msg, 'uint8'); 
% zmq.core.send(socket, data);
% response = zmq.core.recv(socket);
% 
% % get q,q0
% STR = decodeResponseFromSL(response);
% qInit = STR.robot.traj.q;
% qdInit = STR.robot.traj.qd;
% tInit = STR.robot.traj.time;
% 
% % Send to desired starting posture
% Qinit = [qInit;qdInit];
% dt = 0.002;
% tf = 1.0;
% p = generatePoly3rd(Qinit,Q0,dt,tf);
% % change last velocity and acc command to zero
% % p(8:end,end) = 0.0;
% timeSteps = size(p,2);
% ts = repmat(-1,1,timeSteps); % start immediately
% poly = [p;ts];
% poly = poly(:);
% poly = typecast(poly,'uint8');
% % 1 is for clear
% % 2 is for push back
% N = typecast(uint32(timeSteps),'uint8');
% poly_zmq = [uint8(1), uint8(2), N, poly', uint8(0)];
% data = typecast(poly_zmq, 'uint8');
% zmq.core.send(socket, data);
% response = zmq.core.recv(socket);
% pause(2.0);

%% Clear the ball positions
bufferLength = 1e6; %bytes
msg = [uint8(5), uint8(0)];
data = typecast(msg,'uint8');
zmq.core.send(socket,data);
response = zmq.core.recv(socket,bufferLength);

%% Load lookup table

% load the savefile
savefile = 'LookupTable.mat';
load(savefile,'X','Y');

%% Trajectory generation

ballObs = [];
ballTime = [];
ballRaw = [];
ballFilt = [];
curTime = 0.0;
lastBallPos = zeros(3,1);
lastBallTime = -inf;
j = 0;
firsttime = true;
numTrials = 0;
minBall2Predict = 5;
predictTime = 1.0;

table.DIST = dist_to_table;
table.LENGTH = table_length;
table.Z = table_z;
table.WIDTH = table_width;

% finite state machine
WAIT = 0;
PREDICT = 1;
FINISH = 2;
stage = WAIT;
Tret = 1.0;

while numTrials < 5

    observeBallFromSL;    
    filterBallFromSL;
    
    toly = 0.4;
    if size(ballRaw,2) > minBall2Predict && stage == WAIT && ...
            filter.x(2) > dist_to_table - table_length/2 && ...
            filter.x(2) < dist_to_table - table_length/2 + toly
        % otherwise predict
        %{
        dtPred = 0.01;
        [ballPred,~,numBounce,time2PassTable] = ...
            predictBall(dtPred,predictTime,filter,table);
        checkBounceOnOppTable(filter,table);
        %}
        stage = PREDICT;
    end

    % HIT THE BALL IF VALID
    if stage == PREDICT   
        stage = FINISH;
        numTrials = numTrials + 1;
        hitBallFromSL;
    end       
    
end

%% Plot the results

% seperate x into hit and return segments
% Nhit = floor(T/dt);
% 
% % load ball.mat
% load('ball.mat');
% [x,xd,o] = wam.calcRacketState([q;qd]);
% 
% scatter3(ballRaw(1,:),ballRaw(2,:),ballRaw(3,:),'r');
% scatter3(ballFilt(1,:),ballFilt(2,:),ballFilt(3,:),'y');
% scatter3(ballPred(1,:),ballPred(2,:),ballPred(3,:),'b');
% scatter3(x(1,1:Nhit),x(2,1:Nhit),x(3,1:Nhit),'r');
% scatter3(x(1,Nhit+1:end),x(2,Nhit+1:end),x(3,Nhit+1:end),'k');
% hold off;
 
%% Disconnect from zmq and SL
disconnectFromSL(socket,address,context);
%}