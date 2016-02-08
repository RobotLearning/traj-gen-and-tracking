%% Sending polynomials to SL with the VHP method

% First get the ball and time and status 
% Start running the EKF and predict time2reach VHP
% Continue getting ball estimates until time2reach < maxTime

clc; clear; close all;

%% Load table values

% load table parameters
loadTennisTableValues;

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

cleanup = onCleanup(@() disconnectFromSL(socket,address,context));

while true
    % Get the ball positions
    p = generatePoly3rd(Q0,Qf,dt,T);
end

start_time = -1.0; % start immediately
ts = repmat(start_time,1,time_steps);
poly = [p;ts];
poly = poly(:);
poly = typecast(poly,'uint8');

% 1 is for clear
% 2 is for push back
N = typecast(uint32(time_steps),'uint8');
poly_zmq = [uint8(1), uint8(2), N, poly', uint8(0)];
data = typecast(poly_zmq, 'uint8');
zmq.core.send(socket, data);
response = zmq.core.recv(socket);

