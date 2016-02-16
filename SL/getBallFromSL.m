%% Get ball positions from SL via matlab-zmq

obj = onCleanup(@() disconnectFromSL(socket,address,context));
clc; clear all; close all;

%% Create the socket
host = 'localhost'; 
port = '7646';
address = sprintf('tcp://%s:%s',host,port);
context = zmq.core.ctx_new();
socket  = zmq.core.socket(context, 'ZMQ_REQ');
zmq.core.connect(socket, address);

%% Clear the ball positions
% bufferLength = 1e6; %bytes
% msg = [uint8(5), uint8(0)];
% data = typecast(msg,'uint8');
% zmq.core.send(socket,data);
% response = zmq.core.recv(socket,bufferLength);

%% GET BALL POSITIONS

ballObs = [];
ballTime = [];
j = 1;

for i = 1:1000

    msg = [uint8(4), typecast(uint32(1),'uint8'),uint8(0)];
    data = typecast(msg,'uint8');
    zmq.core.send(socket,data);
    response = zmq.core.recv(socket);
    STR = decodeResponseFromSL(response);
    ballObs = STR.ball.pos;
    ballTime = STR.ball.time;
    
    if ~isempty(ballTime)
        % get number of observations
        numObs = length(ballTime);
        ball(:,j:j+numObs-1) = [ballObs;ballTime]; 
        j = j + numObs;
    end
end

scatter3(ball(1,:),ball(2,:),ball(3,:));

% disconnect from zmq and SL
disconnectFromSL(socket,address,context);