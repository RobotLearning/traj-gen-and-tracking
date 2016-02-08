%% Test sending polynomials to SL

clc; clear; close all;

host = 'localhost'; 
port = '7646';
address = sprintf('tcp://%s:%s',host,port);
context = zmq.core.ctx_new();
socket  = zmq.core.socket(context, 'ZMQ_REQ');
zmq.core.connect(socket, address);

msg = [uint8(3), typecast(uint32(1),'uint8'), uint8(0)];
data = typecast(msg, 'uint8'); 
zmq.core.send(socket, data);
response = zmq.core.recv(socket);

% get q,q0
STR = decodeResponseFromSL(response);
q0 = STR.robot.traj.q;
qd0 = STR.robot.traj.qd;
t0 = STR.robot.traj.time;

qdes = [1.7967, -0.2344, -0.0941, 1.7391, -1.5738, 0.0445, 0.2699];
qdotdes = zeros(7,1);
time_steps = 8000;
robot_freq = 500;
T = time_steps / robot_freq;
dt = 1/robot_freq;

Q0 = [q0(:);qd0(:)];
Qf = [qdes(:);qdotdes];
p = generatePoly3rd(Q0,Qf,dt,T);
%p = p(1:14,:); %not including qdd for now 

start_time = -1.0;
%t0 = t0 + start_time;
%t = t0+dt:dt:t0+T;
%poly = [p;t];
poly = [p;repmat(start_time,1,time_steps)];
poly = poly(:);
poly = typecast(poly,'uint8');

% 1 is for clear
% 2 is for push back
N = typecast(uint32(time_steps),'uint8');
poly_zmq = [uint8(1), uint8(2), N, poly', uint8(0)];
data = typecast(poly_zmq, 'uint8');
zmq.core.send(socket, data);
response2 = zmq.core.recv(socket);

zmq.core.disconnect(socket, address);
zmq.core.close(socket);
zmq.core.ctx_shutdown(context);
zmq.core.ctx_term(context);