%% Hit incoming ball from SL

tic
disp('Sending trj data');
% If we're training an offline model save optimization result
b0 = filter.x';        

% for debugging
%{
dtPred = 0.01;
[ballPred,ballTime,numBounce,time2PassTable] = ...
    predictBall(dtPred,predictTime,filter,table);
% land the ball on the centre of opponents court
desBall(1) = 0.0;
desBall(2) = table.DIST - 3*table.LENGTH/4;
desBall(3) = table.Z;
time2reach = 0.5; % time to reach desired point on opponents court
fast = true;
racketDes = calcRacketStrategy(desBall,ballPred,ballTime,time2reach,fast);
[qf,qfdot,T] = calcOptimalPoly(wam,racketDes,q0,Tret);
%}

%%{
N = size(X,1);
% find the closest point among Xs
dif = repmat(b0,N,1) - X;
[~,idx] = min(diag(dif*dif'));
b0
X(idx,:)
val = Y(idx,:);
qf = val(1:7)';
qfdot = val(7+1:2*7)';
T = val(end);
%}
q0dot = zeros(7,1);
dt = 0.002;

% remove expected delay 
TexpDelay = 0.2;
T = T - TexpDelay;

%[q,qd,qdd] = generateSpline(0.002,qInit,qdInit,qf,qfdot,T,Tret);
[q,qd,qdd] = generateSpline(dt,q0,q0dot,qf,qfdot,T,Tret);        
%[q,qd,qdd] = wam.checkJointLimits(q,qd,qdd);

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

%filter.initState(-inf(6,1),eps);