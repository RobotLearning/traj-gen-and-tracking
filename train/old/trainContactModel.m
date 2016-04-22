%% Train contact model

clc; clear; close all;

% FOR ALL DATA
% Detect roughly the striking time first
% Then segment the ball motion into incoming/return traj.
% Apply kalman smoothing on both parts of ball traj
% Apply kalman smoothing on robot traj
% Try to refine the striking time
% Build cartesian racket vel and orientation using Jacobian
% Compare with saved cartesian vel. values
% Build a model of contact: (rdot,orient,bdot_in) -> bdot_out
% Compare with the existing model of contact and consider regressing
% on top of it

%% Load table values

% load table parameters
loadTennisTableValues;

%% Initialize Barrett WAM

initializeWAM;

%% initialize EKF
dim = 6;
eps = 1e-6;
C = [eye(3),zeros(3)];

params.C = Cdrag;
params.g = gravity;
params.zTable = table_z;
params.radius = ball_radius;
params.yNet = dist_to_table - table_y;
params.table_length = table_length;
params.table_width = table_width;
% coeff of restitution-friction vector [TO BE LEARNED!]
params.CFTX = CFTX;
params.CFTY = CFTY;
params.CRT = CRT;
params.ALG = 'RK4';

funState = @(x,u,dt) discreteBallFlightModel(x,dt,params);
% very small but nonzero value for numerical stability
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);

filterPre = EKF(dim,funState,mats);
filterAfter = EKF(dim,funState,mats);

%% Load all data

demoFolder1 = '../Desktop/okanKinestheticTeachin_20141210/unifyData';
set1 = 15:65; % dataset of demonstrations
demoFolder2 = '../Desktop/okanKinestheticTeachin_20151124/unifyData';
set2 = 2:64;
dropSet1 = [17,20,22,23,24,27,29,...
            31,32,33,34,35,38,41,48,49,53,54,56]; % bad recordings
% 27,33,49 have outliers near table_z, these could be removed
dropSet2 = [5,6,9,10,14,19,26,27,32,38,42,55,56,57,62];
% 19 does not have good bounce data
set = setdiff(set1,dropSet1);
scale  = 0.001; % recorded in milliseconds

% bouncing data to be used for regression
velBallPreStrike = zeros(3,length(set));
velBallAfterStrike = zeros(3,length(set));

for idx = 2 %1:length(set)

    %% Get rid of outliers
    M = dlmread([demoFolder1,int2str(set(idx)),'.txt']);
    robot = M(:,end-14:end);
    %robot = M(:,end-20:end);
    q = robot(:,2:N_DOFS+1);
    qd = robot(:,N_DOFS+2:N_DOFS+8);
    t = scale * robot(:,1);

    % this is to throw away really bad observations as to not confuse the
    % filter
    zMax = 0.5;
    % camera 1 time indices for reliable observations
    idxCam1 = find(M(:,2) == 1 & M(:,5) >= table_z & M(:,5) < zMax);
    tCam1 = t(idxCam1);
    qCam1 = q(idxCam1,:);
    qdCam1 = qd(idxCam1,:);
    Mball1 = M(:,3:5);
    bx1 = Mball1(idxCam1,1);
    by1 = Mball1(idxCam1,2);
    bz1 = Mball1(idxCam1,3);

    % camera 3 time indices for reliable observations
    idxCam3 = find(M(:,7) == 1 & M(:,10) >= table_z & M(:,10) < zMax);
    tCam3 = t(idxCam3);
    qCam3 = q(idxCam3,:);
    qdCam3 = qd(idxCam3,:);
    Mball3 = M(:,8:10);
    bx3 = Mball3(idxCam3,1);
    by3 = Mball3(idxCam3,2);
    bz3 = Mball3(idxCam3,3);
    
    b = [tCam1,bx1,by1,bz1,qCam1,qdCam1,ones(length(tCam1),1); % camera 1
         tCam3,bx3,by3,bz3,qCam3,qdCam3,3*ones(length(tCam3),1)]; % camera 3
    [b,idxB] = sortrows(b);
    
    %% Remove same ball positions
    j = 1; % time from ballgun to net 
    tol = 1e-3;
    idxDiffBallPos = [];
    for i = 1:size(b,1)-1
        if norm(b(i,2:4) - b(i+1,2:4)) > tol
            idxDiffBallPos(j) = i;
            j = j+1;
        end
    end
    b = b(idxDiffBallPos,:);
    
    % get the robot positions that correspond to ball positions in t
    Q = b(:,5:5+2*N_DOFS-1); % robot joint angles and vel. sorted via ball observ.
    
    %% Estimate striking time roughly
    % Get cartesian coordinates of robot trajectory
    [racketPos,racketVel,orient] = wam.kinematics(Q');

    % find closest points in space as a first estimate
    diffBallRobot = b(:,2:4)-racketPos';
    distBallRobot = diag(diffBallRobot * diffBallRobot');
    [minDist,idxStrike] = min(distBallRobot);
    tStrike = b(idxStrike,1);
    fprintf('%d. Striking time est. for demo %d: %f\n',idx,set(idx),tStrike);
    
    %% Fit a Kalman smoother till strike
    bTillStrike = b(1:idxStrike,1:4);
    obsLen = size(bTillStrike,1);
    % initialize the state with first observation and first derivative est.
    p0 = bTillStrike(1,2:4);
    v0 = (bTillStrike(5,2:4) - bTillStrike(1,2:4))/(bTillStrike(5,1)-bTillStrike(1,1));
    filterPre.initState([p0(:);v0(:)],eps);
    u = zeros(1,obsLen);
    tTillStrike = bTillStrike(:,1);
    [xEKFSmoothPre, VekfSmoothPre] = filterPre.smooth(tTillStrike,bTillStrike(:,2:4)',u);
    ballEKFSmoothPreStrike = C * xEKFSmoothPre;
    
    %% Update the striking time estimate
    
    windowTime = 0.2;
    dt = 0.01;
    timePossibleStrikes =  b(tTillStrike > timeStrike - windowTime/2);
    
    diff = obj.pos - racketPos;
    distToRacketPlane = racketNormal'*diff;
    distOnRacketPlane = sqrt(diff'*diff - distToRacketPlane^2);
    
    racketRot = quat2Rot(orient);
    racketNormal = racketRot(:,3);
    
    idxStrike = find(distOnRacketPlane < racket_radius & distNextToRacketPlane < ball_radius,1);
    if distOnRacketPlane < racket_radius && distNextToRacketPlane < 0
    
    %% Fit again a Kalman smoother from strike onwards
    
    bAfterStrike = b(b(:,1) >= tStrike & b(:,end) == 3,1:4);
    obsLen = size(bAfterStrike,1);
    % initialize the state with first observation and first derivative est.
    p0 = ballEKFSmoothPreStrike(:,end);
    % Calculate outgoing velocity using conservation of momentum
    velBallPreStrike(:,idx) = xEKFSmoothPre(4:end,end);
    % Get the racket prenormal
    quatRacketStrike = orient(:,idxStrike);
    R = quat2Rot(quatRacketStrike);
    racketNormal = R(:,3); % of unit length
    % Get the racket velocity
    velRacket = racketVel(:,idxStrike);
    velBallPreStrikeNormal = velBallPreStrike(:,idx)' * racketNormal;
    velBallAlongRacket = velBallPreStrike(:,idx) - velBallPreStrikeNormal .* racketNormal;
    velRacketAlongNormal = velRacket' * racketNormal;
    velBallAfterStrikeNormal = velRacketAlongNormal + ...
                     CRR * (velRacketAlongNormal - velBallPreStrikeNormal);
    % assuming velocity along racket stays the same
    velBallAfterStrike(:,idx) = velBallAlongRacket + velBallAfterStrikeNormal .* racketNormal;
    v0 = velBallAfterStrike(:,idx);
    filterAfter.initState([p0(:);v0(:)],VekfSmoothPre(:,:,end));
    u = zeros(1,obsLen);
    tAfterStrike = bAfterStrike(:,1);
    [xEKFSmoothAfter, VekfSmoothAfter] = filterAfter.smooth(tAfterStrike,...
                                             bAfterStrike(:,2:4)',u);
    ballEKFSmoothAfterStrike = C * xEKFSmoothAfter;
    
    %% Refine the striking time
    
    % find the times when the ball path and racket path are within
    % radius cms. and then take the first one
    % use prediction models for the first 
    
    % Predict till both parts coincide
    xPre = xEKFSmoothPre(:,end);
    ballPre = C * xPre;
    xAfter = xEKFSmoothAfter(:,1);
    ballAfter = C * xAfter;
    dt = 0.01;
    i = 1;
    filterPre.initState(xPre,VekfSmoothPre(:,:,end));
    filterAfter.initState(xAfter,VekfSmoothAfter(:,:,1));
    tol = 1e-3;
    
    % FIRST PREDICT BACK TILL TSTRIKE
    % THEN PREDICT TOGETHER WITH BOTH FILTERS
    
    while norm(ballPre - ballAfter,2) > tol
        filterAfter.predict(-dt,0);
        ballAfter = C * filterAfter.x;
        filterPre.predict(dt,0);
        ballPre = C * filterPre.x;
        i = i+1;
    end
    
    % Refine tStrike
    tStrike = tStrike + (i-1)*dt;
    % Get the velocities
    
    % update the velocity before strike
    velBallPreStrike(:,idx) = filterPre.x(4:end);
    % update the velocity after strike
    velBallAfterStrike(:,idx) = filterAfter.x(4:end);
    
end

%% Regression
ballEKFSmooth = [ballEKFSmoothPreStrike, ballEKFSmoothAfterStrike];

figure(4);
scatter3(bx1,by1,bz1,'r');
hold on;
scatter3(bx3,by3,bz3,'b');
%scatter3(ballEKF(1,:),ballEKF(2,:),ballEKF(3,:));
%scatter3(ballEKFSmooth(1,:),ballEKFSmooth(2,:),ballEKFSmooth(3,:));
title('Filtered ball trajectory');
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
%legend('cam1','cam3','EKF smooth');
fill3(T(:,1),T(:,2),T(:,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
hold off;