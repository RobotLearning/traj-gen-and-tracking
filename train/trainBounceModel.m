%% Train bounce model

clc; clear; close all;

% FOR ALL DATA
% Detect when it bounces
% Kalman smoothing on the first part
% Detect when it hits the racket
% Smooth on also this part
% Check if bouncing points are close
% Get velocities before and after impact
% Train the bounce model

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
params.yNet = dist_to_table - table_y;
params.table_length = table_length;
params.table_width = table_width;
% coeff of restitution-friction vector [TO BE LEARNED!]
params.CFTX = CFTX;
params.CFTY = CFTY;
params.CRT = CRT;
params.ALG = 'RK4';

% initial restitution-friction matrix 
M0 = diag([CFTX, CFTY, -CRT]);
M0full = diag([1,1,1,CFTX,CFTY,-CRT]);

funState = @(x,u,dt) discreteBallFlightModel(x,dt,params);
% very small but nonzero value for numerical stability
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);
filter = EKF(dim,funState,mats);

%% Load all data

demoFolder1 = '../Desktop/okanKinestheticTeachin_20141210/unifyData';
set1 = 15:65; % dataset of demonstrations
demoFolder2 = '../Desktop/okanKinestheticTeachin_20151124/unifyData';
set2 = 2:64;
dropSet1 = [17,20,22,23,24,27,29,...
            31,32,33,34,35,38,41,48,49,53,54,56]; % bad recordings
% 27,33,49 have outliers near table_z, these could be removed
dropSet2 = [5,6,9,10,14,19,26,27,32,38,42,55,56,57,60,62];
% 60 caused NaN for some reason in filtering
% 19 does not have good bounce data
%set = setdiff(set1,dropSet1);
set = setdiff(set2,dropSet2);
scale  = 0.001; % recorded in milliseconds

% bouncing data to be used for regression
velBeforeBounce = zeros(3,length(set));
velAfterBounce = zeros(3,length(set));

for idx = 1:length(set)

    %% Get rid of outliers
    M = dlmread([demoFolder2,int2str(set(idx)),'.txt']);
    %M = dlmread([demoFolder1,int2str(set(idx)),'.txt']);
    robot = M(:,end-20:end);
    %robot = M(:,end-14:end);
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
    
    %% Get ball positions till bounce

    % find the last ball position before first bounce
    % detect first bounce
        
    tol = 5e-2;    
    idxBallBounce = find(b(:,4) <= table_z + ball_radius + tol);
    idxJump = find(diff(idxBallBounce) ~= 1);
    if isempty(idxJump)
        [~,idxBounceCandidate] = min(b(idxBallBounce,4));
        idxBallBounce = idxBallBounce(idxBounceCandidate)-1;
    else
        idxBallBounce = idxBallBounce(1:idxJump(1));
        [~,idxBounceCandidate] = min(b(idxBallBounce,4));
        idxBallBounce = idxBallBounce(idxBounceCandidate)-1;
    end
    bTillTable = b(1:idxBallBounce,1:4);
    
    %% Fit a Kalman smoother
    obsLen = size(bTillTable,1);
    % initialize the state with first observation and first derivative est.
    p0 = bTillTable(1,2:4);
    v0 = (bTillTable(5,2:4) - bTillTable(1,2:4))/(bTillTable(5,1)-bTillTable(1,1));
    filter.initState([p0(:);v0(:)],eps);
    u = zeros(1,obsLen);
    tTillTable = bTillTable(:,1);
    [xEKFSmooth, VekfSmooth] = filter.smooth(tTillTable,bTillTable(:,2:4)',u);
    ballEKFSmoothPreBounce = C * xEKFSmooth;
    
    %% Predict if necessary to get velocity before bounce
    zBeforeBounce = ballEKFSmoothPreBounce(3,end);
    dt = 0.001;
    tol = 1e-2;
    i = 1;
    while zBeforeBounce >= table_z + ball_radius + tol
        filter.predict(dt,0);
        ballEKFSmoothPreBounce(:,obsLen+i) = C * filter.x;
        zBeforeBounce = filter.x(3);
        i = i + 1;
    end
    tBounce = bTillTable(end,1) + dt * (i-1);
    
    % Get y coordinate at bounce
    yAtBounce = filter.x(2);
    % Get velocity before bounce
    velBeforeBounce(:,idx) = filter.x(4:6);
    
    %% Estimate striking time
    % Get cartesian coordinates of robot trajectory
    [x,xd,o] = wam.kinematics(Q');

    % find closest points in space as a first estimate
    diffBallRobot = b(:,2:4)-x';
    distBallRobot = diag(diffBallRobot * diffBallRobot');
    [minDist,idxStrike] = min(distBallRobot);
    tStrike = b(idxStrike,1);
    fprintf('%d. Striking time est. for demo %d: %f\n',idx,set(idx),tStrike);
    
    %% Get ball positions after bounce till strike
    % idxNet = 1;
    % while b(idxNet,3) <= dist_to_table - table_y
    %     idxNet = idxNet + 1;
    % end
    idxDiffBallPos = [];
    j = 1; % time from bounce till strike
    % in case returning ball data is not available
    if size(b,1) == idxStrike, idxStrike = idxStrike - 1; end
    for i = idxBallBounce:idxStrike
        if norm(b(i,2:4) - b(i+1,2:4)) > tol
            idxDiffBallPos(j) = i;
            j = j+1;
        end
    end
    bAfterBounceTillStrike = b(idxDiffBallPos,:);
    
    %% Fit again a Kalman smoother
    obsLen = size(bAfterBounceTillStrike,1);
    % initialize the state with first observation and first derivative est.
    p0 = ballEKFSmoothPreBounce(1:3,end);
    v0 = M0 * xEKFSmooth(4:end,end);
    filter.initState([p0(:);v0(:)],M0full*VekfSmooth(:,:,end)*M0full);
    u = zeros(1,obsLen);
    tFromTableTillStrike = bAfterBounceTillStrike(:,1);
    [xEKFSmooth, VekfSmooth] = filter.smooth(tFromTableTillStrike,...
                                         bAfterBounceTillStrike(:,2:4)',u);
    ballEKFSmoothAfterBounce = C * xEKFSmooth;
    
    %% Predict to get velocity after bounce
    zAfterBounce = xEKFSmooth(3,1);
    yAfterBounce = xEKFSmooth(2,1);
    dt = 0.001;
    i = 1;
    filter.initState(xEKFSmooth(:,1),VekfSmooth(:,:,1));
    while zAfterBounce > table_z + ball_radius + tol && yAfterBounce > yAtBounce
        filter.predict(-dt,0);
        ballEKFSmoothAfterBounce = [C * filter.x, ballEKFSmoothAfterBounce];
        zAfterBounce = filter.x(3);
        yAfterBounce = filter.x(2);
    end
    
    % Get velocity before bounce
    yAtBounce2 = filter.x(2);
    velAfterBounce(:,idx) = filter.x(4:6);
    ballEKFSmooth = [ballEKFSmoothPreBounce, ballEKFSmoothAfterBounce];
    
    
end

%% Regression

% For comparing/plotting dataset1 with dataset2 together
%{
velBeforeBounce1 = velBeforeBounce;
velAfterBounce1 = velAfterBounce;
coeff1 = coeff;
load('velExp2.mat');
coeff2 = coeff;
velBeforeBounceTotal = [velBeforeBounce1, velBeforeBounce];
velAfterBounceTotal = [velAfterBounce1, velAfterBounce];
velBeforeBounce = velBeforeBounceTotal;
velAfterBounce = velAfterBounceTotal;
%}

coeff = zeros(1,3);
for i = 1:3
    x = velBeforeBounce(i,:)';
    y = velAfterBounce(i,:)';
    figure;
    scatter(x,y,'r');
    hold on;
    title('Velocity after vs. before');
    grid on;
    axis equal;
    coeff(i) = x \ y;
    plot(x,coeff(i)*x,'b');
    hold off;
end

% figure(4);
% scatter3(bx1,by1,bz1,'r');
% hold on;
% scatter3(bx3,by3,bz3,'b');
% %scatter3(ballEKF(1,:),ballEKF(2,:),ballEKF(3,:));
% scatter3(ballEKFSmooth(1,:),ballEKFSmooth(2,:),ballEKFSmooth(3,:));
% title('Filtered ball trajectory');
% grid on;
% axis equal;
% xlabel('x');
% ylabel('y');
% zlabel('z');
% legend('cam1','cam3','EKF smooth');
% fill3(T(:,1),T(:,2),T(:,3),[0 0.7 0.3]);
% fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
% hold off;