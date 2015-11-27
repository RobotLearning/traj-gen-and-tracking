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

table_z = floor_level - table_height;
table_x = table_center + table_width/2;
table_y = table_length/2;

T1 = [table_center - table_x; 
    dist_to_table - table_length; 
    table_z];
T2 = [table_center + table_x;
    dist_to_table - table_length;
    table_z];
T3 = [table_center + table_x;
    dist_to_table;
    table_z];
T4 = [table_center - table_x;
    dist_to_table;
    table_z];
T = [T1,T2,T3,T4]';

net1 = [table_center - table_x;
        dist_to_table - table_y;
        table_z];
net2 = [table_center + table_x;
        dist_to_table - table_y;
        table_z];
net3 = [table_center + table_x;
        dist_to_table - table_y;
        table_z + net_height];
net4 = [table_center - table_x;
        dist_to_table - table_y;
        table_z + net_height];

net = [net1,net2,net3,net4]';

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

% initial restitution-friction matrix 
M0 = diag([CFTX, CFTY, CRT]);
M0full = diag([1,1,1,CFTX,CFTY,CRT]);

funState = @(x,u,dt) symplecticFlightModel(x,dt,params);
% very small but nonzero value for numerical stability
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);
filter = EKF(dim,funState,mats);

%% Load all data

demoFolder = '../Desktop/okanKinestheticTeachin_20141210/unifyData';
set = 15:65; % dataset of demonstrations
dropSet1 = [17,20,22,23,24,29,31,32,34,35,38,41,48,53,54,56]; % bad recordings
for i = 1:length(dropSet1)
    set = set(set ~= dropSet1(i));
end
scale  = 0.001; % recorded in milliseconds

% bouncing data to be used for regression
bounce = [];
idxBnc = 0;

for i = 15 %set

    M = dlmread([demoFolder,int2str(i),'.txt']);
    t = scale * M(:,end-14);

    % this is to throw away really bad observations as to not confuse the
    % filter
    zMax = 0.5;
    % camera 1 time indices for reliable observations
    idxCam1 = find(M(:,2) == 1 & M(:,5) >= table_z & M(:,5) < zMax);
    tCam1 = t(idxCam1);
    Mball1 = M(:,3:5);
    bx1 = Mball1(idxCam1,1);
    by1 = Mball1(idxCam1,2);
    bz1 = Mball1(idxCam1,3);

    % camera 3 time indices for reliable observations
    idxCam3 = find(M(:,7) == 1 & M(:,10) >= table_z & M(:,10) < zMax);
    tCam3 = t(idxCam3);
    Mball3 = M(:,8:10);
    bx3 = Mball3(idxCam3,1);
    by3 = Mball3(idxCam3,2);
    bz3 = Mball3(idxCam3,3);
    
    b = [tCam1,bx1,by1,bz1,ones(length(tCam1),1); % camera 1
         tCam3,bx3,by3,bz3,3*ones(length(tCam3),1)]; % camera 3
    [b,idxB] = sortrows(b);
    
    robot = M(:,end-14:end);
    q = robot(:,2:N_DOFS+1);
    qd = robot(:,N_DOFS+2:N_DOFS+8);
    Q = [q';qd'];
    Qb = Q(:,idxB); % robot joint angles and vel. sorted via ball observ.
    
    %% Get ball positions till net
    % idxNet = 1;
    % while b(idxNet,3) <= dist_to_table - table_y
    %     idxNet = idxNet + 1;
    % end
    tol = 1e-3;
    idxBallBounce = find(b(:,4) <= table_z + ball_radius + tol);
    idxBallBounce = idxBallBounce(1);
    j = 1; % time from ballgun to net 
    for i = 1:idxBallBounce
        if norm(b(i,2:4) - b(i+1,2:4)) > tol
            idxDiffBallPos(j) = i;
            j = j+1;
        end
    end
    bTillTable = b(idxDiffBallPos,:);
    
    %% Fit a Kalman smoother
    obsLen = size(bTillTable,1);
    % initialize the state with first observation and first derivative est.
    p0 = bTillTable(1,2:4);
    v0 = (bTillTable(5,2:4) - bTillTable(1,2:4))/(bTillTable(5,1)-bTillTable(1,1));
    filter.initState([p0(:);v0(:)],eps);
    u = zeros(1,obsLen);
    tTillTable = bTillTable(:,1);
    [xEKFSmooth, VekfSmooth] = filter.smooth(tTillTable,bTillTable(:,2:4)',u);
    ballEKFSmooth = C * xEKFSmooth;
    
    %% Predict if necessary to get velocity before bounce
    zBeforeBounce = ballEKFSmooth(3,end);
    dt = 0.001;
    i = 1;
    while zBeforeBounce >= table_z + ball_radius + tol
        filter.predict(dt,0);
        ballEKFSmooth(:,obsLen+i) = C * filter.x;
        zBeforeBounce = filter.x(3);
        i = i + 1;
    end
    tBounce = bTillTable(end,1) + dt * (i-1);
    
    % Get velocity before bounce
    velBeforeBounce = filter.x(4:6);
    
    %% Estimate striking time
    % Get cartesian coordinates of robot trajectory
    x = wam.kinematics(Qb);

    % find closest points in space as a first estimate
    diffBallRobot = b(:,2:4)-x';
    distBallRobot = diag(diffBallRobot * diffBallRobot');
    [minDist,idxStrike] = min(distBallRobot);
    tStrike = b(idxStrike,1);
    fprintf('Striking time est. for demo %d: %f\n',i,tStrike);
    
    %% Get ball positions after net till strike
    % idxNet = 1;
    % while b(idxNet,3) <= dist_to_table - table_y
    %     idxNet = idxNet + 1;
    % end
    idxDiffBallPos = [];
    j = 1; % time from bounce till strike
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
    p0 = ballEKFSmooth(1:3,end);
    v0 = M0 * xEKFSmooth(4:end,end);
    filter.initState([p0(:);v0(:)],M0full*VekfSmooth(:,:,end)*M0full);
    u = zeros(1,obsLen);
    tFromTableTillStrike = bAfterBounceTillStrike(:,1);
    [xEKFSmooth, VekfSmooth] = filter.smooth(tFromTableTillStrike,...
                                         bAfterBounceTillStrike(:,2:4)',u);
    ballEKFSmooth = [ballEKFSmooth, C * xEKFSmooth];
    
    %% Predict to get velocity before bounce
    zAfterBounce = ballEKFSmooth(3,end);
    dt = 0.001;
    i = 1;
    filter.initState(xEKFSmooth(:,1),VekfSmooth(:,:,1));
    while zAfterBounce >= table_z + ball_radius + tol
        filter.predict(-dt,0);
        ballEKFSmooth(:,end+1) = C * filter.x;
        zAfterBounce = filter.x(3);
        i = i + 1;
    end
    
    % Get velocity before bounce
    velAfterBounce = filter.x(4:6);
    
    
end

%% Plot smoothing
figure(3);
scatter3(bx1,by1,bz1,'r');
hold on;
scatter3(bx3,by3,bz3,'b');
%scatter3(ballEKF(1,:),ballEKF(2,:),ballEKF(3,:));
scatter3(ballEKFSmooth(1,:),ballEKFSmooth(2,:),ballEKFSmooth(3,:));
title('Filtered ball trajectory');
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
legend('cam1','cam3','EKF smooth');
fill3(T(:,1),T(:,2),T(:,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
hold off;