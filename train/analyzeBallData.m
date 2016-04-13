%% AnalyzeBallData

% analyzes the following prediction scenarios:
% 
% 1. prediction with nonlinear least squares [NLS] + 
% nonlinear model till bounce and striking time
% 
% 2. prediction with (total) linear least squares
% 
% 3. prediction with Extended Kalman Filter with initial 
% ball position and velocity calculated with NLS

clc; clear; close all;

%% Load table values

% load table parameters
loadTennisTableValues;
% Initialize Barrett WAM
initializeWAM;

%% Load last data

demoFolder1 = '../Desktop/okanKinestheticTeachin_20141210/unifyData';
set1 = 15:65; % dataset of demonstrations
demoFolder2 = '../Desktop/okanKinestheticTeachin_20151124/unifyData';
set2 = 2:64;
dropSet1 = [17,20,22,23,24,29,31,32,34,35,38,41,48,53,54,56]; % bad recordings
dropSet2 = [5,6,9,10,14,26,27,32,38,42,55,56,57];
% what's wrong with sample 31? problem on the ekf-smoother bouncing 
% what's wrong with sample 54? symplecticFlightModel does not terminate
%set = setdiff(set1,dropSet1);
demo = demoFolder1;

% adapt between different datasets
if strcmp(demo,demoFolder1)
    dist_to_table = -0.80;
    loadTennisTableValues();
    IDX = 14;
else 
    dist_to_table = -1.15;
    IDX = 20;
end

choose = 16; %set(end);
M = dlmread([demo,int2str(choose),'.txt']);
scale = 0.001; % recorded in milliseconds
t = scale * M(:,11);
% Second dataset includes cartesian positions so idx is 20 then
robot = M(:,end-IDX:end);
dof = N_DOFS;
q = robot(:,2:dof+1);
qd = robot(:,dof+2:dof+8);
Q = [q';qd'];
% Get cartesian coordinates of robot trajectory
[x,xd,o] = wam.kinematics(Q);

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

% Estimate striking time
camInfo = [tCam1,bx1,by1,bz1;
            tCam3,bx3,by3,bz3];

% find closest points in space as a first estimate
distBallRobotSqr = (bx1-x(1,idxCam1)').^2 + (by1-x(2,idxCam1)').^2 + (bz1-x(3,idxCam1)').^2;
[minDist,idx] = min(distBallRobotSqr);
tStrike = camInfo(idx,1);
fprintf('Striking time est. for demo %d: %f\n',choose,tStrike);

%% initialize EKF
dim = 6;
eps = 1e-6;
C = [eye(3),zeros(3)];

Cdrag = 0.1476;
gravity = -10.37;

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
params.ALG = 'RK4';

funState = @(x,u,dt) discreteBallFlightModel(x,dt,params);
% very small but nonzero value for numerical stability
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);
filter = EKF(dim,funState,mats);

%% Order the ball data w.r.t. time and throw away same values

useCam1 = true;

if useCam1
    b = [tCam1,bx1,by1,bz1,ones(length(tCam1),1); % camera 1
         tCam3,bx3,by3,bz3,3*ones(length(tCam3),1)]; % camera 3
else
    b = [tCam3,bx3,by3,bz3,ones(length(tCam3),1)];
end
b = sortrows(b);

j = 1; % time from ballgun to net 
tol = 1e-3;
idxDiffBallPos = []; % index for diff ball positions
for i = 1:size(b,1)-1
    if norm(b(i,2:4) - b(i+1,2:4)) > tol
        idxDiffBallPos(j) = i;
        j = j+1;
    end
end

b = b(idxDiffBallPos,:);

% find the last ball position before net 
idxBallNet = find(b(:,3) <= dist_to_table - table_length/2);
% if there is a returning trajectory
idxJump = find(diff(idxBallNet) ~= 1);
if ~isempty(idxJump)
    idxBallNet = idxBallNet(1:idxJump(1));
end
bTillNet = b(idxBallNet,1:4);
tPredict = tStrike - bTillNet(end,1);
tPredictNLS = tStrike - bTillNet(1,1);
bTillNet(:,1) = bTillNet(:,1) - bTillNet(1,1);

% Get ball positions till bounce

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
bTillTable(:,1) = bTillTable(:,1) - bTillTable(1,1);
bBounce = bTillTable(end,1:4);
tBounce = bBounce(1) - bTillNet(end,1);
tBounceNLS = bBounce(1) - bTillNet(1,1);

% checking something
%bTillNet(:,1) = 1/60 * (1:size(bTillNet,1))';

%% Estimate initial state and predict with NLS

% % make sure t starts from 0
bTillNet(:,1) = bTillNet(:,1) - bTillNet(1,1);
guessBallInitVel = 1.0 * [-0.80; 7.0; 2.0];
ballInitNLS = [bTillNet(1,2:4)';guessBallInitVel];
% Run nonlinear least squares to estimate ballInit better
x0 = [ballInitNLS(:);Cdrag;gravity];
ballData = bTillNet;
len = size(ballData,1);
ballFun = @(b0,C,g) predictNextBall(b0,C,g,ballData,len);
fnc = @(x) ballFun(x(1:6),x(7),x(8));
tic;
x = lsqnonlin(fnc,x0);
toc
ballInitNLS = x(1:6);
Cdrag_est = x(7);
gravity_est = x(8);
p0 = ballInitNLS(1:3);
v0 = ballInitNLS(4:6);
% 
params.C = Cdrag;
params.g = gravity;
funState2 = @(x,u,dt) discreteBallFlightModel(x,dt,params);
filter2 = EKF(dim,funState2,mats);
filter2.initState([p0(:);v0(:)],eps);
dt = 0.01;
N_pred_NLS = tPredictNLS/dt;
N_bounce = round(tBounceNLS/dt);
for i = 1:N_pred_NLS
    filter2.linearize(dt,0);
    filter2.predict(dt,0);
    ballNLS(:,i) = C * filter2.x;
    if i == N_bounce
        % save position to compare later
        bBouncePredict = ballNLS(:,i);
    end
end

%% LS prediction based on projectile motion

%%{
tPredictLS = tPredictNLS;
len = size(bTillNet,1);
M = [ones(len,1),bTillNet(:,1)];
Y = bTillNet(:,2:3);
Yz = bTillNet(:,4) - 0.5*gravity*bTillNet(:,1).^2;
Y = [Y,Yz];
tol = 0.05;
ballInitLS(:,1) = tls(M,Y(:,1)); %pinv(M,tol)*Y;
ballInitLS(:,2) = tls(M,Y(:,2));
ballInitLS(:,3) = tls(M,Y(:,3));

dt = 0.01;
N_pred_LS = round(tPredictLS/dt);
tLS = dt * (1:N_pred_LS);
Mstar = [ones(N_pred_LS,1),tLS(:)];
ballLS = Mstar * ballInitLS + [zeros(N_pred_LS,2),0.5*gravity*tLS(:).^2];
ballLLS = ballLS';
%}

%% Predict the ball with EKF till striking time

obsLen = size(bTillNet,1);
% initialize the state with first observation and first derivative est.
%p0 = bTillNet(1,2:4);
%guessBallInitVel = 1.0 * [-0.80; 7.0; 2.0];
%v0 = guessBallInitVel;
% p0 = ballInit(1,:);
% v0 = ballInit(2,:);
initVar = 1;
filter.initState([p0(:);v0(:)],initVar);
for i = 1:obsLen-1
    dt = bTillNet(i+1,1) - bTillNet(i,1);
    filter.linearize(dt,0);
    filter.update(bTillNet(i,2:4)',0);
    ballEKF(:,i) = C * filter.x;
    filter.predict(dt,0);
end
filter.update(bTillNet(obsLen,2:4)',0);
ballEKF(:,obsLen) = C * filter.x;

% Now we start predicting ball till striking time
dt = 0.01;
N_pred = tPredict/dt;
N_bounce = round(tBounce/dt);
% after bounce there should be less topspin effects
gravityChange = false;

for i = 1:N_pred
    %filter.linearize(dt,0);
    filter.predict(dt,0);
    ballEKF(:,obsLen+i) = C * filter.x;
    if i == N_bounce
        % save position to compare later
        bBouncePredict = ballEKF(:,obsLen+i);
    end
    % try changing constants
    if filter.x(6) > 0 && ~gravityChange
        gravity = -9.80;
        params.g = gravity;
        funAfterBounce = @(x,u,dt) discreteBallFlightModel(x,dt,params);
        filter.f = funAfterBounce;
        gravityChange = true;
    end
        
end

%% Smoothen the ball path EKF smoother to get x0

filter.initState([p0(:);v0(:)],eps);
u = zeros(1,size(bTillNet,1));
tTillNet = bTillNet(:,1);
[xEKFSmooth, VekfSmooth] = filter.smooth(tTillNet,bTillNet(:,2:4)',u);
ballEKFSmooth = C * xEKFSmooth;


%% Plot predictions
figure(3);
s1 = scatter3(bx1,by1,bz1,'r');
hold on;
s2 = scatter3(bx3,by3,bz3,'b');
%s3 = scatter3(ballEKF(1,:),ballEKF(2,:),ballEKF(3,:));
s3 = scatter3(ballNLS(1,:),ballNLS(2,:),ballNLS(3,:),'y');
%scatter3(ballEKFSmooth(1,:),ballEKFSmooth(2,:),ballEKFSmooth(3,:));
title('Filtered ball trajectory');
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
fill3(T(:,1),T(:,2),T(:,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
s1.MarkerEdgeColor = s1.CData; % due to a bug in MATLAB R2015b
s2.MarkerEdgeColor = s2.CData;
s3.MarkerEdgeColor = s3.CData;
legend('cam1','cam3','filter');
hold off;

disp('Difference between actual bouncing pt and predicted state');
err(1) = bBounce(2) - bBouncePredict(1);
err(2) = bBounce(3) - bBouncePredict(2);
err(3) = bBounce(4) - bBouncePredict(3);
fprintf('e1 = %f, e2 = %f, e3 = %f.\n',err(1),err(2),err(3));
fprintf('Norm error = %f.\n',norm(err));

disp('Difference between actual hitting pt and predicted state');
% err(1) = bx1(idx) - ballEKF(1,end);
% err(2) = by1(idx) - ballEKF(2,end);
% err(3) = bz1(idx) - ballEKF(3,end);
err(1) = bx1(idx) - ballNLS(1,end);
err(2) = by1(idx) - ballNLS(2,end);
err(3) = bz1(idx) - ballNLS(3,end);
fprintf('e1 = %f, e2 = %f, e3 = %f.\n',err(1),err(2),err(3));
fprintf('Norm error = %f.\n',norm(err));


% x-y-z vs t
% tAfterNet = tTillNet(end) + dt * (1:N_pred);
% tTotal = [tTillNet; tAfterNet(:)];
% figure(4);
% subplot(3,1,1);
% plot(tTotal, ballEKF(1,:));
% title('xaxis');
% xlabel('t');
% ylabel('x');
% subplot(3,1,2);
% plot(tTotal, ballEKF(2,:));
% title('yaxis');
% xlabel('t');
% ylabel('y');
% subplot(3,1,3);
% plot(tTotal, ballEKF(3,:));
% title('zaxis');
% xlabel('t');
% ylabel('z');