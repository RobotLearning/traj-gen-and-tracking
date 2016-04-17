%% AnalyzeBallData

% analyzes the following prediction scenarios:
% 
% prediction with Extended Kalman Filter with initial 
% ball position and velocity calculated with NLS / initial data pts.
%
% LLS-based polynomial filtering

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
% bad recordings
dropSet1 = [17,20,22,23,24,29,31,32,34,35,38,41,48,53,54,56]; 
dropSet2 = [5,6,9,10,14,26,27,32,38,42,55,56,57];
% 19 doesnt fit well
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

choose = 15; %set(end);
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
b1 = M(idxCam1,3:5);

% camera 3 time indices for reliable observations
idxCam3 = find(M(:,7) == 1 & M(:,10) >= table_z & M(:,10) < zMax);
tCam3 = t(idxCam3);
b3 = M(idxCam3,8:10);

% Estimate striking time
% find closest points in space as a first estimate

diffRobot2Ball = b1 - x(:,idxCam1)';
[minDist,idxStrike] = min(diag(diffRobot2Ball*diffRobot2Ball'));
tStrike = tCam1(idxStrike,1);
fprintf('Striking time est. for demo %d: %f\n',choose,tStrike);

%% initialize EKF

params.C = Cdrag;
params.g = gravity;
params.zTable = table_z;
params.yNet = dist_to_table - table_y;
params.table_length = table_length;
params.table_width = table_width;
params.radius = ball_radius;
% coeff of restitution-friction vector
params.CFTX = CFTX;
params.CFTY = CFTY;
params.CRT = CRT;
params.ALG = 'RK4';

funState = @(x,u,dt) discreteBallFlightModel(x,dt,params);
% very small but nonzero value for numerical stability
eps = 1e-6; dim = 6;
C = [eye(3),zeros(3)];
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);
filter = EKF(dim,funState,mats);

%% Order the ball data w.r.t. time and throw away same values

useCam1 = false;

if useCam1 || strcmp(demo,demoFolder1)
    b = [tCam1,b1; % camera 1
         tCam3,b3]; % camera 3
else
    b = [tCam3,b3];
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
tPredictTillLastBall = tStrike - bTillNet(end,1);

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
tBounce = bTillTable(end,1) - bTillNet(end,1);

% % make sure t starts from 0
bTillNet(:,1) = bTillNet(:,1) - bTillNet(1,1);

%% Estimate initial state with NLS

guessBallInitVel = 1.0 * [-0.80; 7.0; 2.0];
ballInitNLS = [bTillNet(1,2:4)';guessBallInitVel];
% Run nonlinear least squares to estimate ballInit better
%x0 = [ballInitNLS(:);Cdrag;gravity];
x0 = ballInitNLS(:);
ballData = bTillNet;
sampleSize = size(ballData,1);
ballFun = @(b0,C,g) predictNextBall(b0,C,g,ballData,sampleSize);
%fnc = @(x) ballFun(x(1:6),x(7),x(8));
fnc = @(x) ballFun(x(1:6),Cdrag,gravity);
tic;
x = lsqnonlin(fnc,x0);
toc
ballInitNLS = x(1:6);
% Cdrag_est = x(7)
% gravity_est = x(8)
p0 = ballInitNLS(1:3);
v0 = ballInitNLS(4:6);

%% LS prediction based on projectile motion on last 12 samples

%{
sampleSize = 12; %size(bTillNet,1);
M = [ones(sampleSize,1),...
    bTillNet(end-sampleSize+1:end,1),bTillNet(end-sampleSize+1:end,1).^2];
Y = bTillNet(end-sampleSize+1:end,2:4); %2:3
% Yz = bTillNet(:,4) - 0.5*gravity*bTillNet(:,1).^2;
% Y = [Y,Yz];
tol = 0.00;
betaPreBounce = M \ Y;
%ballInitLS = pinv(M,tol)*Y;

% find the time till bounce

m = betaPreBounce(:,3);
a = m(3);
b = m(2);
c = m(1) - table_z - ball_radius;
tBouncePredLS = (-b - sqrt(b^2 - 4*a*c))/(2*a);

dt = 1/60;
N_pred_LS = round(tBouncePredLS/dt);
tLS = dt * (1:N_pred_LS);
Mstar = [ones(N_pred_LS,1),tLS(:),tLS(:).^2];
ballLS = Mstar * betaPreBounce;% + [zeros(N_pred_LS,2),0.5*gravity*tLS(:).^2];
ballLLS = ballLS';

% predict till strike
M = diag([CFTX,CFTY,-CRT]);
% transform ballInitLS
mat = [1, tBouncePredLS;
       0, 1];
   
c = betaPreBounce(1,:);
b = betaPreBounce(2,:);
a = betaPreBounce(3,:);
y = [b*tBouncePredLS + c; (2*a*tBouncePredLS + b)*M - 2*a*tBouncePredLS];
betaAfterBounce = mat \ y;
betaAfterBounce = [betaAfterBounce; a];

N_pred_LS = round((tPredictLS(end) - tBouncePredLS)/dt);
tLS = tBouncePredLS + dt * (1:N_pred_LS);
MstarAfter = [ones(N_pred_LS,1),tLS(:),tLS(:).^2];
ballLSAfter = MstarAfter * betaAfterBounce;
ballLLS = [ballLLS,ballLSAfter'];
%}

%% Predict the ball with EKF till striking time

obsLen = size(bTillNet,1);
% initialize the state with first observation and first derivative est.
% p0 = bTillNet(1,2:4);
% v0 = (bTillNet(3,2:4) - bTillNet(1,2:4))/(bTillNet(3,1)-bTillNet(1,1));
initVar = 0.1;
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
dt = 1/60;
N_strike = round(tPredictTillLastBall/dt);
N_bounce = round(tBounce/dt);
% after bounce there should be less topspin effects
spinChange = true;

for i = 1:N_bounce
    %filter.linearize(dt,0);
    filter.predict(dt,0);
    ballEKF(:,obsLen+i) = C * filter.x;
    % try changing constants
    if filter.x(6) > 0 && spinChange
        gravity = -9.80;
        params.g = gravity;
        funAfterBounce = @(x,u,dt) discreteBallFlightModel(x,dt,params);
        filter.f = funAfterBounce;
        spinChange = false;
    end
        
end

% Smoothen the ball path EKF smoother to get x0
% 
% filter.initState([p0(:);v0(:)],eps);
% u = zeros(1,size(bTillNet,1));
% tTillNet = bTillNet(:,1);
% [xEKFSmooth, VekfSmooth] = filter.smooth(tTillNet,bTillNet(:,2:4)',u);
% ballEKFSmooth = C * xEKFSmooth;


%% Plot predictions

[~,idxMaxY] = max(b1(:,2));

figure;
s1 = scatter3(b1(1:idxStrike,1),b1(1:idxStrike,2),b1(1:idxStrike,3),'r');
hold on;
s2 = scatter3(b3(1:idxMaxY,1),b3(1:idxMaxY,2),b3(1:idxMaxY,3),'b');
s3 = scatter3(ballEKF(1,:),ballEKF(2,:),ballEKF(3,:));
%s4 = scatter3(ballLLS(1,:),ballLLS(2,:),ballLLS(3,:),'y');
title('Filtered ball trajectory');
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
fill3(T(1:4,1),T(1:4,2),T(1:4,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
s1.MarkerEdgeColor = s1.CData; % due to a bug in MATLAB R2015b
s2.MarkerEdgeColor = s2.CData;
s3.MarkerEdgeColor = s3.CData;
%s4.MarkerEdgeColor = s4.CData;
legend('cam1','cam3','EKF-NLS');
hold off;


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