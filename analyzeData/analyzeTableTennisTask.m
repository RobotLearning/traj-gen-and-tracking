%% Look at error in min_acc_task (optimal polynomials)

clc; clear; close all;
% addpath('../saveData/');

load('LookupTable-16-May-2016.mat');
% load('../Desktop/data/realBallData_030516.txt');
% load('../Desktop/data/realJointData_030516.txt');

date = '030516';
moment = '20:1:22';

% ball data
B = load(['../Dropbox/data/realBallData_', date, '.txt']);
% joint data
J = load(['../Dropbox/data/realJointData_', date, '.txt']);

t = J(:,1);
q_des = J(:,2:8);
qd_des = J(:,9:15);
q_act = J(:,16:22);
qd_act = J(:,23:end);

joints = cell(1,7);
cartpos = cell(1,7);
vels = cell(1,7);
accs = cell(1,7);
ff = cell(1,7);
fb = cell(1,7);
utot = cell(1,7);
cartvel = cell(1,7);
for i = 1:7
    cartpos{i} = strcat('x_',int2str(i));
    cartvel{i} = strcat('xd_',int2str(i));
    joints{i} = strcat('q_',int2str(i));
    vels{i} = strcat('qd_',int2str(i));
    accs{i} = strcat('qdd_',int2str(i));
    ff{i} = strcat('uff_',int2str(i));
    fb{i} = strcat('ufb_',int2str(i));
    utot{i} = strcat('utot_',int2str(i));
end

scale = 0.001; % recorded in milliseconds
t = scale * t;

%% load barrett wam and kinematics
% Initialize Barrett WAM
initializeWAM;
% load table parameters
loadTennisTableValues;
% using ball gun or throwing by hand
ballgun = true;

if ballgun
    CRT = 0.96;
    CFTY = 1.20;
    CFTX = 1.20;
    gravity = -11.06;
    Cdrag = 0.1753;
    % post bounce
    %Cdrag_post = 0.1968;
    %gravity_post = -10.83;
end


Q_des = [q_des';qd_des'];
Q_act = [q_act';qd_act'];

%% initialize EKF to predict ball path

params.C = Cdrag;
params.g = gravity;
params.zTable = table_z;
params.yNet = dist_to_table - table_y;
params.table_length = table_length;
params.table_width = table_width;
params.radius = ball_radius;
% coeff of restitution-friction vector
%params.BMAT = BMAT;
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

%% Get reliable ball observations
% this is to throw away really bad observations as to not confuse the
% filter
zMax = 0.5;

% camera 1 time indices for reliable observations
idx1 = find(B(:,2) == 1 & B(:,5) >= table_z & B(:,5) < zMax);
t1 = t(idx1);
% order the ball data w.r.t. time and throw away same values
b1 = [t1,B(idx1,3:5)];
b1 = sortrows(b1);
j = 1; % time from ballgun to net 
tol = 1e-3;
idxDiffBallPos1 = []; % index for diff ball positions
for i = 1:size(b1,1)-1
    if norm(b1(i,2:4) - b1(i+1,2:4)) > tol
        idxDiffBallPos1(j) = i;
        j = j+1;
    end
end
b1 = b1(idxDiffBallPos1,:);

% camera 3 time indices for reliable observations
idx3 = find(B(:,7) == 1 & B(:,10) >= table_z & B(:,10) < zMax);
t3 = t(idx3);
% order the ball data w.r.t. time and throw away same values
b3 = [t3,B(idx3,8:10)];
b3 = sortrows(b3);
j = 1; % time from ballgun to net 
tol = 1e-3;
idxDiffBallPos3 = []; % index for diff ball positions
for i = 1:size(b3,1)-1
    if norm(b3(i,2:4) - b3(i+1,2:4)) > tol
        idxDiffBallPos3(j) = i;
        j = j+1;
    end
end
b3 = b3(idxDiffBallPos3,:);

%% Find the trajectory corresponding to trial

trial = 2;
% note: we can also check for zero desired velocity

N = size(q_des,1);
q0 = q_des(1,:);
tol = 0;
diff2q0 = zeros(N,1);
for i = 1:N
    diff2q0(i) = (q_des(i,:) - q0)*(q_des(i,:) - q0)';
end

ind = find(diff2q0 == 0);
startTrjIdx = find(diff(ind) > 1);
startTrjIdx = [startTrjIdx'; startTrjIdx'+1];
startTrjIdx = startTrjIdx(:);
startTrjIdx = ind(startTrjIdx);

numTrials = length(startTrjIdx)/2;

idxPlot = startTrjIdx(2*trial-1):startTrjIdx(2*trial);

Q_des_plot = Q_des(:,idxPlot);
Q_act_plot = Q_act(:,idxPlot);
q_act_plot = Q_act_plot(1:7,:);
qd_act_plot = Q_act_plot(8:end,:);
q_des_plot = Q_des_plot(1:7,:);
qd_des_plot = Q_des_plot(8:end,:);
t_plot = t(idxPlot);

[x_des,xd_des,o_des] = wam.kinematics(Q_des_plot);
[x_act,xd_act,o_act] = wam.kinematics(Q_act_plot);

%% Get the number of incoming balls 

if trial >= numTrials
    warning('choose trial smaller than max');
    tFinish = ceil(t(startTrjIdx(end)));
else
    tFinish = t(startTrjIdx(2*trial+1));
end

% get the time of first ball 
time2Arrive2Net = 0.5;
tStart = t_plot(1) - time2Arrive2Net;
idxStart3 = b3(:,1) >= tStart & b3(:,1) <= tFinish;
idxStart1 = b1(:,1) >= tStart & b1(:,1) <= tFinish;

%% Predict ball using estimate (SL uses it for lookup table)

b3est = [t3,B(idx3,12:17)];
b3est = sortrows(b3est);
b3est = b3est(idxDiffBallPos3,:);

ballEst = b3est(idxStart3,:);
% remove first 12 entries
ballEst = ballEst(12+1:end,:);

toly = 0.4;

% particular rule implemented in lookup table in SL
y_center = dist_to_table - table_length/2;
idxPred = find(ballEst(:,3) > y_center & ballEst(:,3) < y_center + toly & ...
               ballEst(:,6) > 0.5,1);

ballLookup = ballEst(idxPred,2:end);
tLookup = ballEst(idxPred,1);

% get the indices for plotting
b3_plot = b3(idxStart3,2:4);
t3_plot = b3(idxStart3,1);
b1_plot = b1(idxStart1,2:4);
t1_plot = b1(idxStart1,1);

% remove outliers in b1plot after getting ball closest to robot
[~,idxClosest2Robot] = max(b1_plot(:,2));
b1_plot = b1_plot(1:idxClosest2Robot,:);
t1_plot = t1_plot(1:idxClosest2Robot);

% use ransac to further prune outliers
% outlier detection again
outlierIdx = detectOutlierBalls(t1_plot,b1_plot,1);
inlierIdx = setdiff(1:length(t1_plot),outlierIdx);
b1_plot = b1_plot(inlierIdx,:);
t1_plot = t1_plot(inlierIdx,:);

% find the index at bounce
[~,idxBounce] = min(b3_plot(:,3));
dtPredTillBounce =  t3_plot(idxBounce) - tLookup;
dtPredTillLastBlob3 = t3_plot(end) - tLookup;
dtPredTillLastBlob1 = t1_plot(end) - tLookup;

% predict using ball estimate
dt = 1/60;
initVar = 1;
filter.initState(ballLookup(:),initVar);
filter.linearize(dt,0);
predictHorizon = dtPredTillLastBlob1; % only predict till last blob3
table.DIST = dist_to_table;
table.LENGTH = table_length;
table.Z = table_z;
table.WIDTH = table_width;
[ballPred,ballTime,numBounce,time2PassTable] = ...
            predictBallPath(dt,predictHorizon,filter,table);
        
%% generate desired traj using lookup table and compare

dof = 7;
numLookup = size(X,1);
dist2Lookup = zeros(numLookup,1);
for i = 1:numLookup
    vec = X(i,:) - ballLookup;
    dist2Lookup(i) = vec*vec';
end
[~,idxLookup] = min(dist2Lookup);
paramLookup = Y(idxLookup,:);
qf = paramLookup(1:dof);
qfdot = paramLookup(dof+1:14);
time2hit = paramLookup(end);
time2return = 1.0;
q0dot = zeros(dof,1);
Q0 = [q0(:);q0dot(:)];  
Qf = [qf(:);qfdot(:)];

dt = 0.002;
idxHit = floor(time2hit / dt);


%{
% GET 3RD DEGREE POLYNOMIALS      
pStrike = generatePoly3rd(Q0,Qf,dt,time2hit);
qStrike = pStrike(1:dof,:);
qdStrike = pStrike(dof+1:2*dof,:);
qddStrike = pStrike(2*dof+1:end,:);

pReturn = generatePoly3rd(Qf,Q0,dt,time2return);
qReturn = pReturn(1:dof,:);
qdReturn = pReturn(dof+1:2*dof,:);
qddReturn = pReturn(2*dof+1:end,:);

q_des_lookup = [qStrike,qReturn];
qd_des_lookup = [qdStrike,qdReturn];
qdd_des_lookup = [qddStrike,qddReturn];   

figure('Name','Lookup table vs desired saved values');
for i = 1:7
    subplot(7,2,2*i-1)
    plot(t_plot(1:end-1),q_des_lookup(i,:),t_plot,q_des_plot(i,:));
    legend([joints{i},'des_lookup'],[joints{i},'des']);
    subplot(7,2,2*i)
    plot(t_plot(1:end-1),qd_des_lookup(i,:),t_plot,qd_des_plot(i,:));
    legend([vels{i},'des_lookup'],[vels{i},'des']);
end
%}

%% plot actual and desired values
figure('Name','Actual and Desired Joint pos and vel');
for i = 1:7
    subplot(7,2,2*i-1)
    plot(t_plot,q_act_plot(i,:),t_plot,q_des_plot(i,:));
    legend([joints{i},'act'],[joints{i},'des']);
    subplot(7,2,2*i)
    plot(t_plot,qd_act_plot(i,:),t_plot,qd_des_plot(i,:));
    legend([vels{i},'act'],[vels{i},'des']);
end

%% plot cartesian values

plotIdx = 5;
plotIdxs = 1:plotIdx:size(x_des,2);

figure('Name','Actual vs Desired Cartesian pos');
scatter3(x_des(1,plotIdxs), x_des(2,plotIdxs), x_des(3,plotIdxs),'r');
hold on;
scatter3(x_act(1,plotIdxs), x_act(2,plotIdxs), x_act(3,plotIdxs),'b');
scatter3(x_des(1,idxHit),x_des(2,idxHit),x_des(3,idxHit),200,'filled');
grid on
axis square
legend('des','act','des-hit');
hold off
xlabel('x')
ylabel('y')
zlabel('z')

%% plot cartesian values along with ball observations and table

plotIdx = 10;
plotIdxs = 1:plotIdx:size(x_des,2);
gray = [0.5020 0.5020 0.5020];
red = [1.0000 0.2500 0.2500];
figure('Name','Table tennis performance');
hold on;

s1 = scatter3(b1_plot(:,1),b1_plot(:,2),b1_plot(:,3),'r');
s3 = scatter3(b3_plot(:,1),b3_plot(:,2),b3_plot(:,3),'b');
%predColor = [0.200 0.200 0.200];
sP = scatter3(ballPred(1,:),ballPred(2,:),ballPred(3,:),'k');
scatter3(x_des(1,plotIdxs), x_des(2,plotIdxs), x_des(3,plotIdxs));
scatter3(x_act(1,plotIdxs), x_act(2,plotIdxs), x_act(3,plotIdxs));
scatter3(x_des(1,idxHit),x_des(2,idxHit),x_des(3,idxHit),100,'filled');

grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
legend('cam1','cam3','filter','robot des', 'robot act', 'des-hit');

% draw table and robot 
fill3(T(1:4,1),T(1:4,2),T(1:4,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
% get joints, endeffector pos and orientation for q0
[joint,ee,racket] = wam.drawPosture(q0);
endeff = [joint(end,:); ee];
plot3(joint(:,1),joint(:,2),joint(:,3),'k','LineWidth',10);
plot3(endeff(:,1),endeff(:,2),endeff(:,3),'Color',gray,'LineWidth',5);
fill3(racket(1,:), racket(2,:), racket(3,:),red);
hold off;

%% Give info about cartesian tracking error at hitting time
disp('========TRACKNG ERROR AT DES HITTING TIME=======')
vectoridx = 1:3;
fprintf('xdes[%d] = %f\n', [vectoridx; x_des(:,idxHit)']);
fprintf('xact[%d] = %f\n', [vectoridx; x_act(:,idxHit)']);
fprintf('Diff in cm = %f\n', norm(x_des(:,idxHit) - x_act(:,idxHit)));
fprintf('vdes[%d] = %f\n', [vectoridx; xd_des(:,idxHit)']);
fprintf('vact[%d] = %f\n', [vectoridx; xd_act(:,idxHit)']);
fprintf('Diff in cm = %f\n', norm(xd_des(:,idxHit) - xd_act(:,idxHit)));

R_des = quat2Rot(o_des);
R_act = quat2Rot(o_act);
fprintf('ndes[%d] = %f\n', [vectoridx; R_des(:,3,idxHit)']);
fprintf('nact[%d] = %f\n', [vectoridx; R_act(:,3,idxHit)']);
fprintf('Diff in cm = %f\n', norm(R_des(:,3,idxHit) - R_act(:,3,idxHit)));