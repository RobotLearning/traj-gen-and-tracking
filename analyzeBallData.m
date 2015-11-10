%% Using Extended Kalman Filter to predict ball path after net

clc; clear; close all;

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

%% Load last data

demoFolder = '../Desktop/okanKinestheticTeachin_20141210/unifyData';
set = 15:65; % dataset of demonstrations
dropSet = [17,20,22,23,24,29,56]; % bad recordings
for i = 1:length(dropSet)
    set = set(set ~= dropSet(i));
end
choose = 65; %set(end);
M = dlmread([demoFolder,int2str(choose),'.txt']);
scale = 0.001; % recorded in milliseconds
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
% coeff of restitution-friction vector
params.CFTX = CFTX;
params.CFTY = CFTY;
params.CRT = CRT;

funState = @(x,u,dt) symplecticFlightModel(x,dt,params);
% very small but nonzero value for numerical stability
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);
filter = EKF(dim,funState,mats);

%% Order the ball data w.r.t. time and throw away const values

b = [tCam1,bx1,by1,bz1,ones(length(tCam1),1); % camera 1
     tCam3,bx3,by3,bz3,3*ones(length(tCam3),1)]; % camera 3
b = sortrows(b);
% get the ones till net
% idxNet = 1;
% while b(idxNet,3) <= dist_to_table - table_y
%     idxNet = idxNet + 1;
% end
idxBallTillNet = find(b(:,3) <= dist_to_table - table_y);
N = 1; % time from ballgun to net 
i = 1; % index for diff ball positions
tol = 1e-3;
while idxBallTillNet(N+1) - idxBallTillNet(N) == 1
    if norm(b(N,2:4) - b(N+1,2:4)) > tol
        idxDiffBallPos(i) = N;
        i = i+1;
    end
    N = N + 1;
end
idxBallTillNet = idxBallTillNet(idxDiffBallPos);
bTillNet = b(idxBallTillNet,:);

%% Predict the ball with EKF

obsLen = size(bTillNet,1);
% initialize the state with first observation and first derivative est.
p0 = bTillNet(1,2:4);
v0 = (bTillNet(5,2:4) - bTillNet(1,2:4))/(bTillNet(5,1)-bTillNet(1,1));
filter.initState([p0(:);v0(:)],eps);
for i = 1:obsLen-1
    dt = bTillNet(i+1,1) - bTillNet(i,1);
    filter.linearize(dt,0);
    filter.update(bTillNet(i,2:4)',0);
    ballEKF(:,i) = C * filter.x;
    filter.predict(dt,0);
end
filter.update(bTillNet(obsLen,2:4)',0);
ballEKF(:,obsLen) = C * filter.x;

% Now we start predicting ball till we hit the table
tPredict = 0.8;
dt = 0.01;
N_pred = tPredict/dt;
for i = 1:N_pred
    %filter.linearize(dt,0);
    filter.predict(dt,0);
    ballEKF(:,obsLen+i) = C * filter.x;
end

figure(3);
scatter3(bx1,by1,bz1,'r');
hold on;
scatter3(bx3,by3,bz3,'b');
scatter3(ballEKF(1,:),ballEKF(2,:),ballEKF(3,:));
title('Filtered ball trajectory');
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
legend('cam1','cam3','EKF');
fill3(T(:,1),T(:,2),T(:,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
hold off;