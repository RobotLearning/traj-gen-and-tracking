%% Analyzing ball data and robot joint data from demonstrations

clc; clear; close all

%% Initialize Barrett WAM

initializeWAM;

%% Load table values, ball and robot data

set = 15:65; % dataset of demonstrations
dropSet = [17,20,22,23,24,29,56]; % bad recordings
for i = 1:length(dropSet)
    set = set(set ~= dropSet(i));
end

dof = 7; % seven degrees of freedom 
scale = 1e-3; % recorded in milliseconds
Q = cell(length(set),1); % concatenate to huge q and qd matrices
Qd = cell(length(set),1);
t_strike = cell(length(set),1); % different time profiles

% load the data
i = 20;
    
Ms = dlmread(['../Desktop/okanKinestheticTeachin_20141210/unifyData', ...
             int2str(set(i)),'.txt']);

%% Load table values

% load table parameters
loadTennisTableValues;

table_z = floor_level - table_height;
table_x = table_center + table_width/2;
table_y = table_length/2;

T1 = [table_center - table_x; 
    dist_to_table + table_length; 
    table_z];
T2 = [table_center + table_x;
    dist_to_table + table_length;
    table_z];
T3 = [table_center + table_x;
    dist_to_table;
    table_z];
T4 = [table_center - table_x;
    dist_to_table;
    table_z];
T = [T1,T2,T3,T4]';

net1 = [table_center - table_x;
        dist_to_table + table_y;
        table_z];
net2 = [table_center + table_x;
        dist_to_table + table_y;
        table_z];
net3 = [table_center + table_x;
        dist_to_table + table_y;
        table_z + net_height];
net4 = [table_center - table_x;
        dist_to_table + table_y;
        table_z + net_height];

net = [net1,net2,net3,net4]';

%% Extract joint positions and reliable ball positions

% extract time and the joints - last 14 variables
Mq = Ms(:,end-14:end);

t = Mq(:,1);
t = scale * t;
q = Mq(:,2:dof+1);
qd = Mq(:,dof+2:end);
Q = [q';qd'];

zMax = 0.5;
% camera 1 time indices for reliable observations
idxCam1 = find(Ms(:,2) == 1 & Ms(:,5) >= table_z & Ms(:,5) < zMax);
tCam1 = t(idxCam1);
Mb1 = Ms(:,3:5);
bx1 = Mb1(idxCam1,1);
by1 = Mb1(idxCam1,2);
bz1 = Mb1(idxCam1,3);

% camera 3 time indices for reliable observations
idxCam3 = find(Ms(:,7) == 1 & Ms(:,10) >= table_z & Ms(:,10) < zMax);
tCam3 = t(idxCam3);
Mb3 = Ms(:,8:10);
bx3 = Mb3(idxCam3,1);
by3 = Mb3(idxCam3,2);
bz3 = Mb3(idxCam3,3);

% draw ball every Xth estimate
drawBallIter = 20;
tCam1Draw = tCam1(1:drawBallIter:end);
tCam3Draw = tCam3(1:drawBallIter:end);
precision = 4;
tCam1Cell = num2cell(tCam1Draw,precision);
tCam3Cell = num2cell(tCam3Draw,precision);
for i = 1:length(tCam1Cell)
    tCam1Cell{i} = num2str(tCam1Cell{i});
end
for i = 1:length(tCam3Cell)
    tCam3Cell{i} = num2str(tCam3Cell{i});
end
% annotate some of the ball positions
bx1Draw = bx1(1:drawBallIter:end);
by1Draw = by1(1:drawBallIter:end);
bz1Draw = bz1(1:drawBallIter:end);
bx3Draw = bx3(1:drawBallIter:end);
by3Draw = by3(1:drawBallIter:end);
bz3Draw = bz3(1:drawBallIter:end);

%% Plot cartesian coordinates of robot trajectory

x = wam.kinematics(Q);

for j = 1:3
    joints{j} = ['cart\_', int2str(j)];
    vel{j} = ['cart\_vel\_', int2str(j)];
    acc{j} = ['cart\_acc\_', int2str(j)];
    %err{i} = ['err_j_', int2str(i)];
end

figure;
for j = 1:3
    s(j) = subplot(3,1,j);
    plot(t,x(j,:));
    legend(joints{j});
end
title(s(1),'Robot cartesian trajectory x-y-z vs. t');

%% Estimate striking time

trajBall = [tCam1,bx1,by1,bz1;
            tCam3,bx3,by3,bz3];

% find closest points in space as a first estimate
distBallRobot = [sqrt((bx1-x(1,idxCam1)').^2 + (by1-x(2,idxCam1)').^2 + (bz1-x(3,idxCam1)').^2);
                 sqrt((bx3-x(1,idxCam3)').^2 + (by3-x(2,idxCam3)').^2 + (bz3-x(3,idxCam3)').^2)];
[minDist,idx] = min(distBallRobot);
tStrike = trajBall(idx,1);
disp('Rough estimate of striking time:');
disp(tStrike);

Ttot = 3.0; % total duration in sec
idxT = find(t >= (tStrike - Ttot/2) & t <= (tStrike + Ttot/2));
tCut = t(idxT);
trajRobot = [t,x(1,:)',x(2,:)',x(3,:)'];
% draw joint every Xth estimate
drawRobotIter = 40;

tRobotDraw = tCut(1:drawRobotIter:end);
precision = 4;
tRobotCell = num2cell(tRobotDraw,precision);
for i = 1:length(tRobotCell)
    tRobotCell{i} = num2str(tRobotCell{i});
end
% annotate some of the ball positions
rxDraw = x(1,idxT(1:drawRobotIter:end));
ryDraw = x(2,idxT(1:drawRobotIter:end));
rzDraw = x(3,idxT(1:drawRobotIter:end));


%% Plot cartesian coordinates of robot ball interaction

figure; 
scatter3(bx1,by1,bz1,'r');
hold on;
text(bx1Draw,by1Draw,bz1Draw,tCam1Cell);
scatter3(bx3,by3,bz3,'b');
text(bx3Draw,by3Draw,bz3Draw,tCam3Cell);
scatter3(x(1,idxT),x(2,idxT),x(3,idxT),'k');
text(rxDraw,ryDraw,rzDraw,tRobotCell);
title('Robot and Ball cartesian trajectories');
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
legend('cam1 ball','cam3 ball','robot');
fill3(T(:,1),T(:,2),T(:,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
hold off;