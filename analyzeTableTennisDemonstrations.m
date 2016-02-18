%% Analyzing ball data and robot joint data from demonstrations

clc; clear; close all

%% Initialize Barrett WAM

initializeWAM;

%% Load table values, ball and robot data

%demoFolder = '../Desktop/okanKinestheticTeachin_20141210/unifyData';
demoFolder = '../Desktop/okanKinestheticTeachin_20151124/unifyData';
%set = 15:65; % dataset of demonstrations
set = 2:64;
%dropSet = [17,20,22,23,24,29,56]; % bad recordings
dropSet = [5,6,9,10,14,26,27,32,38,42,55,56,57];
for i = 1:length(dropSet)
    set = set(set ~= dropSet(i));
end

dof = 7; % seven degrees of freedom 
scale = 1e-3; % recorded in milliseconds

% put everything in a structure
demoss = struct();

% load the data
for i = 1:length(set)
    demoss(i).name = ['demo ',int2str(set(i))];
    demoss(i).raw = dlmread([demoFolder,int2str(set(i)),'.txt']);
    % extract time and the joints - last 14 variables
    %robot = demo(i).raw(:,end-14:end);
    robot = demoss(i).raw(:,end-20:end);
    q = robot(:,2:dof+1);
    qd = robot(:,dof+2:dof+8);
    x_SL = robot(:,dof+9:dof+11);
    xd_SL = robot(:,dof+12:dof+14);
    demoss(i).t = scale * robot(:,1);
    demoss(i).Q = [q';qd'];
    demoss(i).x = [x_SL'; xd_SL'];
end

%% Load table values

% load table parameters
loadTennisTableValues;

%% Extract reliable ball positions and estimate striking time

for i = 1:length(set)

    M = demoss(i).raw;
    Q = demoss(i).Q;
    t = demoss(i).t;
    
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

    % Get cartesian coordinates of robot trajectory
    [x,xd,o] = wam.kinematics(Q);
    %{
    % Differentiate x
    xdiff = diff(x')';
    xdiff = [xdiff, xdiff(:,end)];    
    % Compare with x coming from SL
    figure(3);
    view(3);
    plot3(xdiff(1,:),xdiff(2,:),xdiff(3,:));
    hold on;
    plot3(demo(i).x(4,:),demo(i).x(5,:),demo(i).x(6,:));
    grid on;
    axis equal
    xlabel('x');
    ylabel('y');
    zlabel('z');
    legend('xd matlab','xd sl');
    hold off;
    figure(4)
    plot(t,xd(1,:),t,demo(i).x(4,:));
    figure(5)
    plot(t,xd(2,:),t,demo(i).x(5,:));
    figure(6)
    plot(t,xd(3,:),t,demo(i).x(6,:));
    %}

    camInfo = [tCam1,bx1,by1,bz1;
                tCam3,bx3,by3,bz3];

    % find closest points in space as a first estimate
    distBallRobot = [sqrt((bx1-x(1,idxCam1)').^2 + (by1-x(2,idxCam1)').^2 + (bz1-x(3,idxCam1)').^2);
                     sqrt((bx3-x(1,idxCam3)').^2 + (by3-x(2,idxCam3)').^2 + (bz3-x(3,idxCam3)').^2)];
    [minDist,idx] = min(distBallRobot);
    tStrike = camInfo(idx,1);
    fprintf('Striking time est. for demo %d: %f\n',set(i),tStrike);

    %% Segment the signals
    duration = 0.6; % total duration in sec
    strikeDuration = 0.5;
    idxD = find(t >= (tStrike - strikeDuration) & t <= (tStrike + duration - strikeDuration));
    tCutStrike = t(idxD);
    
    % for discrete DMPs
    demoss(i).t_strike = tCutStrike;
    demoss(i).o_strike = o(:,idxD);
    demoss(i).x_strike = x(:,idxD);
    demoss(i).Q_strike = Q(:,idxD);
    
    % for rhythmic DMPs we need to segment differently
    Treturn = t(end) - tStrike;
    idxStrike = find(t >= (tStrike - duration) & t <= tStrike);
    % do not take all of the return times, skip some to make it faster
    jump = floor(Treturn / duration);
    idxReturn = idxStrike(end)+jump:jump:length(t); 
    idxR = [idxStrike', idxReturn];
    idxStrikeReturn = [idxStrike',idxStrike(end)+1:idxStrike(end)+length(idxReturn)];
    tStrikeReturn = t(idxStrikeReturn);
    
    demoss(i).t_rdmp = tStrikeReturn;
    demoss(i).Q_rdmp = Q(:,idxR);

end

%% Plot cartesian coordinates of robot ball interaction for last data

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

% draw joint every Xth estimate
drawRobotIter = 40;

tRobotDraw = tCutStrike(1:drawRobotIter:end);
precision = 4;
tRobotCell = num2cell(tRobotDraw,precision);
for i = 1:length(tRobotCell)
    tRobotCell{i} = num2str(tRobotCell{i});
end
% annotate some of the robot positions
rxDraw = x(1,idxD(1:drawRobotIter:end));
ryDraw = x(2,idxD(1:drawRobotIter:end));
rzDraw = x(3,idxD(1:drawRobotIter:end));

for j = 1:3
    cartpos{j} = ['cart\_', int2str(j)];
    cartvel{j} = ['cart\_vel\_', int2str(j)];
    cartacc{j} = ['cart\_acc\_', int2str(j)];
    %err{i} = ['err_j_', int2str(i)];
end

figure(1);
for j = 1:3
    s(j) = subplot(3,1,j);
    plot(t,x(j,:));
    legend(cartpos{j});
end
title(s(1),'Robot cartesian trajectory x-y-z vs. t');

figure(2);
scatter3(bx1,by1,bz1,'r');
hold on;
text(bx1Draw,by1Draw,bz1Draw,tCam1Cell);
scatter3(bx3,by3,bz3,'b');
text(bx3Draw,by3Draw,bz3Draw,tCam3Cell);
scatter3(x(1,idxD),x(2,idxD),x(3,idxD),'k');
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