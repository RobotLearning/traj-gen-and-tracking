%% Train bounce model

clc; clear; close all;

% Detect when it bounces
% Use all data before bounce
% Use nonlinear least squares to optimize flight model parameters
% C and g and also estimate initial ball pos and velocities

%% Load table values

% load table parameters
loadTennisTableValues;

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
set = setdiff(set2,dropSet2);
scale  = 0.001; % recorded in milliseconds

% use camera 3? 1 might be unreliable
useCam3 = true;
% use first dataset
demoFolder = demoFolder2;

for idx = 1:length(set)

    %% Get rid of outliers
    M = dlmread([demoFolder,int2str(set(idx)),'.txt']);
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
    
    if useCam3
        b1 = [tCam1,bx1,by1,bz1,ones(length(tCam1),1)];
        b3 = [tCam3,bx3,by3,bz3,ones(length(tCam3),1)]; % camera 3
        b = [b1;b3];
    else
        b = [tCam1,bx1,by1,bz1,ones(length(tCam1),1)]; % camera 1
    end
    [b,idxB] = sortrows(b);
    
    % make sure time starts from zero
    b(:,1) = b(:,1) - b(1,1);
    
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
    ballData = bTillTable;
    
    %% Provide good initial ball estimates

    % Fitting a second order polynomial doesnt work well
    %%{
    len = size(ballData,1);
    M = [ones(len,1),ballData(:,1)];
    Y = ballData(:,2:3);
    Yz = ballData(:,4) - 0.5*gravity*ballData(:,1).^2;
    Y = [Y,Yz];
    ballInit = M \ Y;
    ballInit = [ballInit(1,:),ballInit(2,:)];
    %}
    
    % Run nonlinear least squares to estimate ballInit better
    x0 = [ballInit(:); Cdrag; gravity];
    ballFun = @(b0,C,g) predictNextBall(b0,C,g,ballData,len);
    fnc = @(x) ballFun(x(1:6),x(7),x(8));
    x = lsqnonlin(fnc,x0);
    ballInit = x(1:6);
    
    % data used for optimization is saved here
    ballTrain.bInit(:,idx) = ballInit(:);
    ballTrain.data{idx} = ballData;
end

save('train/ballTrain2.mat','ballTrain');

%% Nonlinear least squares to estimate params

%%{
x0 = [ballTrain.bInit(:); Cdrag; gravity];
data = [];
for i = 1:length(set)
    data = [data;ballTrain.data{idx}];
    lenSet(i) = size(ballTrain.data{idx},1);
end
ballFun = @(b0,C,g) predictNextBall(b0,C,g,data,lenSet);
fun = @(x) ballFun(x(1:6*length(set)),x(end-1),x(end));
options = optimoptions('lsqnonlin','Display','iter','MaxFunEvals',3000);
x = lsqnonlin(fun,x0,[],[],options);
bInit = reshape(x(1:end-2),6,length(set));
C = x(end-1)
g = x(end)
%}


%% Test fitted trajectories

idx = 10;
bInit = ballTrain.bInit;
C = Cdrag;
g = gravity;
data = ballTrain.data{idx};

figure;
b0 = bInit(:,idx);
tIdx = data(:,1);
bIdx = data(:,2:4);
bFitIdx = fitBallPath(b0,C,g,tIdx);

scatter3(bIdx(:,1),bIdx(:,2),bIdx(:,3),'r');
hold on;
scatter3(bFitIdx(1,:),bFitIdx(2,:),bFitIdx(3,:),'k');
legend('obs','fit');
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
fill3(T(:,1),T(:,2),T(:,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
hold off;