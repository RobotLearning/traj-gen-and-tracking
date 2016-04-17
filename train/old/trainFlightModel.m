%% Train flight model

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
dropSet1 = [17,20,22,23,24,27,29,...
            31,32,33,34,35,38,41,48,49,53,54,56]; % bad recordings
set = setdiff(set1,dropSet1);
scale  = 0.001; % recorded in milliseconds

% use camera 1? 1 might be unreliable
useCam1 = false;

for idx = 1:length(set)

    %% Get rid of outliers
    M = dlmread([demoFolder1,int2str(set(idx)),'.txt']);
    t = scale * M(:,end-14);

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
    
    if useCam1
        b = [tCam1,b1;
             tCam3,b3];
    else
        b = [tCam3,b3]; % camera3 only
    end
    b = sortrows(b);
    
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
    % make sure time starts from zero
    b(:,1) = b(:,1) - b(1,1);
    
    %% Get ball positions till bounce

    % find the last ball position before first bounce
    % detect first bounce
        
    tol = 20e-2;    
    idxBallBounce = find(b(:,4) <= table_z + ball_radius + tol);
    idxJump = find(diff(idxBallBounce) ~= 1);
    if isempty(idxJump)
        [~,idxBounceCandidate] = min(b(idxBallBounce,4));
        idxBallBounce = idxBallBounce(idxBounceCandidate)-1;
    else
        idxBallBounce = idxBallBounce(1:idxJump(1));
        [~,idxBounceCandidate] = min(b(idxBallBounce,4));
        idxBallBounce = idxBallBounce(idxBounceCandidate);
    end
    bTillTable = b(1:idxBallBounce,1:4);
    ballData = bTillTable;
    
    %% Provide good initial ball estimates
    
    % using polyfit on first 12 balls
    sampleSize = 12; 
    M = [ones(sampleSize,1),ballData(1:sampleSize,1),ballData(1:sampleSize,1).^2];
    Y = ballData(1:sampleSize,2:4);
    beta = M \ Y;
    ballInitPosEst = beta(1,:);
    ballInitVelEst = beta(2,:);
    ballInit = [ballInitPosEst,ballInitVelEst];
    
    % Run nonlinear least squares to estimate ballInit better
    x0 = ballInit(:);
    ballFun = @(b0,C,g) predictNextBall(b0,C,g,ballData,length(ballData));
    fnc = @(x) ballFun(x(1:6),Cdrag,gravity);
    options = optimoptions('lsqnonlin');
    options.Display = 'final';
    options.Algorithm = 'levenberg-marquardt';
    options.MaxFunEvals = 1500;
    [x,err] = lsqnonlin(fnc,x0,[],[],options);
    ballInit = x(1:6);
    
    % data used for optimization is saved here
    ballTrain.bInit(:,idx) = ballInit(:);
    ballTrain.data{idx} = ballData;
    ballTrain.error{idx} = err;
end

save('train/ballTrain.mat','ballTrain');

%% Nonlinear least squares to estimate params

%%{
b0 = ballTrain.bInit(:);
x0 = [b0; Cdrag; gravity];
data = [];
for i = 1:length(set)
    data = [data;ballTrain.data{idx}];
    lenSet(i) = size(ballTrain.data{idx},1);
end
ballFun = @(b0,C,g) predictNextBall(b0,C,g,data,lenSet);
fun = @(x) ballFun(x(1:end-2),x(end-1),x(end));
options = optimoptions('lsqnonlin','Display','iter','MaxFunEvals',3000);
x = lsqnonlin(fun,x0,[],[],options);
C = x(end-1)
g = x(end)
%}


%% Test fitted trajectories

idx = 10;
bInit = ballTrain.bInit;
% C = Cdrag;
% g = gravity;
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
fill3(T(1:4,1),T(1:4,2),T(1:4,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
hold off;