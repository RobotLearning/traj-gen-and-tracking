%% Test spin model

clc; clear; close all; rng(2);
% load table parameters
loadTennisTableValues;
load('ballInitDist1.mat','mu','Sigma');

params.Clift = Clift;
params.C = Cdrag;
params.g = gravity;
params.zTable = table_z;
params.yNet = dist_to_table - table_y;
params.table_length = table_length;
params.table_width = table_width;
params.radius = ball_radius;
% coeff of friction vector
params.mu = 0.15;
params.CRT = CRT;
params.ALG = 'RK4';

ball_func = @(x,u,dt) discreteBallSpinModel(x,dt,params);

% initialize ball with high topspin (3000rpm or more)
w0 = [-50*2*pi;0;0]; %3000rpm topspin
phi0 = rand(3,1);
pv0 = mu + chol(Sigma)*randn(6,1);
x0 = [pv0(1:3);phi0(:);pv0(4:6);w0(:)];

totalTime = 0.8;
dt = 0.01;
N = totalTime/dt;
t = dt:dt:N*dt;
x = zeros(12,N);
x(:,1) = x0;

for i = 1:N-1
    x(:,i+1) = ball_func(x(:,i),0,dt);
end
balls = x(1:3,:)';

%% Try the spin-free model on same initial state

params.CFTX = CFTX;
params.CFTY = CFTY;
params.CRT = CRT;
ball_func_nospin = @(x,u,dt) discreteBallFlightModel(x,dt,params);

x2 = zeros(6,N);
x2(:,1) = [x0(1:3);x0(7:9)];
for i = 1:N-1
    x2(:,i+1) = ball_func_nospin(x2(:,i),0,dt);
end

%% Plot the results

figure;
s1 = scatter3(x(1,:),x(2,:),x(3,:),'b');
hold on;
s2 = scatter3(x2(1,:),x2(2,:),x2(3,:),'r');
title('Predicted ball trajectory');
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
fill3(T(1:4,1),T(1:4,2),T(1:4,3),[0 0.7 0.3]);
fill3(net(:,1),net(:,2),net(:,3),[0 0 0]);
s1.MarkerEdgeColor = s1.CData;
s2.MarkerEdgeColor = s2.CData;
legend('spinning ball','no spin');
hold off;

%% Prediction as we get more ball data with spin-free model

% initialize EKF
spin.flag = false;
spin.Clift = Clift;
spin.est = w0;
filter = initFilterEKF(spin);

% filter balls over the trajectory balls data
ballEsts = filterBallsEKF(t,balls,filter,spin);

idx_start = 13;
idx_end = 60;
rms_pred = calcModelPredErrors(filter,ballEsts,idx_start,idx_end,t,balls,true);      
plotPredErrors(rms_pred);