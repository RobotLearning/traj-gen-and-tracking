%% Test Extended Kalman filter

clc; clear; close all;
seed = 5;
rng(seed);
eps = 1e-4;
N = 20;
dt = 0.02;
x0 = rand(3,1);
xd0 = rand(3,1);

loadTennisTableValues;
table_z = floor_level - table_height;
table_x = table_center + table_width/2;
table_y = table_length/2;

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

%% Testing with the table tennis flight model with air drag

x(:,1) = [x0;xd0];
y(:,1) = x(1:3,1);
yNs(:,1) = y(:,1) + sqrt(eps) * randn(3,1);
% symplectic Euler
for i = 2:N
    x(:,i) = symplecticFlightModel(x(:,i-1),dt,params);
    y(:,i) = x(1:3,i);
    yNs(:,i) = y(:,i) + sqrt(eps) * randn(3,1);
end

% initialize EKF
dim = 6;
C = [eye(3),zeros(3)];
funState = @(x,u,h) symplecticFlightModel(x,h,params);
% very small but nonzero value for numerical stability
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);floor_level
filter = EKF(dim,funState,mats);
filter.initState(x(:,1),eps);
for i = 1:N-1
    filter.linearize(dt,0);
    filter.update(yNs(:,i),0);
    yEKF(:,i) = C * filter.x;
    filter.predict(dt,0);
end
filter.update(yNs(:,N),0);
yEKF(:,N) = C * filter.x;

% Extended Kalman Smoother
t = dt * (1:N);
u = zeros(1,N);
filter.initState(x(:,1),eps);
[xsmooth,~] = filter.smooth(t,yNs,u);
yEKF_smooth = C * xsmooth;

SSE_ekf = trace((yEKF - y)*(yEKF - y)')
SSE_ekf_smooth = trace((yEKF_smooth - y)*(yEKF_smooth - y)')
SSE_ns = trace((yNs - y)*(yNs - y)')

plot3(y(1,:), y(2,:), y(3,:), 'ks-', ...
      yNs(1,:), yNs(2,:), yNs(3,:), 'b*', ...
      yEKF(1,:), yEKF(2,:), yEKF(3,:), 'rx:', ...
      yEKF_smooth(1,:), yEKF_smooth(2,:), yEKF_smooth(3,:), '-rd');
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
axis tight;
legend('actual','observations','EKF', 'EKF smoother');