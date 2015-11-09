%% Test Extended Kalman filter

clc; clear; close all;
seed = 5;
rng(seed);
eps = 1e-6;
N = 20;
dt = 0.02;
x0 = rand(3,1);
xd0 = rand(3,1);
drag = 0.1414;
g = -9.802;

%% Testing with the table tennis flight model with air drag

x(:,1) = [x0;xd0];
y(:,1) = x(1:3,1);
yNs(:,1) = y(:,1) + sqrt(eps) * randn(3,1);
% symplectic Euler
for i = 2:N
    x(:,i) = symplecticFlightModel(x(:,i-1),dt,drag,g);
    y(:,i) = x(1:3,i);
    yNs(:,i) = y(:,i) + sqrt(eps) * randn(3,1);
end

% initialize EKF
dim = 6;
C = [eye(3),zeros(3)];
funState = @(x,u,h) symplecticFlightModel(x,h,drag,g);
% very small but nonzero value for numerical stability
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);
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

SSE(1) = trace((yEKF - y)*(yEKF - y)');
SSE(2) = trace((yNs - y)*(yNs - y)')

plot3(y(1,:), y(2,:), y(3,:), 'ks-', ...
      yNs(1,:), yNs(2,:), yNs(3,:), 'b*', ...
      yEKF(1,:), yEKF(2,:), yEKF(3,:), 'rx:');
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
axis tight;
legend('actual','observations','EKF');