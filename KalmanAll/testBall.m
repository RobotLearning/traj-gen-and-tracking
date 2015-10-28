% Make a ball move in the 2D plane
% State = (x y xdot ydot). We only observe (x y).

% This code was used to generate Figure 15.9 of "Artificial Intelligence: a Modern Approach",
% Russell and Norvig, 2nd edition, Prentice Hall, 2003.

% X(t+1) = F X(t) + noise(Q)
% Y(t) = H X(t) + noise(R)

function testBall

clc; clear; close all

numx = 4; % state size
numy = 2; % observation size
A = [1 0 1 0; 
     0 1 0 1; 
     0 0 1 0; 
     0 0 0 1]; 
C = [1 0 0 0; 
     0 1 0 0];
Q = 0.1 * eye(numx);
R = 1 * eye(numy);
x0 = [10 10 1 0]';
V0 = 10 * eye(numx);

seed = 9;
rng(seed);
T = 15;
x(:,1) = x0;
y(:,1) = C * x(:,1) + chol(R) * randn(numy, 1);

numu = 1;
B = zeros(numx, numu);
u = zeros(numu, T-1);

for t = 2:T
  x(:,t) = A*x(:,t-1) + B*u(:,t-1) + chol(Q)*randn(numx,1);
  y(:,t) = C*x(:,t)  + chol(R)*randn(numy,1);
end

[xMyFilt, VMyFilt] = myKalman(y, A, C, Q, R, x0, V0);

[xfilt, Vfilt, VVfilt, loglik] = kalman_filter(y, A, C, Q, R, x0, V0);
[xsmooth, Vsmooth] = kalman_smoother(y, A, C, Q, R, x0, V0);

dobs = y([1 2],:) - x([1 2],:);
mse_obs = sqrt(sum(sum(dobs.^2)))

dfilt = x([1 2],:) - xfilt([1 2],:);
mse_filt = sqrt(sum(sum(dfilt.^2)))

dMyfilt = x([1 2],:) - xMyFilt([1 2],:);
mse_myfilt = sqrt(sum(sum(dMyfilt.^2)))

dsmooth = x([1 2],:) - xsmooth([1 2],:);
mse_smooth = sqrt(sum(sum(dsmooth.^2)))


figure(1)
plot(x(1,:), x(2,:), 'ks-', y(1,:), y(2,:), 'g*', ...
     xfilt(1,:), xfilt(2,:), 'rx:', xMyFilt(1,:), xMyFilt(2,:), 'b*');
xlabel('x');
ylabel('y');

end

function [x, V] = myKalman(y,A,C,Q,R,x0,V0)

mats.A = A;
mats.B = 0;
mats.C = C;
mats.D = 0;
mats.O = Q;
mats.M = R;

N = size(y,2);
filter = KF(length(x0),mats);
filter.initState(x0,V0);
for i = 1:N-1
    filter.update(y(:,i),0);
    x(:,i) = filter.x;
    filter.predict(0);
end
filter.update(y(:,N),0);
x(:,N) = filter.x;
V = [];

end