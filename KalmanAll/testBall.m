% Make a ball move in the 2D plane
% State = (x y xdot ydot). We only observe (x y).

% This code was used to generate Figure 15.9 of "Artificial Intelligence: a Modern Approach",
% Russell and Norvig, 2nd edition, Prentice Hall, 2003.

% x(t+1) = A x(t) + noise(Q)
% y(t) = C x(t) + noise(R)

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

[xfilt, Vfilt, VVfilt, loglik] = kalman_filter(y, A, C, Q, R, x0, V0);
[xsmooth, Vsmooth] = kalman_smoother(y, A, C, Q, R, x0, V0);

[xMyFilt, ~, ~, ~] = myKalman(y, A, C, Q, R, x0, V0);
[xMySmooth, ~] = myKalmanSmoother(y, A, C, Q, R, x0, V0);


dobs = y([1 2],:) - x([1 2],:);
mse_obs = sqrt(sum(sum(dobs.^2)))

dfilt = x([1 2],:) - xfilt([1 2],:);
mse_filt = sqrt(sum(sum(dfilt.^2)))

dMyfilt = x([1 2],:) - xMyFilt([1 2],:);
mse_myfilt = sqrt(sum(sum(dMyfilt.^2)))

dsmooth = x([1 2],:) - xsmooth([1 2],:);
mse_smooth = sqrt(sum(sum(dsmooth.^2)))

dMysmooth = x([1 2],:) - xMySmooth([1 2],:);
mse_mysmooth = sqrt(sum(sum(dMysmooth.^2)))


figure(1)
plot(x(1,:), x(2,:), 'ks-', y(1,:), y(2,:), 'g*', ...
     xfilt(1,:), xfilt(2,:), 'rx:', xsmooth(1,:), xsmooth(2,:), 'b*');
xlabel('x');
ylabel('y');

end

function [x, V] = myKalmanSmoother(y,A,C,Q,R,x0,V0)

[x,V,xpred,Vpred] = myKalman(y,A,C,Q,R,x0,V0);

N = size(y,2);
for i = N-1:-1:1
    % Rauch recursion 
    H = V(:,:,i) * A' * inv(Vpred(:,:,i));
    x(:,i) = x(:,i) + H * (x(:,i+1) - xpred(:,i));
    V(:,:,i) = V(:,:,i) + H * (V(:,:,i+1) - Vpred(:,:,i)) * H';
end


end

function [x, V, xpred, Vpred] = myKalman(y,A,C,Q,R,x0,V0)

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
    V(:,:,i) = filter.P;
    filter.predict(0);
    xpred(:,i) = filter.x;
    Vpred(:,:,i) = filter.P;
end
filter.update(y(:,N),0);
x(:,N) = filter.x;
V(:,:,N) = filter.P;

end