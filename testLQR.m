%% Test script for LQR

clc; clear; close all;

% 3D Linear Time Invariant test for LQR
% Example taken from http://www.egr.msu.edu/classes/me851/jchoi/lecture/Lect_14.pdf

dimx = 3;
dimu = 1;

% simulation variables
t0 = 0;
tf = 4;
h = 0.02;
t = t0:h:tf;
N = length(t)-1;

% system and weighting matrices
Q = eye(dimx);
R = eye(dimu);
% continuous time matrices
A = [0 1 0; 0 0 1; -0.4 -4.2 -2.1];
B = [0;0;1];
% trick to get discrete time versions
Mat = [A, B; zeros(dimu, dimx + dimu)];
MD = expm(h * Mat);
Ad = MD(1:dimx,1:dimx);
Bd = MD(1:dimx,dimx+1:end);

% optimal feedback matrix
MODE.N = Inf;
MODE.LTI = true;
K = LQR(Q,R,Ad,Bd,MODE);

% simulate system
x = zeros(dimx,N+1);
u = zeros(dimu,N);
x0 = ones(dimx,1);
x(:,1) = x0;

for i = 1:N
    if length(size(K)) == 3 % finite horizon
        u(:,i) = K(:,:,i)*x(:,i);
    else % infinite horizon
        u(:,i) = K*x(:,i);
    end
        
    x(:,i+1) = Ad*x(:,i) + Bd*u(:,i);
end

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
% plotting outcome
figure;
plot(t,x1);
legend('system state x1');
figure;
plot(t,x2);
legend('system state x2');
figure;
plot(t,x3);
legend('system state x3');
figure;
plot(t(1:end-1),u);
legend('control input');

%% Trajectory tracking test

close all;
% track the ramp trajectory
s1 = sin(2*pi*t);
s2 = t;
s3 = t;
s = [s1;s2;s3];
s0 = s(:,1);

% form the time varying matrices Abar and Bbar
Abar = zeros(dimx+1,dimx+1,N);
Bbar = zeros(dimx+1,dimu,N);
for i = 1:N
    Abar(:,:,i) = [Ad, Ad*s(:,i) - s(:,i+1); ...
                   zeros(1,dimx), 0];
    Bbar(:,:,i) = [Bd; 0];
end
            
MODE.N = N;
MODE.LTI = false;
Q = 100*diag([1,0,0,0]);
R = eye(dimu);

K = LQR(Q,R,Abar,Bbar,MODE);

% simulate system
u = zeros(dimu,N);
x0 = zeros(dimx,1);
e0 = x0 - s0;
ebar(:,1) = [e0; 1];
xtest = zeros(dimx,length(t));
xtest(:,1) = x0;

for i = 1:N
    if length(size(K)) == 3 % finite horizon
        u(:,i) = K(:,:,i)*ebar(:,i);
    else % infinite horizon
        u(:,i) = K*ebar(:,i);
    end
    
    %xtest(:,i+1) = Ad * xtest(:,i) + Bd * u(:,i);
    ebar(:,i+1) = Abar(:,:,i)*ebar(:,i) + Bbar(:,:,i)*u(:,i);
end

% extract the signals from error
x = ebar(1:3,:) + s;
x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
% plotting outcome
figure;
plot(t,x1,'-',t,s1,'-.');
legend('system state x1', 'desired trajectory s1');
figure; 
plot(t,x2,'-',t,s2,'-.');
legend('system state x2','desired trajectory s2');
figure;
plot(t,x3,'-',t,s3,'-.');
legend('system state x3','desired trajectory s1');
figure;
plot(t(1:end-1),u);
legend('control input');