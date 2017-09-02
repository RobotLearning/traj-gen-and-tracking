%% Test script for LQR

clc; close all; clear classes; 

% 3D Linear Time Invariant test for LQR
% Example taken from http://www.egr.msu.edu/classes/me851/jchoi/lecture/Lect_14.pdf

dimx = 3;
dimu = 1;
dimy = 1;
% simulation variables
t0 = 0;
tf = 10;
h = 0.02;
t = t0:h:tf;
N = length(t)-1;

% system and weighting matrices
Q = 100*eye(dimy);
R = eye(dimu);
% continuous time matrices
A = [0 1 0; 0 0 1; -0.4 -4.2 -2.1];
B = [0;0;1];
C = [1 0 0];
% trick to get discrete time versions
Mat = [A, B; zeros(dimu, dimx + dimu)];
MD = expm(h * Mat);
Ad = MD(1:dimx,1:dimx);
Bd = MD(1:dimx,dimx+1:end);

% optimal feedback matrix
lqr = LQR(Q,R,Q,A,B,C,N,h);
K = lqr.computeFinHorizonLTI();

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
s = t.^2;
%s = sin(pi/6*t);
s0 = s(1);
% form sbar 
sbar = C'*((C*C')\s);

% form the time varying matrices Abar and Bbar
Abar = zeros(dimx+1,dimx+1,N);
Bbar = zeros(dimx+1,dimu,N);
for i = 1:N
    Abar(:,:,i) = [Ad, Ad*sbar(:,i) - sbar(:,i+1); ...
                   zeros(1,dimx), 1];
    Bbar(:,:,i) = [Bd; 0];
end
            
Q = 100*eye(dimy);
R = 0.01*eye(dimu);

% optimal feedback matrix
lqr = LQR(Q,R,Q,A,B,C,N,h);
[K,uff] = lqr.computeFinHorizonTracking(sbar);
Kbar = lqr.computeFinHorizonTrackingStateFeedback(sbar);

% simulate system
x0 = zeros(dimx,1);
x = zeros(dimx,N);
x(:,1) = x0;
u = zeros(dimu,N);
for i = 1:N
    u(:,i) = K(:,:,i)*x(:,i) + uff(:,i);
    x(:,i+1) = Ad*x(:,i) + Bd*u(:,i);
end

x1 = x(1,:);

%simulate system
u = zeros(dimu,N);
x0 = zeros(dimx,1);
e0 = x0 - s0;
ebar(:,1) = [e0; 1];
xtest = zeros(dimx,length(t));
xtest(:,1) = x0;
for i = 1:N
    if length(size(Kbar)) == 3 % finite horizon
        u(:,i) = Kbar(:,:,i)*ebar(:,i);
    else % infinite horizon
        u(:,i) = Kbar*ebar(:,i);
    end
    
    %xtest(:,i+1) = Ad * xtest(:,i) + Bd * u(:,i);
    ebar(:,i+1) = Abar(:,:,i)*ebar(:,i) + Bbar(:,:,i)*u(:,i);
end

% extract the signals from error
xx = ebar(1:3,:) + sbar;
x11 = xx(1,:);

% plotting outcome
figure;
plot(t,x1,'-',t,s,'-.');
legend('system state x1', 'desired trajectory s');
figure;
plot(t,x11,'-',t,s,'-.');
legend('system state x11', 'desired trajectory s');
figure;
plot(t(1:end-1),u);
legend('control input');

%% Test error form of LQR

% optimal feedback matrix
lqr = LQR(Q,R,Q,A,B,C,N,h);
[K,uff] = lqr.computeErrorForm(sbar);

% simulate system
x0 = ones(dimx,1);
x = zeros(dimx,N);
x(:,1) = x0;
u = zeros(dimu,N);
for i = 1:N
    u(:,i) = K(:,:,i)*(x(:,i)-sbar(:,i)) + uff(:,i);
    x(:,i+1) = Ad*x(:,i) + Bd*u(:,i);
end

x1e = x(1,:);
% plotting outcome
figure;
plot(t,x1e,'-',t,s,'-.');
legend('system state x1e', 'desired trajectory s');
title('Error form of LQR');

%% Test LQR with DMP for discrete system
%{
clc; clear; close all;
dimx = 2;
dimu = 1;
dimy = 1;

% simulation variables
t0 = 0;
tf = 1;
h = 0.02;
t = t0:h:tf;
N = length(t)-1;
numIt = 100;

% system and weighting matrices
Q = 100;
R = 1 * eye(dimu);
% discrete time matrices
A = [0 1; 1.8 -0.81];
B = [0;1];
%C = [0 1; -1/h 1/h];
C = [0 1];

% create the structures
SIM.discrete = true;
SIM.dimx = dimx;
SIM.dimu = dimu;
SIM.dimy = dimy;
SIM.h = h;
SIM.eps_m = 0;
SIM.int = 'Euler';
PAR.Ad = A;
PAR.Bd = B;
PAR.C = C;
CON = [];
COST.Q = Q;
COST.R = R;
lin = Linear(PAR,CON,COST,SIM);

% track smoothed step
ref = 1./(1+exp(-(100*t-20)/10));

% initialize states and inputs
x0 = [ref(1);ref(1)];
y0 = C * x0;
% create yin with zero velocity
yin = [ref(1);0];   

% create DMP trajectory and execute LQR
goal = ref(:,end);
numbf = 40;
[dmp,s] = lin.dmpTrajectory(t,numbf,goal,yin,ref);
s = s(1,:);
% form sbar 
sbar = C'*((C*C')\s);
% optimal feedback matrix
lqr = LQR(Q,R,Q,A,B,C,N,h,true);
[K,uff] = lqr.computeErrorForm(sbar);

% simulate system
x0 = [ref(1);ref(1)];
x = zeros(dimx,N);
x(:,1) = x0;
y0 = C * x0;
y = zeros(dimy,N);
u = zeros(dimu,N);
for i = 1:N
    u(:,i) = K(:,:,i)*(x(:,i)-sbar(:,i)) + uff(:,i);
    x(:,i+1) = A*x(:,i) + B*u(:,i);
    y(:,i) = C*x(:,i);
end
y(:,i+1) = C*x(:,i+1);

% plotting outcome
figure;
plot(t,y(1,:),'-',t,s(1,:),'-.');
legend('system state', 'desired trajectory s');
title('Error form of LQR');

% [y,us] = lin.observeWithDMPFeedback(dmp,traj,x0);
% traj.addPerformance(us,y,lin.COST,'LQR'); 
% 
% lin.plot_inputs(traj);
% lin.plot_outputs(traj);
%}