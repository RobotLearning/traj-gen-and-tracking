%% Test script for ILC (Arimoto-style and Model-based ILC tests)

%# store breakpoints
tmp = dbstatus;
save('tmp.mat','tmp')

%# clear all
close all
clear classes %# clears even more than clear all
clc

%# reload breakpoints
load('tmp.mat')
dbstop(tmp)

%# clean up
clear tmp
delete('tmp.mat')

%% ILC Example 
% taken from the ILC survey by Bristow et al.

dimx = 2;
dimu = 1;

% simulation variables
t0 = 0;
tf = 50;
h = 1;
t = t0:h:tf;
N = length(t)-1;
numIt = 100;

% system and weighting matrices
Q = 1;
R = eye(dimu);
% discrete time matrices
A = [0 1; 1.8 -0.81];
B = [0;1];
C = [0 1];

% Create the 2-norm cost function
cost.fnc = @(y,s) diag((y-s)'*Q*(y-s));

% track smoothed step
r = 1./(1+exp(-(2*t-20)/10));

% initialize states and inputs
x0 = [r(1);r(1)];
x = zeros(dimx,N+1);
x(:,1) = x0;
y(:,1) = C*x0;
u0 = zeros(dimu,1);
u = zeros(dimu,N);
u(:,1) = u0;

% evolve model with zero input
for j = 1:N
    x(:,j+1) = A*x(:,j) + B*u(:,j);
    y(:,j+1) = C*x(:,j+1);
end

% form a trajectory
trj = Trajectory(t,r,u,[]);
trj.addPerformance(u,y,cost,'zeros');

% create the simpler ilc
ilc = bILC(trj);

% Perform ILC updates
for i = 1:numIt
    % Update the controls
    u = ilc.feedforward(trj,y);
    % Evolve system with both ILC inputs
    for j = 1:N
        x(:,j+1) = A*x(:,j) + B*u(:,j);
        y(:,j+1) = C*x(:,j+1);
    end
    % Add performances
    trj.addPerformance(u,y,cost,ilc);

end

figure(1);
plot(1:numIt,ilc.error);
title('Squared-2-Norm of ILC error');
figure(2);
plot(t,y,'.-',t,r,'-');
title('Last iteration result');
legend('ILC trajectory','Reference');

%% More complicated example
% Example taken from http://www.egr.msu.edu/classes/me851/jchoi/lecture/Lect_14.pdf

close all; clear; clc;
dimx = 3;
dimu = 1;
dimy = 2;

% simulation variables
t0 = 0;
tf = 1;
h = 0.02;
t = t0:h:tf;
N = length(t)-1;

% system and weighting matrices
Q = 100*eye(dimy);
R = 0.01*eye(dimu);
% continuous time matrices
A = [0 1 0; 0 0 1; -0.4 -4.2 -2.1];
B = [0;0;1];
C = [1 0 0; 0 1 0];
% trick to get discrete time versions
Mat = [A, B; zeros(dimu, dimx + dimu)];
MD = expm(h * Mat);
Ad = MD(1:dimx,1:dimx);
Bd = MD(1:dimx,dimx+1:end);

% Create the 2-norm cost function
cost.fnc = @(y,s) diag((y-s)'*Q*(y-s));

% track the reference trajectory
%r = sin(pi/6*t);
r = t.^2;
r = [r; 2*t];

% initialize states and inputs
x0 = [r(:,1);0];
x = zeros(dimx,N+1);
x(:,1) = x0;
y(:,1) = C*x0;
u0 = zeros(dimu,1);
u = zeros(dimu,N);
u(:,1) = u0;

% evolve model with zero input
for j = 1:N
    x(:,j+1) = Ad*x(:,j) + Bd*u(:,j);
    y(:,j+1) = C*x(:,j+1);
end

% form a trajectory
trj = Trajectory(t,r,u,[]);
trj.addPerformance(u,y,cost,'zeros');

% create the simpler ilc
%ilc = bILC(trj);

% create the structures
SIM.discrete = false;
SIM.dimx = dimx;
SIM.dimu = dimu;
SIM.dimy = dimy;
SIM.h = h;
SIM.eps = 0;
SIM.int = 'Euler';
PAR.A = A;
PAR.B = B;
PAR.C = C;
CON = [];
COST.Q = Q;
COST.R = R;
lin = Linear(PAR,CON,COST,SIM);
% create the model based ilc
ilc = mILC(lin,trj);
ilc.u_last = u;
num_trials = 1;

% Perform ILC updates
for i = 1:num_trials
    % Update the controls
    u = ilc.feedforward(trj,y);
    % Evolve system with both ILC inputs
    for j = 1:N
        x(:,j+1) = Ad*x(:,j) + Bd*u(:,j);
        y(:,j+1) = C*x(:,j+1);
    end
    % Add performances
    trj.addPerformance(u,y,cost,ilc);

end

figure(1);
plot(1:num_trials,ilc.error);
title('Squared-2-Norm of ILC error');
figure(2);
plot(t,y(1,:),'.-',t,r(1,:),'-');
title('Last iteration result');
legend('ILC trajectory','Reference');