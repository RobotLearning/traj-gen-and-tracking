%% Testing (general) Inverse Kinematics function on Barrett WAM

clc; clear; close all; dbstop if error;

%% Initialize Barrett WAM

initializeWAM;

%% Initialize arm posture

% initialize the arm with zero velocity on the right hand side
q = [1.8; -0.2; -0.1; 1.8; -1.57; 0.1; 0.3] + 0.1*randn(N_DOFS,1);
qd = zeros(N_DOFS,1);

n = 5;
s2 = 0.01;
qs = repmat(q,1,n) + sqrt(s2)*randn(N_DOFS,n);
qds = repmat(qd,1,n);

[x,xd,o] = wam.kinematics([qs;qds]);

q2 = wam.invKinematics(x,o,q);
[x2,xd2,o2] = wam.kinematics([q2;qds]);

%scatter3(x(1,:),x(2,:),x(3,:))