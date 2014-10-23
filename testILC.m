%% Test script for ILC (Arimoto-style and Model-based ILC tests)

clc; close all; clear classes; 

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
s = 1./(1+exp(-(2*t-20)/10));

% initialize states and inputs
x0 = [0;s(1)];
x = zeros(dimx,N+1);
x(:,1) = x0;
y(:,1) = C*x0;
u0 = zeros(dimu,1);
u = zeros(dimu,N);
u = u; u2 = u;
u(:,1) = u0;
u2(:,1) = u0;

% evolve model
for j = 1:N
    x(:,j+1) = A*x(:,j) + B*u(:,j);
    y(:,j+1) = C*x(:,j+1);
end

% form a trajectory
trj = Trajectory(t,s,u,[]);
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
plot(t,y,'.-',t,s,'-');
title('Last iteration result');
legend('ILC trajectory','Reference');