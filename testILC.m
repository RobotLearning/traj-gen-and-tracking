%% Test script for ILC (Arimoto-style)

clc; close all; clear classes; 

%% Nonmonotonic ILC Example 
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

% track smoothed step
s = 1./(1+exp(-(2*t-20)/10));

% initialize states and inputs
x0 = [0;s(1)];
x = zeros(dimx,N+1);
x(:,1) = x0;
y(:,1) = C*x0;
u0 = zeros(dimu,1);
u = zeros(dimu,N);
u1 = u; u2 = u;
u1(:,1) = u0;
u2(:,1) = u0;

% evolve model
for j = 1:N
    x(:,j+1) = A*x(:,j) + B*u(:,j);
    y(:,j+1) = C*x(:,j+1);
end
x1 = x; x2 = x;
y1 = y; y2 = y;

% form a trajectory
trj = Trajectory(t,[],s,u);

% Create two ilc controllers
bilc = bILC(trj);

% Create the 2-norm cost function
cost.fnc = @(y,s) diag((y-s)'*Q*(y-s));

% Perform ILC updates
for i = 1:numIt
    % Update the controls
    u1 = bilc.feedforward(trj,y1);
    % Evolve system with both ILC inputs
    for j = 1:N
        x1(:,j+1) = A*x1(:,j) + B*u1(:,j);
        y1(:,j+1) = C*x1(:,j+1);
    end
    % Add performances
    trj.addPerformance(u1,y1,cost,bilc);

end

figure(1);
plot(1:numIt,bilc.error);
title('Squared-2-Norm of ILC error');
legend('Basic PD-type ILC');
figure(2);
plot(t,y1,'.-',t,s,'-');
title('Last iteration result');
legend('bILC trajectory','Reference');

%% More complicated example

% clc; close all; clear classes; 
% dimx = 3;
% dimu = 1;
% 
% % simulation variables
% t0 = 0;
% tf = 0.1;
% h = 0.02;
% t = t0:h:tf;
% N = length(t)-1;
% 
% % system and weighting matrices
% Q = 1*diag([1,0,0]);
% R = eye(dimu);
% % continuous time matrices
% A = [0 1 0; 0 0 1; -0.4 -4.2 -2.1];
% B = [0;0;1];
% 
% % create the structures
% SIM.dimx = dimx;
% SIM.dimu = dimu;
% SIM.h = h;
% SIM.int = 'Euler';
% PAR.A = A;
% PAR.B = B;
% COST.Q = Q;
% COST.R = R;
% model = Linear(PAR,[],COST,SIM);
% 
% % track sin and ramp
% s = [sin(2*pi*t);t.^2;t];
% x0 = zeros(dimx,1);
% 
% % execute LQR
% [x,unom,K] = model.lqr(t,x0,s);
% 
% % Update the trajectory class
% trj = Trajectory(t,[],s,unom);
% [x2,trj2] = model.evolveWithFeedback(trj,x0,K);
% trj.addPerformance(unom,x,model.COST,'LQR');
% 
% model.plot_controls(trj);
% model.plot_states(trj);
% 
% % Create an ilc controller
% ilc = bILC(trj);
% for i = 1:5
%     % Update the controls
%     u = ilc.feedforward(trj,x);
%     % TODO: test the evolution here!
%     x = model.evolve(t,x0,u);
%     
%     trj.addPerformance(u,x,model.COST,ilc);
% 
%     model.plot_controls(trj);
%     model.plot_states(trj);
% end