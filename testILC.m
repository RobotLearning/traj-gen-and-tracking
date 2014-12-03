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

%% Does the same blowup phenomenon occur also for model-free regression WILC?

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
C = [0 1];

% create the structures
SIM.discrete = true;
SIM.dimx = dimx;
SIM.dimu = dimu;
SIM.dimy = dimy;
SIM.h = h;
SIM.eps = 0;
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

% create trajectory and execute LQR
%traj = lin.generateInputs(t,ref);
%[y,us] = lin.observeWithFeedback(traj,x0);

% create DMP trajectory and execute LQR
[traj,dmp] = lin.generateDMP(t,yin,ref);
[y,us] = lin.observeWithDMPFeedback(dmp,traj,x0);
traj.addPerformance(us,y,lin.COST,'LQR'); 

% or create zero input
% us = zeros(dimu,N);
% traj = Trajectory(t,ref,us,[]);
% y = lin.observe(t,x0,us);
% traj.addPerformance(us,y,lin.COST,'zeros'); 

% Create an ilc controller
% create the simpler ilc
%ilc = bILC(traj);
ilc = wILC(lin,traj);
num_trials = 10;

for i = 1:num_trials
    % update the weights of the dmp
    %us = ilc.feedforward(traj,y);
    ilc.feedforward(dmp,traj,y);
    % get the measurements
    %y = lin.observe(t,x0,us);
    % get the measurements
    [y,us] = lin.observeWithDMPFeedback(dmp,traj,x0);
    traj.addPerformance(us,y,lin.COST,ilc);

end

lin.plot_inputs(traj);
lin.plot_outputs(traj);

figure;
plot(1:num_trials,ilc.error);
title('Squared-2-Norm of ILC error');
legend(ilc.name);

%% More complicated example
% % Example taken from http://www.egr.msu.edu/classes/me851/jchoi/lecture/Lect_14.pdf
% 
% close all; clear; clc;
% dimx = 3;
% dimu = 1;
% dimy = 2;
% 
% % simulation variables
% t0 = 0;
% tf = 1;
% h = 0.02;
% t = t0:h:tf;
% N = length(t)-1;
% 
% % system and weighting matrices
% Q = 100*eye(dimy);
% R = 0.01*eye(dimu);
% % continuous time matrices
% A = [0 1 0; 0 0 1; -0.4 -4.2 -2.1];
% B = [0;0;1];
% C = [1 0 0; 0 1 0];
% % trick to get discrete time versions
% Mat = [A, B; zeros(dimu, dimx + dimu)];
% MD = expm(h * Mat);
% Ad = MD(1:dimx,1:dimx);
% Bd = MD(1:dimx,dimx+1:end);
% 
% % Create the 2-norm cost function
% cost.fnc = @(y,s) diag((y-s)'*Q*(y-s));
% 
% % track the reference trajectory
% %r = sin(pi/6*t);
% r = t.^2;
% r = [r; 2*t];
% 
% % initialize states and inputs
% x0 = [r(:,1);0];
% x = zeros(dimx,N+1);
% x(:,1) = x0;
% y(:,1) = C*x0;
% u0 = zeros(dimu,1);
% u = zeros(dimu,N);
% u(:,1) = u0;
% 
% % evolve model with zero input
% for j = 1:N
%     x(:,j+1) = Ad*x(:,j) + Bd*u(:,j);
%     y(:,j+1) = C*x(:,j+1);
% end
% 
% % form a trajectory
% trj = Trajectory(t,r,u,[]);
% trj.addPerformance(u,y,cost,'zeros');
% 
% % create the simpler ilc
% %ilc = bILC(trj);
% 
% % create the structures
% SIM.discrete = false;
% SIM.dimx = dimx;
% SIM.dimu = dimu;
% SIM.dimy = dimy;
% SIM.h = h;
% SIM.eps = 0;
% SIM.int = 'Euler';
% PAR.A = A;
% PAR.B = B;
% PAR.C = C;
% CON = [];
% COST.Q = Q;
% COST.R = R;
% lin = Linear(PAR,CON,COST,SIM);
% % create the model based ilc
% ilc = mILC(lin,trj);
% ilc.u_last = u;
% num_trials = 1;
% 
% % Perform ILC updates
% for i = 1:num_trials
%     % Update the controls
%     u = ilc.feedforward(trj,y);
%     % Evolve system with both ILC inputs
%     for j = 1:N
%         x(:,j+1) = Ad*x(:,j) + Bd*u(:,j);
%         y(:,j+1) = C*x(:,j+1);
%     end
%     % Add performances
%     trj.addPerformance(u,y,cost,ilc);
% 
% end
% 
% figure(1);
% plot(1:num_trials,ilc.error);
% title('Squared-2-Norm of ILC error');
% figure(2);
% plot(t,y(1,:),'.-',t,r(1,:),'-');
% title('Last iteration result');
% legend('ILC trajectory','Reference');