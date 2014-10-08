%% Test script for ILC (Arimoto-style)

clc; close all; clear classes; 
dimx = 3;
dimu = 1;

% simulation variables
t0 = 0;
tf = 0.1;
h = 0.02;
t = t0:h:tf;
N = length(t)-1;

% system and weighting matrices
Q = 1*diag([1,0,0]);
R = eye(dimu);
% continuous time matrices
A = [0 1 0; 0 0 1; -0.4 -4.2 -2.1];
B = [0;0;1];

% create the structures
SIM.dimx = dimx;
SIM.dimu = dimu;
SIM.h = h;
SIM.int = 'Euler';
PAR.A = A;
PAR.B = B;
COST.Q = Q;
COST.R = R;
model = Linear(PAR,[],COST,SIM);

% track sin and ramp
s = [sin(2*pi*t);t.^2;t];
x0 = zeros(dimx,1);

% execute LQR
[x,unom,K] = model.lqr(t,x0,s);

% Update the trajectory class
trj = Trajectory(t,[],s,unom);
[x2,trj2] = model.evolveWithFeedback(trj,x0,K);
trj.addPerformance(unom,x,model.COST,'LQR');

model.plot_controls(trj);
model.plot_states(trj);

% Create an ilc controller
ilc = bILC(trj);
for i = 1:5
    % Update the controls
    u = ilc.feedforward(trj,x);
    % TODO: test the evolution here!
    x = model.evolve(t,x0,u);
    
    trj.addPerformance(u,x,model.COST,ilc);

    model.plot_controls(trj);
    model.plot_states(trj);
end