%% Example to test wILC algorithm on linear systems with drift
%
% TODO: this does not run now!
%

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

%% Show that regression on states and state derivatives are equivalent

t0 = 0;
tf = 1;
h = 0.02;
t = t0:h:tf;
t = t(:);
N = length(t)-1;
s = 1./(1+exp(-(100*t-20)/10));
s = s(:);

% form the regression matrix
bfs = N+1;
Psi = zeros(N+1,bfs);
hh = ones(bfs,1) * bfs^(1.5);
c = linspace(t(1),t(end),bfs);
for i = 1:bfs
    Psi(:,i) = exp(-hh(i) * (t - c(i)).^2);
end

% usual regression
w = Psi \ s(:);
s1 = Psi * w;

% regression on derivatives
D = (diag(-ones(1,N),-1) + eye(N+1))/h;
S = h*tril(ones(N+1));
w = Psi \ (D * s(:));

sd = Psi * w;
s2 = S * sd;

%% wILC is equivalent to bILC + feedback

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

% create DMP trajectory and execute LQR
[traj,dmp] = lin.generateDMP(t,x0,ref);
[y,~] = lin.observeWithFeedbackErrorForm(traj,x0);
traj.addPerformance([],y,lin.COST,'LQR'); 

% create trajectory and execute LQR
%traj = lin.generateInputs(t,ref,1);
%[y,us] = lin.observeWithFeedbackErrorForm(traj,x0);
%traj.addPerformance(traj.s,y,lin.COST,'LQR'); 

% Create an ilc controller
% create the simpler ilc
%ilc = wILC(traj,lin,'dmp');
ilc = wILC(traj,lin,'wd');
num_trials = 0;

for i = 1:num_trials
    
    traj2 = ilc.feedforward(traj,[],y);
    % update the weights of the dmp
    %ilc.feedforward(traj,dmp,y);
    % get the measurements
    [y,~] = lin.observeWithFeedbackErrorForm(traj2,x0);
    %[y,~] = lin.observeWithFeedbackErrorForm(traj,x0,dmp);
    traj.addPerformance([],y,lin.COST,ilc);

end

lin.plot_outputs(traj);

%% Example taken from http://www.egr.msu.edu/classes/me851/jchoi/lecture/Lect_14.pdf

dimx = 3;
dimu = 1;
dimy = 1;

% simulation variables
t0 = 0;
tf = 1;
h = 0.02;
t = t0:h:tf;
N = length(t)-1;

% system and weighting matrices
Q = 100*eye(dimy);
%Q = diag([100,0,0]);
R = 0.01*eye(dimu);
% continuous time matrices
A = [0 1 0; 0 0 1; -0.4 -4.2 -2.1];
B = [0;0;1];
% partial observation model
%C = eye(dimy);
%C = [1 0 0; 0 1 0];
C = [1 0 0];

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

% track the sin trajectory
ref = t.^2;
%ref = [ref; 2*t];
%ref = [ref; 0, 2*ones(1,length(t)-1)];
x0 = zeros(dimx,1);
y0 = C * x0;
% create yin with zero velocity
yin = [y0;0];      

% create DMP trajectory and execute LQR
[traj,dmp] = lin.generateDMP(t,yin,ref);
[y,~] = lin.observeWithFeedbackErrorForm(traj,x0);
traj.addPerformance([],y,lin.COST,'LQR'); 

lin.plot_outputs(traj);

% Create an ilc controller
ilc = wILC(traj,lin,'t');
%ilc = wILC(traj,lin,'dmp');
num_trials = 10;

for i = 1:num_trials
    % update the weights of the dmp
    traj2 = ilc.feedforward(traj,[],y);
    %dmp = ilc.feedforward(traj,dmp,y);     
    % get the measurements
    [y,~] = lin.observeWithFeedbackErrorForm(traj2,x0);
    %[y,~] = lin.observeWithFeedbackErrorForm(traj,x0,dmp);
    traj.addPerformance([],y,lin.COST,ilc);
end

lin.plot_outputs(traj);