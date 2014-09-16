%% Simulate trajectories for the two-wheeled robot kinematical model

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

%% Set system parameters and constraints

% parameter values of the experimental setup
% radius of the robot
PAR.wheel1.radius = 0.2; % meter
PAR.wheel2.radius = 0.2;
PAR.length = 0.8;

% constraints on the system dynamics
CON.state.x.max = 5; 
CON.state.x.min = -5;
CON.state.y.max = 5;
CON.state.y.min = -5;
CON.state.theta.max = 2*pi;
CON.state.theta.min = -2*pi;
CON.wheel1.u.max = 200;
CON.wheel1.u.min = -200;
CON.wheel2.u.max = 200;
CON.wheel2.u.min = -200;
CON.wheel1.udot.max = 100;
CON.wheel1.udot.min = -100;
CON.wheel2.udot.max = 100;
CON.wheel2.udot.min = -100;

% Simulation Values 
% dimension of the x vector
SIM.dimx = 3;
% dimension of the control input
SIM.dimu = 2;
% time step h 
SIM.h = 0.01;
% noise and initial error
SIM.eps = 0.03;
SIM.eps_M = 0.005;
% integration method
SIM.int = 'Euler';

% cost structure
% only penalize positions
COST.Q = diag([1,1,0]);

% initialize model
TW = TwoWheeledCar(PAR,CON,COST,SIM);

%% Create desired trajectory 

dim_u = SIM.dimu;
dim_x = SIM.dimx;
h = SIM.h;
tfin = 1;
t = h:h:tfin;
N = length(t) - 1;
% create a pre-nominal sine-curve 
s(1,:) = t;
s(2,:) = sin(2*pi*t);
s(3,1:end-1) = atan2(diff(s(2,:)),diff(s(1,:)));
s(3,end) = s(3,1); % since trajectory is periodical after t = 1

Traj = TW.trajectory(t,s);

%% Create desired trajectory and linearized model

% Construct linearized model around the desired trajectory
% structure containing matrices returned
STR = construct_model(t,s,PAR,Traj.unom,CON);

% lifted-domain noise covariance
M = SIM.eps_M * eye(N*dim_x); % covariance of process x measurement noise
Omega0 = SIM.eps * eye(N*dim_x); % initial covariance of d noise

%% Iterative Learning Control

w = ones(1,N+1); % equal weighting
% lsq-cost for trajectory deviation
lsq_cost = @(x1,x2,w) norm(x1([1,2],:)-x2([1,2],:),2); 

d0_hat = 0;
P0 = eye(N*dim_x);

%%%%%%%%%%%%%%%%% CONSTRUCT SCALING MATRIX S %%%%%%%%%%%%%
% penalize only deviations given by the matrix
Sx = eye(dim_x); % scaling matrix
Sw = diag([1, 1, 0]); % weighing trajectory deviation by these values
w_l = 0.0; % put weight on last w_l second 
% emphasizing the performance objective of reaching upright position
last = w_l * round(1/h);
Sx1 = cell(1,N - last);
[Sx1{:}] = deal(Sx * Sw);
Sx1 = blkdiag(Sx1{:});
Sx2 = cell(1,last);
% termination matrix
Tw = diag(ones(1,dim_x));
mat = Tw * Sx * Sw; 
[Sx2{:}] = deal(mat);
Sx2 = blkdiag(Sx2{:});
S = blkdiag(Sx1,Sx2);

%%%%%%%%%%%%%%%% PASS TO STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%
% Pass necessary parameters through structure
STR.dim_x = dim_x;
STR.dim_u = dim_u;
STR.COV.Omega = Omega0;
STR.COV.M = M;
STR.CON = CON;
STR.PAR = PAR;
STR.d0 = d0_hat * ones(N*dim_x,1);
STR.P0 = P0;
STR.w = w;
STR.lsq_cost = lsq_cost;
STR.S = S;

u_trj = Traj.unom;
%%%%%%%%%%%%% ITERATIVE LEARNING CONTROL LOOP %%%%%%%%%%%%%%%%%
fun = @robotTwoWheelsError;
[err, u_app] = ILC(t,s,u_trj,fun,STR,10);
% errors can be reused in some other code