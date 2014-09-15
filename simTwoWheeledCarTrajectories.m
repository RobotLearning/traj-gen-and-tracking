%% Set system parameters

clc; clear, close all;

%%%%%%%%%%%%% Define System and Constraint Parameters %%%%%%%%
% parameter values of the experimental setup
% radius of the robot
PAR.R1 = 0.2; % meter
PAR.R2 = 0.2;
PAR.d = 0.8;
% constraints on the system dynamics
CON.x_cnstr = 5; %0.5 or 1
CON.y_cnstr = 5;
CON.w1_cnstr = 100;
CON.w2_cnstr = 100;

%%%%%%%%%%%%%%%%%%%%%%% Simulation Values %%%%%%%%%%%%%%%%%
% dimension of the x vector
dim = 3;
% dimension of the control input
dim_u = 2;
% time step h 
h = 0.01;
tfinal = 1;
t = 0:h:tfinal;
N = tfinal/h + 1; 
Nu = N-1; % one less than the length of t
% noise and initial error
eps = 0.03;
eps_M = 0.005; % covar of the lifted-domain noise
d0_hat = 0;
P0 = eye(N*dim);
% state trajectories
x = zeros(dim, N);
x_0 = 0; y_0 = 0; phi_0 = pi/6; % initial conditions

%% Create desired trajectory and linearized model

% generate trajectory from script
traj_gen;

% make room for learning
CON.w1_cnstr = 150; 
CON.w2_cnstr = 150; 

% Construct linearized model around the desired trajectory
% structure containing matrices returned
STR = construct_model(t,x,PAR,u_trj,CON);

% lifted-domain noise covariance
M = eps_M * eye(N*dim); % covariance of process x measurement noise
Omega0 = eps * eye(N*dim); % initial covariance of d noise

%% Iterative Learning Control

w = ones(1,N+1); % equal weighting
% lsq-cost for trajectory deviation
lsq_cost = @(x1,x2,w) norm(x1([1,2],:)-x2([1,2],:),2); 

%%%%%%%%%%%%%%%%% CONSTRUCT SCALING MATRIX S %%%%%%%%%%%%%
% penalize only deviations given by the matrix
Sx = eye(dim); % scaling matrix
Sw = diag([1, 1, 0]); % weighing trajectory deviation by these values
w_l = 0.0; % put weight on last w_l second 
% emphasizing the performance objective of reaching upright position
last = w_l * round(1/h);
Sx1 = cell(1,N - last);
[Sx1{:}] = deal(Sx * Sw);
Sx1 = blkdiag(Sx1{:});
Sx2 = cell(1,last);
% termination matrix
Tw = diag(ones(1,dim));
mat = Tw * Sx * Sw; 
[Sx2{:}] = deal(mat);
Sx2 = blkdiag(Sx2{:});
S = blkdiag(Sx1,Sx2);

%%%%%%%%%%%%%%%% PASS TO STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%
% Pass necessary parameters through structure
STR.dim = dim;
STR.dim_u = dim_u;
STR.COV.Omega = Omega0;
STR.COV.M = M;
STR.CON = CON;
STR.PAR = PAR;
STR.d0 = d0_hat * ones(N*dim,1);
STR.P0 = P0;
STR.w = w;
STR.lsq_cost = lsq_cost;
STR.S = S;

%%%%%%%%%%%%% ITERATIVE LEARNING CONTROL LOOP %%%%%%%%%%%%%%%%%
fun = @robotTwoWheelsError;
[err, u_app] = ILC(t,x,u_trj,fun,STR,10);
% errors can be reused in some other code