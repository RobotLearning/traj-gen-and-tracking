%% Generic script to initialize Barrett WAM 

N_DOFS = 7;
% Simulation Values 
% system matrices are continous
SIM.discrete = false;
% learn in cartesian space
SIM.cartesian = false;
% dimension of the x vector
SIM.dimx = 2*N_DOFS;
% dimension of the output y
SIM.dimy = 2*N_DOFS;
% dimension of the control input
SIM.dimu = N_DOFS;
% time step h 
SIM.h = 0.002; % 500 Hz recorded data
% measurement noise covariance (eps_m * eye)
SIM.eps_m = 0e-10;
% integration method
SIM.int = 'Symplectic Euler';
% reference trajectory in joint space?
SIM.jref = true;

% load (nominal) parameter values for robot dynamics
%loadActualBarrettValues;
loadNominalBarrettValues;

% We observe all joints and all joint velocities
PAR.links = links;
PAR.link0 = link0;
PAR.eff = eff;
PAR.basec = basec;
PAR.baseo = baseo;
PAR.uex = uex;
PAR.uex0 = uex0;
PAR.C = eye(SIM.dimy,SIM.dimx);

% form constraints
CON = [];

% cost struc
Q1 = 1*diag([ones(1,4),1*ones(1,3),1*ones(1,4),1*ones(1,3)]);
% only penalize positions
Q2 = 1*diag([ones(1,4),1*ones(1,3),0.0*ones(1,4),0.0*ones(1,3)]);
COST.Q = Q1;
COST.R = 0.05 * eye(SIM.dimu);

% initialize model
wam = BarrettWAM(PAR,CON,COST,SIM);

% PD control defined here
PD = zeros(N_DOFS,2*N_DOFS);
PD(1,1) = -200;
PD(1,N_DOFS+1) = -7.0;
PD(2,2) = -300;
PD(2,N_DOFS+2) = -15.0;
PD(3,3) = -100;
PD(3,N_DOFS+3) = -5.0;
PD(4,4) = -50;
PD(4,N_DOFS+4) = -2.5;
PD(5,5) = -10;
PD(5,N_DOFS+5) = -0.3;
PD(6,6) = -10;
PD(6,N_DOFS+6) = -0.3;
PD(7,7) = -2.5;
PD(7,N_DOFS+7) = -0.075;

%% Initialize arm posture

% initialize the arm with zero velocity on the right hand side
q0 = [1.8; -0.2; -0.1; 1.8; -1.57; 0.1; 0.3];
qd0 = zeros(N_DOFS,1);
Q0 = [q0; qd0];
[x0,xd0] = wam.kinematics(Q0);
X0 = [x0;xd0];