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
SIM.int = 'RK4'; %'Symplectic Euler';
% reference trajectory in joint space?
SIM.jref = true;

% load (nominal) parameter values for robot dynamics
loadActualBarrettValues;
%loadNominalBarrettValues;

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
MAX_VEL = 10;
MAX_ACC = 200;
SLACK = 0.05;
CON.q.max = [2.60; 2.00; 2.80; 3.10; 1.30; 1.60; 2.20] - SLACK;
CON.q.min = [-2.60; -2.00; -2.80; -0.90; -4.80; -1.60; -2.20] + SLACK;
CON.qd.max = MAX_VEL * ones(7,1);
CON.qd.min = -MAX_VEL * ones(7,1);
CON.qdd.max = MAX_ACC * ones(7,1);
CON.qdd.min = -MAX_ACC * ones(7,1);
CON.u.max = [75; 125; 39; 30; 3; 4; 1];
CON.u.min = -CON.u.max; % ONLY Q and QD constraints are used in opt.

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

% initialize the arm with zero velocity on the right hand side
q_right = [1.0; -0.2; -0.1; 1.8; -1.57; 0.1; 0.3];
q_left = [-1.0; 0.0; 0.0; 1.5; -1.57; 0.1; 0.3];
q_centre = [0.0; 0.0; 0.0; 1.5; -1.75; 0.0; 0.0];

try
    switch robot_side
        case 'RIGHT'
            q0 = q_right;
        case 'LEFT'
            q0 = q_left;
        case 'CENTRE'
            q0 = q_centre;
        otherwise
            error('Pose not identified!');
    end
catch
    warning('Robot side was not set. Init to right');
    q0 = q_right;
end
qd0 = zeros(7,1);
Q0 = [q0;qd0];

% to help with inverse kin
wam.regressOnFinalJointsFromDemo();

% take points from the workspace boundary
wam.buildWorkspace();