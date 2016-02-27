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
MAX_VEL = 200;
MAX_ACC = 200;
SLACK = 0.05;
CON.q.max = [2.60; 2.00; 2.80; 3.10; 1.30; 1.60; 2.20] - SLACK;
CON.q.min = [-2.60; -2.00; -2.80; -0.90; -4.80; -1.60; -2.20] + SLACK;
CON.qd.max = MAX_VEL * ones(7,1);
CON.qd.min = -MAX_VEL * ones(7,1);
CON.qdd.max = MAX_ACC * ones(7,1);
CON.qdd.min = -MAX_ACC * ones(7,1);
CON.u.max = [75; 125; 39; 30; 3; 4; 1];
CON.u.min = -CON.u.max;

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
q0 = [1.0; -0.2; -0.1; 1.8; -1.57; 0.1; 0.3];
qd0 = zeros(7,1);
Q0 = [q0;qd0];

% to help with inverse kin
wam.regressOnFinalJointsFromDemo();

% construct robot workspace by using mesh
% tic;
% numPt = 5;
% for i = 1:N_DOFS
%     mesh{i} = linspace(CON.q.min(i),CON.q.max(i),numPt);
% end
% [x1,x2,x3,x4,x5,x6,x7] = ndgrid(mesh{:});
% Qworkspace = [x1(:),x2(:),x3(:),x4(:),x5(:),x6(:),x7(:)];
% [xw,~,~] = wam.kinematics([Qworkspace';zeros(7,numPt^7)]);
% Del = delaunay(xw(1,:),xw(2,:),xw(3,:));
% trisurf(Del,xw(1,:),xw(2,:),xw(3,:));
% save('BarrettWorkspace.mat','Del');
% toc
% 
% load('BarrettWorkspace.mat','Del');

%{
% Search for an initial posture with less jacobian condition number
numIter = 100;
s2 = 0.1;
qs = repmat(q0,1,numIter) + sqrt(s2)*randn(N_DOFS,numIter);
% choose posture around q0 minimizing jacobian

for i = 1:numIter
    J = wam.calcJacobian(qs(:,i));
    c(i) = cond(J);
end

[jacCond,ind] = min(c);
fprintf('Initial Jacobian condition: %f.\n',jacCond);
q0 = qs(:,ind);
%}