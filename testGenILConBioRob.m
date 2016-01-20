%% Test generalization of ILC on BioRob

clc; clear; close all; dbstop if error;

folder = [pwd, '/ILC_ProMP'];
data_folder = [folder '/common_data'];

%% Load BioRob nominal model

N_DOFS = 4;
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
SIM.h = 0.02; % 50 Hz recorded data
% measurement noise covariance (eps_m * eye)
SIM.eps_m = 0e-10;
% integration method
SIM.int = 'Symplectic Euler';
% reference trajectory in joint space?
SIM.jref = false;

% load (nominal) parameter values for robot dynamics
loadNominalBioRobValues;

% We observe all joints and all joint velocities
PAR.links = links;
PAR.eff = eff;
PAR.basec = basec;
PAR.baseo = baseo;
PAR.uex = uex;
PAR.uex0 = uex0;
PAR.C = eye(SIM.dimy,SIM.dimx);

% form constraints
deg = pi/180;
CON.link1.q = [-160 160]*deg;
CON.link2.q = [-360 360]*deg;
CON.link3.q = [-360 360]*deg;
CON.link4.q = [-360 360]*deg;

% cost struc
Q1 = 1*diag([ones(1,4),1*ones(1,3),1*ones(1,4),1*ones(1,3)]);
% only penalize positions
Q2 = 1*diag([ones(1,4),1*ones(1,3),0.0*ones(1,4),0.0*ones(1,3)]);
COST.Q = Q1;
COST.R = 0.05 * eye(SIM.dimu);

% initialize model
rob = BioRob(PAR,CON,COST,SIM);

%% Create or load trajectories

dt = SIM.h;
numRefs = 5;
tf = 3.0;
l = 0.25;
X0 = [0.2;-0.1];
A = [0.3;0.1];
B = A + [0;l];
C = A + [l;l];
D = A + [l;0];
squareVertices = [A,B,C,D];
y = [X0(1)*ones(1,4);
     squareVertices(1,:)];
z = [X0(2)*ones(1,4);
     squareVertices(2,:)];
phi = [pi/3;0];

% TODO: repair inverse kinematics, dynamics, implement forward kinematics

% for i = 1:numRefs
%     trj(i) = rob.genReachingTraj(tf,y(:,i),z(:,i),phi);
% end

% ===========================
% Create trajectories without plotting them.
% Find the IK and ID for different reaching positions.
% ==========================
if ~exist('ref_reaching_all.mat')
    refnames = {'A', 'B', 'C', 'D', 'E'};
    for m = 1:numel(refnames)
        refname       = refnames{m};
        ref_file_name = ['ref_' refname '.mat'];
        refc{m}       = create_reference(bioRob, 0.02, refname);
        ref           = refc{m};
        save([data_folder ref_file_name], 'ref');
    end
    % store the reaching positions as a mat file.
    save([data_folder 'ref_reaching_all.mat'], 'refc');
else
    % ===========================
    % Load and plot all reference trajectories
    % ==========================
    load([data_folder '/ref_reaching_all.mat']);
end

h = rob.plot_outputs(refc);
numRefs = numel(refc);

%% Learn torques with ILC for each trajectory

% run this for loop for each of the reaching positions
for k = 1:numRefs 
    
    prm.plot.refAnimation=1;
    
    fprintf('\n\n\n');
    fprintf('k: %g out of %g\n',k, nRefs);
    close all
    ref = refc{k};
    h = prepare_figures(dfv, ref, [], []);

    % animate the ideal tracking of the reference with zero errror
    if prm.plot.refAnimation
        figure(h.anim);
        title('Ideal tracking of the reference');
        animate_bioRob(bioRob, ref.q, 25);
    end

    % load controller parameters.
    % Without a second input argument, the feedforward inputs of the first
    % iteration will be set to zero.
    % Those feedforward inputs are in prm.ilc.u.
    prm = load_bioRob_controller(ref);
    
    % change to zero if you dont want to see animation after each iteration
    prm.plot.animation = 0; 
    
    % main learning loop    
    nMax = 30; % number of ILC iterations
    [iter, prm] = ilc_main_learning_loop([], nMax, prm, bioRob, h, ref);
   
    if 1 % enable saving
        refname = ref.name;
        filename = ['ilcsol_' refname '_' num2str( size(iter,2), 2 ) '.mat'];
        fprintf('Saving ILC solution with name: %s\n', filename)
        save([data_folder filename], 'iter', 'prm');
    end
   
end