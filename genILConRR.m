%% Generalizing wILC inside the convex hull of goal points
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

%% Define constants and parameters

% Simulation Values 
% system is continous
SIM.discrete = false;
% learn in cartesian space
SIM.cartesian = false;
% dimension of the x vector
SIM.dimx = 4;
% dimension of the output y
SIM.dimy = 4;
% dimension of the control input
SIM.dimu = 2;
% time step h 
SIM.h = 0.01;
% noise and initial error
SIM.eps = 3e-10;
SIM.eps_d = 3e-10;
% integration method
SIM.int = 'Euler';
% trajectory in joint space?
SIM.jref = false;

% constants
gravity = 9.81;
% joint parameters
m1 = 1; %mass of first link, kg
m2 = 0.5; %mass of second link, kg
l1 = 0.50; %length of first link, m
l2 = 0.40; %length of second link, m
l_c1 = 0.25; %distance of first link's center of gravity to prev. joint, m
l_c2 = 0.20; %dist. of second link's c.oZ.g. to prev. joint, m
I1 = (1/12)*m1*l1^2; %assume thin rod moment of inertia around c.o.g.
I2 = (1/12)*m2*l2^2; %kg m^2
% motor parameters
J_a1 = 0.100; % actuator inertia of link 1
J_g1 = 0.050; % gear inertia of link1
J_m1 = J_a1 + J_g1;
J_a2 = 0.080; % actuator inertia of link 2
J_g2 = 0.040; % gear inertia of link2
J_m2 = J_a2 + J_g2;
% gear ratio typically ranges from 20 to 200 or more - comment from book
r_1 = 20;
r_2 = 20;
 
% pass it as a parameter structure
PAR.const.g = gravity;
PAR.link1.mass = m1;
PAR.link2.mass = m2;
PAR.link1.length = l1;
PAR.link2.length = l2;
PAR.link1.centre.dist = l_c1;
PAR.link2.centre.dist = l_c2;
PAR.link1.inertia = I1;
PAR.link2.inertia = I2;
PAR.link1.motor.inertia = J_m1;
PAR.link2.motor.inertia = J_m2;
PAR.link1.motor.gear_ratio = r_1;
PAR.link2.motor.gear_ratio = r_2;

% we should observe joint angles only
PAR.C = eye(SIM.dimy,SIM.dimx);

% form constraints
CON = [];

% cost structure
% only penalize positions
COST.Q = 100*diag([1,1,0,0]);
COST.R = 1 * eye(SIM.dimu);

% initialize model
rr = RR(PAR,CON,COST,SIM);

%% Apply ILC for DMPs

% create reference
h = SIM.h;
tin = 0; tfin = 1;
t = tin:h:tfin;
y_des = 0.4 + 0.4 * t;
x_des = 0.6 - 0.2 * t;
ref = [x_des; y_des]; % displacement profile 

%traj = rr.generateInputs(t,ref);
bfs = 100;
% trajectory generated in joint space
[traj,dmpOrigin] = rr.generateInputsWithDMP(t,bfs,ref); 

%save dmp weights for later use
w_origin = zeros(length(dmpOrigin),bfs);
for i = 1:length(dmpOrigin)
    w_origin(i,:) = dmpOrigin(i).FORCE.w;
end

% modify reference to obtain different dmps
% center of the region
g0 = traj.s(:,end);
q0 = traj.s(:,1);
% radius of region
r = 0.2;
% mesh of points
diff_g = r * [1 1 0 0; 1 -1 0 0; -1 1 0 0; -1 -1 0 0];

num_trials = 5;
dmp = cell(1,4);

for i = 1:4
    % take a point from the mesh
    g(:,i) = g0 + diff_g(i,:)';
    [dmp{i},jNew] = adaptDMP(q0,g(:,i),dmpOrigin,w_origin);
    [~,ref] = rr.kinematics(jNew);

    traj = rr.generateInputs(t,ref);
    qin = traj.s(:,1);
    % Generate feedback with LQR
    rr.generateFeedback(traj);

    % Evolve system dynamics and animate the robot arm
    % add zero velocity or perturb initial velocity
    %qin(3:4) = 0;
    % observe output
    %qact = rr.observeWithFeedbackErrorForm(traj,qin);
    qact = rr.observeWithFeedbackErrorForm(traj,qin,dmp{i});
    % add performance to trajectory
    traj.addPerformance([],qact,rr.COST,'ID + LQR');

    % Plot the controls and animate the robot arm
    rr.plot_outputs(traj);
    %rr.animateArm(qact(1:2,:),ref);

    % Start learning with ILC

    %ilc = wILC(traj,rr,'t');
    ilc = wILC(traj,rr,'dmp');

    for j = 1:num_trials
        % get next inputs
        dmp{i} = ilc.feedforward(traj,dmp{i},qact);
        %traj2 = ilc.feedforward(traj,[],qact);   
        % get the measurements
        qact = rr.observeWithFeedbackErrorForm(traj,qin,dmp{i});
        %qact = rr.observeWithFeedbackErrorForm(traj2,qin);
        traj.addPerformance([],qact,rr.COST,ilc);
        % Plot the controls and animate the robot arm
        %rr.animateArm(qact(1:2,:),ref);
    end

    % Plot the controls and animate the robot arm
    rr.plot_outputs(traj);
    rr.animateArm(qact(1:2,:),ref);
    
    % record weights of dmps
    p = [];
    for j = 1:length(dmp{i})
        p = [p; dmp{i}(j).Y0(:); dmp{i}(j).FORCE.w(:); dmp{i}(j).goal(:)];
    end
    wbar(:,i) = p; % record the parameters here
        
end

%% Sample a dmp from the convex hull and show generalization ability

% sample from the convex hull
vec = rand(4,1);
a = vec/sum(vec);
gnew = g * a;

[dmp,jNew] = adaptDMP(q0,gnew,dmpOrigin,w_origin);
[~,ref] = rr.kinematics(jNew);

traj = rr.generateInputs(t,ref);
qin = traj.s(:,1);
% Generate feedback with LQR
rr.generateFeedback(traj);

qact = rr.observeWithFeedbackErrorForm(traj,qin,dmp);
% add performance to trajectory
traj.addPerformance([],qact,rr.COST,'ID + LQR');

wnew = wbar * a;
% deal to individual dmps
wnew = reshape(wnew,length(wnew)/length(dmp),length(dmp));

for i = 1:length(dmp)

    % set/update goal and initial state
    dmp(i).setInitState(wnew(1:2,i));
    % initial states of DMPs
    dmp(i).setGoal(wnew(end,i));
    % set weights of dmp
    dmp(i).FORCE.w = wnew(3:end-1,i);
    [~,Q] = dmp(i).evolve();
    q(i,:) = Q(1,:);
end

ilc = wILC(traj,rr,'dmp');
inp = q(:);
dim = SIM.dimx;
inp = inp(1:end-dim);

qact = rr.observeWithFeedbackErrorForm(traj,qin,dmp);
% add performance to trajectory
traj.addPerformance(inp,qact,rr.COST,ilc);

for j = 1:num_trials
    % get next inputs
    dmp = ilc.feedforward(traj,dmp,qact);
    %traj2 = ilc.feedforward(traj,[],qact);   
    % get the measurements
    qact = rr.observeWithFeedbackErrorForm(traj,qin,dmp);
    %qact = rr.observeWithFeedbackErrorForm(traj2,qin);
    traj.addPerformance([],qact,rr.COST,ilc);
    % Plot the controls and animate the robot arm
    %rr.animateArm(qact(1:2,:),ref);
end

% Plot the controls and animate the robot arm
rr.plot_outputs(traj);
rr.animateArm(qact(1:2,:),ref);