%% Test Putting with full actuated RR arm

%clc; clear; close all;

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
% measurement noise covariance
SIM.eps_m = 3e-5;
% integration method
SIM.int = 'Euler';
% trajectory in joint space?
SIM.jref = false;

% constants
g = 9.81;
% joint parameters
m1 = 1.0; %mass of first link, kg
m2 = 0.5; %mass of second link, kg
l1 = 0.50; %length of first link, m
l2 = 0.40; %length of second link, m
l_c1 = 0.25; %distance of first link's center of gravity to prev. joint, m
l_c2 = 0.20; %dist. of second link's c.oZ.g. to prev. joint, m
I1 = (1/12)*m1*l1^2; %assume thin rod moment of inertia around c.o.g.
I2 = (1/12)*m1*l2^2; %kg m^2
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
PAR.const.g = g;
PAR.link1.mass = m1;
PAR.link2.mass = m2;
PAR.link1.length = l1;
PAR.link2.length = l2;
PAR.link1.centre.dist = l_c1;
PAR.link2.centre.dist = l_c2;
PAR.link1.inertia = I1;
PAR.link2.inertia = I2;
PAR.link1.motor.inertia = J_m1 + 0.05 * randn;
PAR.link2.motor.inertia = J_m2 + 0.05 * randn;
PAR.link1.motor.gear_ratio = r_1;
PAR.link2.motor.gear_ratio = r_2;

% TODO: we should observe joint angles only
PAR.C = eye(SIM.dimy,SIM.dimx);

% form constraints
CON = [];

% cost structure
% only penalize positions
COST.Q = 100*diag([1,1,0,0]);
COST.R = 1 * eye(SIM.dimu);

%% Calculate desired final velocity of the arm endeffector

% distance to the hole
dist = 0.45;
% hole radius
hole_rad = 0.05;
% total length of ball trajectory
l = dist + hole_rad;
% kinetic friction constant
mu_k = 0.6;
% mass of the ball
m_b = 1;
% gravity is already defined
% initial desired velocity of the ball
vin = sqrt(2*mu_k*g*l);
% calculate vf using conservation of momentum
vf = m_b * vin / (m1 + m2);

%% Generate inputs for a desired trajectory

h = SIM.h;
tin = 0; tfin = 1;
t = tin:h:tfin;
amp = pi/6;
sigma = 0.5;
mean = (tfin + tin)/2;
% trajectory consists of a compact gaussian
theta = amp * exp(-0.5*((t-mean).^2)/(sigma^2));
theta = theta - theta(1);
rad = l1 + l2;
y_des = -rad * sin(theta);
x_des = rad * cos(theta);
ref = [x_des; y_des]; % displacement profile 
% modify ref slightly to avoid singularity in inverse kinematics 
ref = ref - 1e-6;
% modify time profile to meet desired final velocity 
% in y_des due to 90 deg rotation
scale = vf/((y_des(end)-y_des(end-1))/h); 
h = h/scale;
t = t/scale;

% initialize RR model
SIM.h = h;
rr = RR(PAR,CON,COST,SIM);

% create trajectory in cartesian space
%traj = rr.generateInputs(t,ref);
bfs = 100;
% trajectory generated in joint space
[traj,dmp] = rr.generateInputsWithDMP(t,bfs,ref); 

%% Evolve system dynamics and animate the robot arm

q0 = traj.s(:,1);
% add nonzero velocity
%q0(3:4) = q0(3:4) + 0.1*rand(2,1);
% observe output
qact = rr.observe(t,q0,traj.unom);
% add performance to trajectory
traj.addPerformance(traj.unom,qact,rr.COST,'Computed Torque - IDM');

% Plot the controls and animate the robot arm
rr.plot_inputs(traj);
rr.plot_outputs(traj);
%rr.animateArm(qact(1:2,:),ref);

%% Start learning with ILC

num_trials = 10;

ilc = mILC(rr,traj);
%ilc = aILC(rr,traj);
if exist('ilcClass','var')
    if strcmp(ilcClass,'bILC')
        ilc = eval([ilcClass,'(traj,rr)']);
    else
        ilc = eval([ilcClass,'(rr,traj)']);
    end
end

for i = 1:num_trials
    % get next inputs
    u = ilc.feedforward(traj,qact);
    % evolve complete system
    qact = rr.observe(t,q0,u);
    % add performance to trajectory
    traj.addPerformance(u,qact,rr.COST,ilc);
    % Plot the controls and animate the robot arm
    %rr.animateArm(qact(1:2,:),ref);
end

% Plot the controls and animate the robot arm
rr.plot_inputs(traj);
rr.plot_outputs(traj);
%rr.animateArm(qact(1:2,:),ref);

%% ILC interaction with DMP

% num_trials = 10;
% ilc = mILC(rr,traj);
% uadd = zeros(size(traj.unom));
% e0diff = zeros(size(q0));
% q0new = q0;
% 
% for i = 1:num_trials
%     % change initial condition slightly
%     q0last = q0new;
%     q0new = q0 + 0.0 * randn(length(q0),1);
%     trajModified = rr.generateInputsWithDMP(t,bfs,ref,q0new);
%     trajNew = Trajectory(traj.t,trajModified.s,trajModified.unom,traj.K);
%     % adapt the dmps accordingly
%     % get the next inputs normally as in standard ILC
%     e0diff =  q0last - q0new;
%     uadd = ilc.feedforwardDMP(trajNew,qact,uadd,e0diff);
%     % adjust for the IDM change
%     trajNew.unom = uadd + trajNew.unom;
%     qact = rr.evolve(t,q0new,trajNew.unom);
%     trajNew.addPerformance(trajNew.unom,qact,rr.COST,ilc);
%     % Plot the controls and animate the robot arm
%     %rr.animateArm(qact(1:2,:),ref);
% end
% 
% % Plot the controls and animate the robot arm
% rr.plot_inputs(trajNew);
% rr.plot_outputs(trajNew);
% rr.animateArm(qact(1:2,:),ref);

%% Modify the animation
% 
% y = rr.kinematics(qact);
% shift = [0; 0.5];
% R = [0 1; -1 0];
% width = 0.025;
% % initial position of the ball
% ball = R * (y(1:2,end) + [0;width]);
% h1 = scatter(ball(1,1)+shift(1),ball(2,1)+shift(2),100,[.6 .6 .6],'filled','LineWidth',4);
% % % trajectory of the ball
% % balltrj_x = [y(1,end)+width y(1,end)+width];
% % balltrj_y = [y(2,end)+width/2 y(2,end) + l];
% % h2 = line(balltrj_x, balltrj_y, 'LineStyle', '-.', 'color', [.4 .4 .4],'LineWidth',1);
% % position the hole
% hole = [y(1,end)+width/2 y(1,end)+width/2; y(2,end)+l-hole_rad y(2,end)+l+hole_rad];
% hole = R * hole;
% h3 = line(hole(1,:)+shift(1), hole(2,:)+shift(2), 'LineStyle', '-', 'color', [0 0 0],'LineWidth',4);
% % print as grayscale eps 
% % print(gcf,'-depsc','putting.eps');