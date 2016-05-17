%% Initialize the RRR robot
% Simulation Values 
% system is continous
SIM.discrete = false;
% learn in cartesian space
SIM.cartesian = false;
% dimension of the x vector
SIM.dimx = 6;
% dimension of the output y
SIM.dimy = 3;
% dimension of the control input
SIM.dimu = 3;
% time step h 
SIM.h = 0.01;
% measurement noise covariance
SIM.eps_m = 3e-10;
% integration method
SIM.int = 'Euler';
% trajectory in joint space?
SIM.jref = false;

% constants
g = 9.81;
% joint parameters
m1 = 1; %mass of first link, kg
m2 = 0.5; %mass of second link, kg
m3 = 0.5; 
l1 = 0.5; %length of first link, m
l2 = 0.5; %length of second link, m
l3 = 0.5;
l_c1 = 0.25; %distance of first link's center of gravity to prev. joint, m
l_c2 = 0.20; %dist. of second link's c.oZ.g. to prev. joint, m
l_c3 = 0.5;
I1 = (1/12)*m1*l1^2; %assume thin rod moment of inertia around c.o.g.
I2 = (1/12)*m2*l2^2; %kg m^2
I3 = (1/12)*m3*l3^2;
% motor parameters
J_a1 = 0.100; % actuator inertia of link 1
J_g1 = 0.050; % gear inertia of link1
J_m1 = J_a1 + J_g1;
J_a2 = 0.080; % actuator inertia of link 2
J_g2 = 0.040; % gear inertia of link2
J_m2 = J_a2 + J_g2;
J_a3 = 0.12;
J_g3 = 0.06;
J_m3 = J_a3 + J_g3;
% gear ratio typically ranges from 20 to 200 or more - comment from book
r_1 = 20;
r_2 = 20;
r_3 = 40;
 
% pass it as a parameter structure
PAR.const.g = g;
PAR.link1.mass = m1;
PAR.link2.mass = m2;
PAR.link3.mass = m3;
PAR.link1.length = l1;
PAR.link2.length = l2;
PAR.link3.length = l3;
PAR.link1.centre.dist = l_c1;
PAR.link2.centre.dist = l_c2;
PAR.link3.centre.dist = l_c3;
PAR.link1.inertia = I1;
PAR.link2.inertia = I2;
PAR.link3.inertia = I3;
PAR.link1.motor.inertia = J_m1;
PAR.link2.motor.inertia = J_m2;
PAR.link3.motor.inertia = J_m3;
PAR.link1.motor.gear_ratio = r_1;
PAR.link2.motor.gear_ratio = r_2;
PAR.link3.motor.gear_ratio = r_3;

PAR.C = eye(SIM.dimy,SIM.dimx);

% form constraints
CON.link1.u.min = -Inf;
CON.link1.u.max = Inf;
CON.link2.u.min = -Inf;
CON.link2.u.max = Inf;
CON.link3.u.min = -Inf;
CON.link3.u.max = Inf;
CON.link1.udot.min = -Inf;
CON.link1.udot.max = Inf;
CON.link2.udot.min = -Inf;
CON.link2.udot.max = Inf;
CON.link3.udot.min = -Inf;
CON.link3.udot.max = Inf;

% cost structure
% only penalize positions
COST.Q = 100*diag([1,1,1,0,0,0]);
COST.R = 1 * eye(SIM.dimu);

% initialize model
rrr = RRR(PAR,CON,COST,SIM);

% initialize the arm with zero velocity on the right hand side
q0 = [0.3;-0.3;-0.4];
qd0 = zeros(3,1);
Q0 = [q0;qd0];