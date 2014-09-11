%% Simulate trajectories for the planar RR arm

clc; clear; close all;

%% Define constants and parameters

% constants
g = 9.81;

% joint parameters
m1 = 1; %mass of first link, kg
m2 = 0.5; %mass of second link, kg
l1 = 0.50; %length of first link, m
l2 = 0.40; %length of second link, m
l_c1 = 0.25; %distance of first link's center of gravity to prev. joint, m
l_c2 = 0.20; %dist. of second link's c.o.g. to prev. joint, m
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
PAR.link1.motor.inertia = J_m1;
PAR.link2.motor.inertia = J_m2;
PAR.link1.motor.gear_ratio = r_1;
PAR.link2.motor.gear_ratio = r_2;

% form constraints
CON.link1.q.max = Inf;
CON.link1.q.min = -Inf;
CON.link1.u.max = Inf;
CON.link1.u.min = -Inf;
CON.link2.q.max = Inf;
CON.link2.q.min = -Inf;
CON.link2.u.max = Inf;
CON.link2.u.min = -Inf;

% cost structure
COST.Q = eye(2);

% initialize model
rr = RRplanar(PAR,CON,COST);

%% Generate a desired trajectory

h = 0.01;
y_des = 0.4:h:0.6;
x_des = 0.6 * ones(1,length(y_des));
t = h * 1:length(y_des);
X_des = [x_des;y_des];
q = RRplanarInverseKinematics(X_des,PAR);
% check for correctness
%[x1,x2] = RRplanarKinematics(q,PAR);
qd = diff(q')' / h; 
qdd = diff(qd')' / h; 

% keep velocity and acceleration vectors the same length as displacements
qd(:,end+1) = qd(:,end);
qdd(:,end+1) = qdd(:,end);
qdd(:,end+1) = qdd(:,end);

% get the desired inputs
ud = zeros(size(q));
for i = 1:size(q,2)
    ud(:,i) = RRplanarDynamics(q(:,i),qd(:,i),qdd(:,i),PAR);
end

%% Evolve system dynamics

% change parameters as a disturbance model
PAR.const.g = g; 
PAR.link1.mass = 1.0 * m1;
PAR.link2.mass = 1.0 * m2;
PAR.link1.length = 1.0 * l1;
PAR.link2.length = 1.0 * l2;
PAR.link1.centre.dist = l_c1;
PAR.link2.centre.dist = l_c2;
PAR.link1.inertia = 1.0 * I1;
PAR.link2.inertia = 1.0 * I2;
PAR.link1.motor.inertia = 1.0 * J_m1;
PAR.link2.motor.inertia = 1.0 * J_m2;
PAR.link1.motor.gear_ratio = 1.0 * r_1;
PAR.link2.motor.gear_ratio = 1.0 * r_2;

% TODO: add a nonzero friction matrix B

Q = [q; qd];
Q_a = Q(:,1);
fun = @(q,u) RRplanarInverseDynamics(q,u,PAR);
for i = 1:size(q,2)-1
    Q_a(:,i+1) = evolveDynamics(h,Q_a(:,i),ud(:,i),fun);
end

% get the cartesian coordinates of the actual trajectory followed
q_a = Q_a(1:2,:);
[x1,x2] = RRplanarKinematics(q_a,PAR);

%% Plot the controls and animate the robot arm

figure(1);
subplot(1,2,1);
plot(t,ud(1,:),'LineWidth',2);
title('Control input to first joint');
xlabel('Time (s)');
ylabel('Scaled voltages u_1 = r_1 * K_m_1 / R_1 * V_1');
subplot(1,2,2);
plot(t,ud(2,:),'LineWidth',2);
title('Control input to second joint');
ylabel('Scaled voltages u_2 = r_2 * K_m_2 / R_2 * V_2');
xlabel('Time (s)');

animateRR(x1,x2,X_des);

%% Start learning with ILC